library(RMTRCode2)
library(dplyr)
library(ggplot2)

if(exists(".Random.seed")){
  oldseed <- .Random.seed
}
# > runif(1) * 1E8
# [1] 4966376
set.seed(4966376)

ranSeeds <- runif(5) * 1E8

numBasal <- 3; numConsum <- 2; numEnviron <- 3
EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 30 # Number of steps to reach next event to smooth.

divide_time_by <- 1E3
preferred_rows_per_event <- 1.5

egPool <- RMTRCode2::LawMorton1996_species(numBasal, numConsum,
                                           seed = ranSeeds[1])

egInteractions <- CreateEnvironmentInteractions(
  egPool, numEnviron, LawMorton1996_CommunityMat,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 10),
  EnvironmentSeeds = ranSeeds[2:4]
)

egMatrix <- Matrix::bdiag(egInteractions$Mats)

eigs <- eigen(egMatrix, only.values = TRUE)$values
rate <- max(abs(eigs))

egEvents <- CreateAssemblySequence(
  (numBasal + numConsum), numEnviron,
  # Coupon Collector's Problem, Limit Theorem Probability Bound
  ArrivalEvents = ceiling((numBasal + numConsum) * numEnviron * (
    log((numBasal + numConsum) * numEnviron) + 5
  )),
  ArrivalRate = rate,
  ArrivalFUN = RMTRCode2::ArrivalFUN_Example2,
  ExtinctFUN = RMTRCode2::ExtinctFUN_Example2,
  HistorySeed = ranSeeds[5])

egDynamics <- PerCapitaDynamics_Type1(egPool$ReproductionRate, egMatrix,
                                      NumEnvironments = numEnviron)

egEventFunction <- EliminationAndNeutralEvents(
  EventsAndSeed = egEvents, Species = nrow(egPool),
  PerCapitaDynamics = egDynamics,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity
)

PerIslandDistances <- c(1e0, 1e5, Inf) #10^7
Environments <- numEnviron # Synonym, but borrowed from different code blocks.

resultsSpace <- lapply(PerIslandDistances, function(PerIslandDistance) {
  # Ring Dynamics
  DistanceMatrix <- Matrix::bandSparse(
    Environments, k = c(-1, 1),
    diagonals = list(rep(PerIslandDistance, Environments - 1),
                     rep(PerIslandDistance, Environments - 1))
  )
  DistanceMatrix[Environments, 1] <- PerIslandDistance
  DistanceMatrix[1, Environments] <- PerIslandDistance

  DispersalMatrices <- RMTRCode2::CreateDispersalMatrix(
    EnvironmentDistances = DistanceMatrix,
    SpeciesSpeeds = rep(1, nrow(egPool))
  )

  results <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
    Pool = egPool, NumEnvironments = numEnviron,
    InteractionMatrices = egInteractions,
    Events = egEvents,
    PerCapitaDynamics = egDynamics,
    DispersalMatrix = DispersalMatrices,
    EliminationThreshold = EliminationThreshold,
    ArrivalDensity = ArrivalDensity,
    MaximumTimeStep = MaximumTimeStep,
    BetweenEventSteps = BetweenEventSteps,
    Verbose = FALSE
  )

  # We can't handle all of the data that we are going to be looking at;
  # a small sample had ~180k rows for ~5.3k events = ~34 rows per event.
  # To reduce it, we will divide time up so that there are about the
  # preferred number of rows per event.

  bythin <- floor((nrow(results$Abundance)
                   / nrow(results$Events))
                  / preferred_rows_per_event)

  results$Abundance <- results$Abundance[seq(from = 1,
                                             to = nrow(results$Abundance),
                                             by = bythin), ]

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- results$Abundance[, -1] <
    results$Parameters$EliminationThreshold & results$Abundance[, -1] > 0
  results$Abundance[, -1][toEliminate] <- 0
  results$Abundance[, 1] <- results$Abundance[, 1] / divide_time_by

  results$Space <- PerIslandDistance

  return(results)
})

diversitiesSpace <- lapply(resultsSpace, function(loaded) {
  # Convert to Binary, since we are only interested in richness.
  loaded$Abundance[, -1] <- loaded$Abundance[, -1] > 0
  # loaded$Abundance is now prepared.
  # Now, we calculate richnesses.

  ### Alpha Diversity: ##################################################
  diversity_alpha <- lapply(
    1:loaded$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      richness <- rowSums(env)
      species <- apply(
        env, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      data.frame(Time = time,
                 Richness = richness,
                 Species = species,
                 Environment = i,
                 Space = loaded$Space,
                 stringsAsFactors = FALSE)
    },
    abund = loaded$Abundance,
    numSpecies = (ncol(loaded$Abundance) - 1) / loaded$NumEnvironments
  )

  diversity_alpha <- dplyr::bind_rows(diversity_alpha)

  ### Gamma Diversity: ##################################################
  diversity_gamma <- diversity_alpha %>% dplyr::group_by(
    Time
  ) %>% dplyr::summarise(
    Mean = mean(Richness),
    SpeciesTotal = toString(sort(unique(unlist(strsplit(paste(
      Species, collapse = ", "), split = ", ", fixed = TRUE))))),
    Gamma = unlist(lapply(
      strsplit(
        SpeciesTotal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Space = mean(Space),
    .groups = "drop"
  ) %>% tidyr::pivot_longer(
    cols = c(Mean, Gamma),
    names_to = "Aggregation",
    values_to = "Richness"
  ) #%>% dplyr::select(-SpeciesTotal)

  ### Beta Diversity (Jaccard, Space): ##################################
  diversity_beta <- apply(
    loaded$Abundance,
    MARGIN = 1, # Rows
    function(row, envs) {
      time <- row[1]
      # Vegan complains about rows with all 0's.
      # The warning is generic, so we cannot silence it specifically.
      dists <- suppressWarnings(vegan::vegdist(
        method = "jaccard",
        x = matrix(row[-1] > 0, nrow = envs, byrow = TRUE)
      ))

      dataf <- expand.grid(
        Env1 = 1:envs,
        Env2 = 1:envs
      ) %>% dplyr::filter(
        Env1 < Env2
      ) %>% dplyr::mutate(
        Time = time,
        Jaccard = dists,
        Space = loaded$Space
      )

      return(dataf)
    },
    envs = loaded$NumEnvironments
  )

  ### Return Diversities: ###############################################
  return(list(
    alpha = diversity_alpha,
    beta = diversity_beta,
    gamma = diversity_gamma
  ))
}
)

# If we are interested in between system Jaccards, we'll need to compute those
# separately.
jaccardsSpace <- lapply(diversitiesSpace, function(div) {
  div$gamma %>% dplyr::filter(Aggregation == "Gamma") %>% dplyr::transmute(
    Time = Time,
    S1 = grepl("1", SpeciesTotal, fixed = TRUE),
    S2 = grepl("2", SpeciesTotal, fixed = TRUE),
    S3 = grepl("3", SpeciesTotal, fixed = TRUE),
    S4 = grepl("4", SpeciesTotal, fixed = TRUE),
    S5 = grepl("5", SpeciesTotal, fixed = TRUE)
  )
})
jaccardsSpace <- lapply(jaccardsSpace[[1]]$Time, function(time) {
  dists <- vegan::vegdist(
    method = "jaccard", binary = TRUE,
    do.call(rbind,
            lapply(jaccardsSpace, function(js) {
              js %>% dplyr::filter(
                Time == time
              ) %>% dplyr::select(
                -Time
              ) %>% as.matrix(
                nrow = 1
              )
            }))
  )
  dataf <- expand.grid(
    Env1 = seq_along(jaccardsSpace),
    Env2 = seq_along(jaccardsSpace)
  ) %>% dplyr::filter(
    Env1 < Env2
  ) %>% dplyr::mutate(
    Time = time,
    Jaccard = as.numeric(dists)
  )

  return(dataf)
}) %>% dplyr::bind_rows()

# Formatting:
#  x = Time,
#  y = Value, which includes Alpha, Beta, Gamma, Island 1, Island 2, Island 3
#  group = Type, which includes Species Number or Environment (Combination).
#  Facet = Alpha, Beta, Gamma, Island 1, Island 2, Island 3

plottingData <- dplyr::bind_rows(
  lapply(resultsSpace, function(results) {
    data.frame(results$Abundance) %>% tidyr::pivot_longer(
      dplyr::starts_with("X"),
      names_to = "SpeciesEnvironment", values_to = "Value"
    ) %>% dplyr::mutate(
      SpeciesEnvironment = as.numeric(substring(SpeciesEnvironment, 2)),
      Species = ((SpeciesEnvironment - 1) %% nrow(egPool)) + 1,
      Environment = ((SpeciesEnvironment - 1) %/% nrow(egPool)) + 1
    ) %>% dplyr::transmute(
      Time = time,
      Value = log10(Value),
      Type = Species, # paste("Species", Species),
      Island = Environment,
      Facet = paste("Log10 Abundance,", "Space", results$Space)
      #, Island", Environment)
    ) %>% dplyr::filter(Island == 1)
  }),
  lapply(diversitiesSpace, function(diversities) {
    diversities$alpha %>% dplyr::group_by(
      Time
    ) %>% dplyr::summarise( # Over Environment and Species
      Time = Time,
      Value = mean(Richness),
      Type = Space, #Environment * 2 - 1, # paste("Island", Environment),
      Facet = "Mean Richness, Alpha",
      .groups = "drop"
    )
  }),

  # dplyr::bind_rows(diversities$beta) %>% dplyr::transmute(
  #   Time = Time,
  #   Value = as.numeric(Jaccard),
  #   Type = as.numeric(factor(paste(Env1, Env2))) * 2 - 1,
  #   Facet = "Jaccard Distance"
  # ),

  jaccardsSpace %>% dplyr::transmute(
    Time = Time,
    Value = Jaccard,
    Type = as.numeric(interaction(Env1, Env2)),
    Facet = "Jaccard Distance"
  ),

  lapply(diversitiesSpace, function(diversities) {
    diversities$gamma %>% dplyr::select(
      -SpeciesTotal
    ) %>% dplyr::filter(
      Aggregation == "Gamma"
    ) %>% dplyr::transmute(
      Time = Time,
      Value = Richness,
      Type = Space, # 1, # "Gamma",
      Facet = "Richness, Gamma"
    )
  })
)

plottingData <- plottingData %>% dplyr::mutate(
  dplyr::across(Facet, factor, levels = c(
    "Log10 Abundance, Space 1",
    "Log10 Abundance, Space 1e+05",
    "Log10 Abundance, Space Inf",
    "Mean Richness, Alpha",
    "Jaccard Distance",
    "Richness, Gamma"
  ))
)

print(ggplot2::ggplot(
  plottingData,
  ggplot2::aes(
    x = Time,
    y = Value,
    color = factor(Type)
  )
) + ggplot2::geom_line(
) + ggplot2::theme_bw(
) + ggplot2::facet_wrap(
  ~ Facet, scales = "free_y"
  #) + ggplot2::scale_color_viridis_d(
) + ggplot2::guides(
  color = "none"
) + ggplot2::labs(
  x = bquote(paste(Time, ",") ~ 10^.(log10(divide_time_by)) ~ units)
) + ggplot2::coord_cartesian(
  ylim = c(0, NA)
))

# From Server_HandleDiversity_ParametersAndPlots2

time_grouping_size <- 100
time_averaging_size <- 10

DiversitiesAlphaGamma <- dplyr::full_join(
  diversities$alpha %>% dplyr::select(-Species) %>% dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) %>% dplyr::group_by(
    Time, #History,
    Environment#, Case
  ) %>% summarise(
    Richness = floor(median(Richness)),
    .groups = "drop"
  ),
  diversities$gamma %>% dplyr::filter(
    Aggregation == "Gamma"
  ) %>% dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) %>% dplyr::group_by(
    Time#, History, Case
  ) %>% summarise(
    Richness = floor(median(Richness)),
    .groups = "drop"
  ),
  by = c("Time"#, "History", "Case",
         #"Pool", "Noise", "Neutral", "Space"
         ),
  suffix = c(", Alpha", ", Gamma")
)

with(list(dag = DiversitiesAlphaGamma ), {

  ggplot2::ggplot(
    dag %>% dplyr::filter(Time > 3),
    ggplot2::aes(x = `Richness, Alpha`,
                 y = `Richness, Gamma`)
  ) + ggplot2::geom_bin2d(
    breaks = seq(-0.5, max(dag$`Richness, Alpha`,
                           dag$`Richness, Gamma`) + 1, 1)
  ) + ggplot2::geom_point(
    data = dag %>% dplyr::filter(Time > 3) %>% group_by(
      #Neutral, Space
    ) %>% summarise(
      `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
      `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
      .groups = "drop"
    ),
    shape = 4, size = 3, color = "red", alpha = 0.8
  ) + ggplot2::scale_fill_viridis_c(
    direction = -1, trans = "log10"
    #) + ggplot2::facet_grid(
    #    Neutral ~ Space
  ) + ggplot2::theme_bw(
  ) + ggplot2::geom_path(
    data = dag %>% dplyr::filter(
      Environment == 2
    ),
    mapping = ggplot2::aes(
      color = Time
    )
  ) + ggplot2::scale_color_viridis_c(
    option = "B"
  ) + ggplot2::geom_abline(
    slope = 1, intercept = 0
  ) + ggplot2::labs(
    title = "Gamma by Alpha",
    #subtitle = paste0("Pool Modifier = ", pbs, ", ",
    #                  "Noise Modifier = ", inm)
  )
}) %>% print()

if(exists("oldseed"))
  set.seed(oldseed)
