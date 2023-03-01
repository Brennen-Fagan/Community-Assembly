# Copy of Viking_HandleDiversity_HelperFunctions.R with
#  Viking_HandleDiversity_HelperFunctionsBC.R added.
#  Viking_HandleOutput_Invadability2Burnout_HiDisp.R
#  Viking_HandleOutput_DiversityBC.R

# Helper functions for Viking_HandleDiversity_ParametersAndPlots3.R (v3+)


# Functions: ###################################################################
## Functions from HandleOutput files.

# Recycling from Viking_HandleOutput_Diversity.R.
thinAndCalculateInvadabilities <- function(loaded, dyn, dis) {
  # We can't handle all of the data that we are going to be looking at;
  # a small sample had ~180k rows for ~5.3k events = ~34 rows per event.
  # To reduce it, we will divide time up so that there are about the
  # preferred number of rows per event.
  bythin <- floor((nrow(loaded$Abundance)
                   / nrow(loaded$Events))
                  / preferred_rows_per_event)

  loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                           to = nrow(loaded$Abundance),
                                           by = bythin), ]

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- loaded$Abundance[, -1] <
    loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
  loaded$Abundance[, -1][toEliminate] <- 0
  loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divide_time_by

  # Identify and grab only the last pre-burn-out time:
  target <- which.max(loaded$Abundance[, 1] > burnout / divide_time_by) - 1

  if (target == 0) {
    stop("Beginning of burn-out not found.")
  }

  # # Convert to Binary, since we are only interested in richness.
  # loaded$Abundance[, -1] <- loaded$Abundance[, -1] > 0

  # loaded$Abundance is now prepared.
  # Now, we calculate richnesses.
  print("prep")

  ### Invadability: #####################################################
  invadability <- CalculateLocalInvadables_KnockOn(
        Abundance = loaded$Abundance[target, -1],
        PerCapitaDynamics = dyn,
        Environments = loaded$NumEnvironments,
        ArrivalDensity = loaded$Parameters$ArrivalDensity,
        DispersalMatrix = dis,
        TimeScale = loaded$ReactionTime
      )

  binary <- rep(FALSE, length(loaded$Abundance[target, -1]))
  binary[invadability$candidates] <- invadability$invadable

  invadabilityMat <- (
    Matrix::drop0(Matrix::Matrix(binary, nrow = 1, ncol = length(binary)))
  )

  # candidateRegional = candidates who are not present on any patch
  # invadableRegional = above candidates, but who can invade any patch.
  # Matrix(Invadability$Invadabilities[[1]]$invadability[,], nrow = 100, ncol = 10)
  # Compare Matrix(Invadability$Invadabilities[[1]]$species, nrow = 100, ncol = 10)
  # Then candidates are then the all 1 rows of
  # Matrix(loaded$Abundance[target, -1] < EliminationThreshold, nrow = 100, ncol = 10)
  # and they are invadable if there are any 1's in there corresponding row of
  # Matrix(Invadability$Invadabilities[[1]]$invadability[,], nrow = 100, ncol = 10)

  theSpecies <- 1:((ncol(loaded$Abundance) - 1) / loaded$NumEnvironments)
  ### Return Invadabilities: ###############################################
  return(list(
    invadability = invadabilityMat, # Sparse
    time = loaded$Abundance[target, 1], # Not Sparse.
    species = rep(
      theSpecies,
      loaded$NumEnvironments
    ), # i.e. 1 2 3 1 2 3 1 2 3, Not Sparse.
    environment = rep(
      1:loaded$NumEnvironments,
      each = ((ncol(loaded$Abundance) - 1) / loaded$NumEnvironments)
    ), # i.e. 1 1 1 2 2 2 3 3 3, Not Sparse.
    speciesRegional = theSpecies,
    invadabilityRegional =
      apply(Matrix(loaded$Abundance[target, -1] <
                     loaded$Parameters$EliminationThreshold,
                   nrow = 100, ncol = 10), 1, all) *
      apply(Matrix(invadabilityMat[,], nrow = 100, ncol = 10), 1, any),
    effectAbundance = invadability$effectAbundance,
    effectRichness = invadability$effectRichness,
    effectEstablish = invadability$effectEstablish
  ))
}

# Recycling from SecondAttempt-Doc-Analysis2-Gallery.Rmd.
thinAndCalculateDiversities <- function(loaded, nspecies) {
  # We can't handle all of the data that we are going to be looking at;
  # a small sample had ~180k rows for ~5.3k events = ~34 rows per event.
  # To reduce it, we will divide time up so that there are about the
  # preferred number of rows per event.
  bythin <- floor((nrow(loaded$Abundance)
                   / nrow(loaded$Events))
                  / preferred_rows_per_event)

  loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                           to = nrow(loaded$Abundance),
                                           by = bythin), ]

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- loaded$Abundance[, -1] <
    loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
  loaded$Abundance[, -1][toEliminate] <- 0
  loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divide_time_by

  # Convert to Binary, since we are only interested in richness.
  loaded$Abundance[, -1] <- loaded$Abundance[, -1] > 0

  # loaded$Abundance is now prepared.
  # Now, we calculate richnesses.
  print("prep")

  ### Alpha Diversity: ##################################################
  diversity_alpha <- lapply(
    1:loaded$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      env_basal <- env[, 1:nspecies[1]]
      env_consumer <- env[, nspecies[1] + 1:nspecies[2]]
      richness <- rowSums(env)
      richness_basal <- rowSums(env_basal)
      richness_consumer <- rowSums(env_consumer)
      species <- apply(
        env, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_basal <- apply(
        env_basal, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_consumer <- apply(
        env_consumer, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0) + nspecies[1]) # + basals doesn't change it, but
          # results in a more intuitive numbering system of species.
        }
      )
      data.frame(Time = time,
                 Richness = richness,
                 Richness_Basal = richness_basal,
                 Richness_Consumer = richness_consumer,
                 Species = species,
                 Species_Basal = species_basal,
                 Species_Consumer = species_consumer,
                 Environment = i,
                 stringsAsFactors = FALSE)
    },
    abund = loaded$Abundance,
    numSpecies = sum(nspecies)
  )

  diversity_alpha <- dplyr::bind_rows(diversity_alpha)

  print("alpha")
  ### Gamma Diversity: ##################################################
  diversity_gamma <- diversity_alpha %>% dplyr::group_by(
    Time
  ) %>% dplyr::summarise(
    Mean = mean(Richness),
    Mean_Basal = mean(Richness_Basal),
    Mean_Consumer = mean(Richness_Consumer),
    Var = var(Richness),
    Var_Basal = var(Richness_Basal),
    Var_Consumer = var(Richness_Consumer),

    SpeciesTotal = toString(sort(unique(unlist(strsplit(paste(
      Species, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Basal = toString(sort(unique(unlist(strsplit(paste(
      Species_Basal, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Consumer = toString(sort(unique(unlist(strsplit(paste(
      Species_Consumer, collapse = ", "), split = ", ", fixed = TRUE))))),

    Gamma = unlist(lapply(
      strsplit(
        SpeciesTotal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Gamma_Basal = unlist(lapply(
      strsplit(
        SpeciesTotal_Basal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Gamma_Consumer = unlist(lapply(
      strsplit(
        SpeciesTotal_Consumer, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    ))
  ) %>% tidyr::pivot_longer(
    cols = c(Mean, Var, Gamma),
    names_to = "Aggregation",
    values_to = "Richness"
  ) %>% dplyr::mutate(
    Basals = case_when(
      Aggregation == "Mean" ~ Mean_Basal,
      Aggregation == "Var" ~ Var_Basal,
      Aggregation == "Gamma" ~ as.numeric(Gamma_Basal),
      TRUE ~ -1
    ),
    Consumers = case_when(
      Aggregation == "Mean" ~ Mean_Consumer,
      Aggregation == "Var" ~ Var_Consumer,
      Aggregation == "Gamma" ~ as.numeric(Gamma_Consumer),
      TRUE ~ -1
    )
  ) %>% dplyr::select(
    -Mean_Basal, -Var_Basal, -Gamma_Basal,
    -Mean_Consumer, -Var_Consumer, -Gamma_Consumer,
    -SpeciesTotal, -SpeciesTotal_Basal, -SpeciesTotal_Consumer
  )

  print("gamma")
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
        Jaccard = dists
      )

      return(dataf)
    },
    envs = loaded$NumEnvironments
  )
  print("beta")

  ### Return Diversities: ###############################################
  return(list(
    alpha = diversity_alpha %>% dplyr::select(
      -Species_Basal, -Species_Consumer
    ),
    beta = diversity_beta,
    gamma = diversity_gamma
  ))
}

## Attributes: ################################################################
extractAttributes <- function(Diversity, idNums) {
  data.frame(
    Pool = toString(Diversity$PoolMod),
    Noise = toString(Diversity$NoiseMod),
    Neutral = toString(Diversity$NeutralMod),
    Space = toString(Diversity$SpaceMod),
    Set = idNums[1],
    CaseNumber = idNums[2],
    History = if(is.na(idNums[3])) {
      # All Histories
      1:10
    } else idNums[3],
    Part = if(is.na(idNums[4])) {
      # All Parts
      1
    } else idNums[4],
    stringsAsFactors = FALSE
  )
}

assignAttributes <- function(i, d, a) {
  d %>% dplyr::bind_rows(
  ) %>% dplyr::mutate(
    Set = a$Set[i],
    Number = a$CaseNumber[i],
    History = a$History[i],
    Pool = a$Pool[i],
    Noise = a$Noise[i],
    Neutral = a$Neutral[i],
    Space = a$Space[i]
  )
}

## Diversity Metrics: ##########################################################
# Primary distinctions are what metric is extracted in function(i, d, a),
# what columns might be dropped, what variables are used to group, and flooring.
extractAlphas <- function(Diversity, Attributes) {
  dplyr::bind_rows(
    lapply(
      seq_along(Diversity$Diversities),
      function(i, d, a) assignAttributes(i, d[[i]]$alpha, a),
      d = Diversity$Diversities,
      a = Attributes
    )
  ) %>% dplyr::select(-Species) %>% dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) %>% dplyr::group_by(
    Time, Environment, Set, Number, History, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    Richness = floor(median(Richness)),
    Richness_Basal = floor(median(Richness_Basal)),
    Richness_Consumer = floor(median(Richness_Consumer)),
    .groups = "drop"
  )
}

extractBetas <- function(Diversity, Attributes) {
  dplyr::bind_rows(
    lapply(
      seq_along(Diversity$Diversities),
      function(i, d, a) assignAttributes(i, d[[i]]$beta, a),
      d = Diversity$Diversities,
      a = Attributes
    )
  ) %>% dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) %>% dplyr::group_by(
    Time, Env1, Env2, Set, Number, History, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    Jaccard = median(Jaccard),
    .groups = "drop"
  )
}

extractGammas <- function(Diversity, Attributes) {
  dplyr::bind_rows(lapply(
    seq_along(Diversity$Diversities),
    function(i, d, a) assignAttributes(i, d[[i]]$gamma, a),
    d = Diversity$Diversities,
    a = Attributes)) %>% dplyr::filter(
      Aggregation == "Gamma"
    ) %>% dplyr::mutate(
      Time = floor(Time * time_grouping_size)/time_grouping_size
    ) %>% dplyr::group_by(
      Time, Set, Number, History, Pool, Noise, Neutral, Space
    ) %>% summarise(
      Richness = floor(median(Richness)),
      Richness_Basal = floor(median(Basals)),
      Richness_Consumer = floor(median(Consumers)),
      .groups = "drop"
    )
}

## Load Function: ##############################################################
loadDiversity <- function(file, verbose = TRUE) {
  # Retrieve the trailing id numbers from before the file extension.
  # 1st: Set, 2nd: CaseNumber, 3rd: History, 4th: Part
  # Note, if last two are not present, all histories are bundled together.
  # If last two are present, each file has a single part of a single history.
  idNums <- suppressWarnings(na.omit(
    as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                             split = "-", fixed = TRUE)[[1]],
                    n = 4))),
    classes = "simpleWarning")

  if (verbose) {
    print(file)
    print(Sys.time())
  }

  tryCatch({
    load(file)

    if (!exists("Diversity")) {stop("Diversity does not exist.")}

    Attributes <- extractAttributes(Diversity, idNums)
    DiversitiesAlpha <- extractAlphas(Diversity, Attributes)
    DiversitiesBeta <- extractBetas(Diversity, Attributes)
    DiversitiesGamma <- extractGammas(Diversity, Attributes)

    return(list(
      Attr = Attributes,
      Alpha = DiversitiesAlpha,
      Beta = DiversitiesBeta,
      Gamma = DiversitiesGamma
    ))
  }, error = function(e) {
    print(paste(paste(idNums, collapse = " "), e, sep = ": "))
    return(NULL)
  })
}

## Average Function: ###########################################################
createAvg <- function(d, target, ps = ribbonQuantiles, cE = centralEst) {
  d %>% dplyr::mutate(
    Time = floor(Time * time_averaging_size)/time_averaging_size
  ) %>%  dplyr::group_by(
    Time, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    Low = unlist(dplyr::across(dplyr::any_of(target),
                               .fns = ~ quantile(.x, p = ps[1], na.rm = TRUE))),
    High = unlist(dplyr::across(dplyr::any_of(target),
                                .fns = ~ quantile(.x, p = ps[2], na.rm = TRUE))),
    Central = unlist(dplyr::across(dplyr::any_of(target),
                                   .fns = ~ cE(.x))),
    .groups = "drop"
  )
}

## Plotting Functions: #########################################################
# Format is [Y-axis]Over[X-axis].
StatisticOverTime <- function(
  i, d, da, levels,
  spacelabs, neutrallabs,
  ytarget, ylabel, title,
  burnin = time_burnin,
  burnout = time_burnout,
  units = time_units
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]

  # Attach labels.
  d <- d %>% dplyr::filter(
    pbs == Pool,
    inm == Noise
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  ) %>% tidyr::unite(
    col = "group", !Time & !dplyr::all_of(ytarget), remove = FALSE
  )

  da <- da %>% dplyr::filter(
    pbs == Pool,
    inm == Noise
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  ) %>% tidyr::unite(
    col = "group", !Time & !Low & !High & !Central, remove = FALSE
  )

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Time,
                 group = group)
  ) + ggplot2::geom_line(
    mapping = ggplot2::aes_string(
      y = ytarget
    ),
    alpha = 0.05
  ) + ggplot2::geom_line(
    data = da,
    mapping = ggplot2::aes(
      y = Central
    ),
    color = "red", size = 1.5
  ) + ggplot2::geom_ribbon(
    data = da,
    mapping = ggplot2::aes(
      ymin = Low,
      ymax = High,
      x = Time
    ),
    alpha = 0.4,
    fill = "blue",
    inherit.aes = FALSE
  ) + ggplot2::geom_vline(
    xintercept = c(burnin, burnout)
  ) + ggplot2::theme_bw(
  ) + ggplot2::facet_grid(
    Immigration + Extinction ~ Space + Dispersal,
    labeller = ggplot2::label_parsed
  ) + ggplot2::labs(
    title = title,
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm),
    ylab = ylabel
  ) + ggplot2::xlab(
    bquote("Time, " * 10^.(log10(time_units)) * " Units")
  )

  return(obj)
}

StatisticOverSpace <- function(
  i, d, levels,
  spacelabs,
  # ytarget, SET TO Richness
  ylabel, title,
  burnin = time_burnin,
  burnout = time_burnout
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]

  # Attach labels.
  d <- d %>% dplyr::filter(
    pbs == Pool,
    inm == Noise,
    Time > burnin, Time < burnout
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = "")
  ) %>% tidyr::unite(
    col = "group", Neutral, Space, remove = FALSE
  )

  xlabs <- d %>% dplyr::select(
    Space, Dispersal
  ) %>% dplyr::distinct() %>% dplyr::arrange(Space) %>% dplyr::mutate(
    Space = as.character(Space),
    Label = paste(Space, Dispersal, sep = "\n")
  )

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = as.character(Space), y = Richness,
                 fill = Neutral,
                 group = group)
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.5), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.5),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
    title = title,
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm),
    y = ylabel
  ) + ggplot2::scale_x_discrete(
    name = "Inter-community distance and unit proportion dispersing",
    breaks = xlabs$Space,
    labels = xlabs$Label
  ) + ggplot2::scale_fill_discrete(
    name = "(Imm., Ext.)\nModifiers"
  )

  return(obj)
}

AvgGammaAlphaOverTime <- function(DAG) {
  DAG %>% dplyr::group_by(
    Time, Set, Number, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    `Richness, Alpha` = mean(`Richness, Alpha`),
    `Richness_Basal, Alpha` = mean(`Richness_Basal, Alpha`),
    `Richness_Consumer, Alpha` = mean(`Richness_Consumer, Alpha`),
    `Richness, Gamma` = mean(`Richness, Gamma`),
    `Richness_Basal, Gamma` = mean(`Richness_Basal, Gamma`),
    `Richness_Consumer, Gamma` = mean(`Richness_Consumer, Gamma`),
    .groups = "drop"
  ) %>% tidyr::pivot_longer(
    cols = c("Richness, Alpha", "Richness, Gamma",
             "Richness_Basal, Alpha", "Richness_Consumer, Alpha",
             "Richness_Basal, Gamma", "Richness_Consumer, Gamma"),
    names_to = "Measurement", values_to = "Value"
  ) %>% ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = Time, y= Value, 
      group = Measurement, color = Measurement
    )
  ) + ggplot2::geom_line(
  )
}

GammaOverAlpha <- function(
  i, d, levels,
  spacelabs, neutrallabs,
  burnin = time_burnin,
  burnout = time_burnout,
  # For each parameter set, there are 10 Pools, and in each...
  # Pool (for which there are 10 Histories),
  # History (for which there are 10 Environments), Environment
  NumberSelect = NULL ,#1,
  HistorySelect = NULL, #1,
  EnvirSelect = NULL #1
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]

  d <- d %>% dplyr::filter(
    Pool == pbs,
    Noise == inm
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  )

  if (!is.null(NumberSelect) &&
      !is.null(HistorySelect) &&
      !is.null(EnvirSelect)) {
    number <- d %>% dplyr::select(
      Space, Neutral, Number, History, Environment
    ) %>% dplyr::distinct() %>% dplyr::arrange(
      Number
    ) %>% dplyr::group_by(
      Space, Neutral
    ) %>% dplyr::filter(
      Number %in% (min(Number) - 1 + NumberSelect)
    ) %>% dplyr::group_by(
      Space, Neutral, Number
    ) %>% dplyr::filter(
      History %in% (min(History) - 1 + HistorySelect)
    ) %>% dplyr::group_by(
      Space, Neutral, Number, History
    ) %>% dplyr::filter(
      Environment %in% (min(Environment) - 1 + EnvirSelect)
    )

    dTraj <- d %>% dplyr::semi_join(number, by = c(
      "Space", "Neutral", "Number", "History", "Environment"
    ))

    dTrajStop <- dTraj %>% dplyr::group_by(
      Number, History, Environment
    ) %>% dplyr::slice_max(Time)

  } else if (!is.null(NumberSelect) ||
             !is.null(HistorySelect) ||
             !is.null(EnvirSelect)) {
    warning("Not all of Selects are null or non-null. Skipping trajectory.")
  }

  dga <- d %>% dplyr::filter(
    Time > burnin, Time < burnout
  )

  obj <- ggplot2::ggplot(
    dga,
    ggplot2::aes(x = `Richness, Alpha`,
                 y = `Richness, Gamma`)
  ) + ggplot2::geom_bin2d(
    breaks = seq(-0.5, max(dga$`Richness, Alpha`,
                           dga$`Richness, Gamma`) + 1, 1)
  ) + ggplot2::geom_point(
    data = dga %>% group_by(
      Immigration, Extinction, Space, Dispersal
    ) %>% summarise(
      `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
      `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
      .groups = "drop"
    ),
    shape = 4, size = 3, color = "red", alpha = 0.8
  ) + ggplot2::scale_fill_viridis_c(
    direction = -1, trans = "log10"
  ) + ggplot2::facet_grid(
    Immigration + Extinction ~ Space + Dispersal,
    labeller = ggplot2::label_parsed
  ) + ggplot2::theme_bw(
  ) + ggplot2::geom_abline(
    slope = 1, intercept = 0
  ) + ggplot2::labs(
    title = if (exists("dTraj")) {
      "Sample trajectories in diversity counts"
    } else {
      "Gamma by Alpha"
    },
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm)
  )

  if (exists("dTraj")) {
    obj <- obj + ggplot2::geom_path(
      data = dTraj,
      mapping = ggplot2::aes(
        group = interaction(Environment, Number, History),
        color = Time
      ), alpha = 0.3
    ) + ggplot2::geom_point(
      data = dTrajStop,
      shape = 1, size = 3, color = "yellow", alpha = 0.8
    ) + ggplot2::scale_color_viridis_c(
      option = "B"
    )
  }

  return(obj)
}

