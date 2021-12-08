# We are going to look for a few stable systems specifically and see how the
# known-independent runs are correlated when the pool is small (and should
# only have a single known steady-state with no transitions between
# steady-states).
# Since we only want stable systems, we can solve for the appropriate
# reproductive rates in order to maintain a system of the given sizes.
# We do search for a random seed for a system that will be best behaved.
library(RMTRCode2)

Species1 <- c(Basal = 1, Consumer = 1)
Species2 <- c(Basal = 1, Consumer = 2)
Species3 <- c(Basal = 3, Consumer = 0)
Species4 <- c(Basal = 10, Consumer = 0)
# Species5 <- c(Basal = 30, Consumer = 0)
Species6 <- c(Basal = 2, Consumer = 0)
Environments <- 10
EventsEach <- Environments * ceiling(10 * (log(10 + 5)))

LMParameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
LMLogBodySize <- c(-2, -1, -1, 0)

PerIslandDistance <- Inf
SpeciesSpeeds <- 1
TimeScaling <- 10
Space <- match.arg("Ring", c("None", "Ring", "Line", "Full"))

EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 30 # Number of steps to reach next event to smooth.

# 1, 9, 10 work but have low (1E-5) characteristic rates.
# 13, 16 bottom out at 1E-4
# 140, 375, 402 are good at 2E-4, but this seems to be about the best we can do.
PoolSeed <- 375
# > runif(1) * 1e8
# [1] 12032489
EnvironmentSeed <- 12032489
# > runif(1) * 1e8
# [1] 28665115
HistorySeed <- 28665115

Species <- list(
  Species1, Species2, Species3, Species4, Species6
)

Pools <- lapply(Species, function(Sp) {
  LawMorton1996_species(
    Basal = Sp[1],
    Consumer = Sp[2],
    Parameters = LMParameters,
    LogBodySize = LMLogBodySize,
    seed = PoolSeed
  )
})

InteractionMatrices <- lapply(
  Pools,
  function(Pool) {
    temp <- CreateEnvironmentInteractions(
      Pool = Pool, NumEnvironments = Environments,
      ComputeInteractionMatrix = LawMorton1996_CommunityMat,
      Parameters = LMParameters,
      EnvironmentSeeds = EnvironmentSeed
    )

    temp$Mats <- lapply(temp$Mats, function(m) m * TimeScaling)

    temp
  }
)

Pools <- lapply(
  seq_along(Pools),
  function(i, Pools, Mats) {
    if ("Consumer" %in% Pools[[i]]$Types) {
      Pools[[i]]$ReproductionRate <- -rowMeans(do.call(
        cbind,
        lapply(Mats[[i]]$Mats,
               function(m) m %*% rep(1, nrow(Pools[[i]])))
      ))

      stopifnot(
        all(unlist(
          lapply(Mats[[i]]$Mats,
                 function(m) - solve(m) %*% Pools[[i]]$ReproductionRate)) > 0)
        )

      Pools[[i]]
    } else {
      Pools[[i]]
    }
  },
  Pools = Pools, Mats = InteractionMatrices
)

CharacteristicRates <- lapply(
  InteractionMatrices,
  function(Mat) {
    max(unlist(lapply(
      Mat$Mats, function(m) abs(eigen(m, only.values = TRUE)$values))))
  }
)
print(CharacteristicRates)

Events <- lapply(
  seq_along(Species), function(
    i, spec, charrate, envs, eventnum, seed
  ) {
    CreateAssemblySequence(
      Species = sum(spec[[i]]),
      NumEnvironments = envs,
      ArrivalEvents = eventnum,
      ArrivalRate = charrate[[i]],
      ArrivalFUN = ArrivalFUN_Example,
      ExtinctEvents = eventnum,
      ExtinctRate = charrate[[i]],
      ExtinctFUN = ExtinctFUN_Example,
      HistorySeed = seed
    )
  },
  envs = Environments,
  seed = HistorySeed,
  eventnum = EventsEach,
  spec = Species,
  charrate = CharacteristicRates
)

IntMats <- lapply(InteractionMatrices, function(m) {
  Matrix::bdiag(m$Mats)
})

PerCapitaDynamics <- lapply(
  seq_along(Pools), function(i, pools, mats, envs) {
    PerCapitaDynamics_Type1(
      pools[[i]]$ReproductionRate, mats[[i]],
      NumEnvironments = envs
    )
  },
  pools = Pools,
  mats = IntMats,
  envs = Environments
)

if (Space == "None") {
  DistanceMatrix <- Matrix::sparseMatrix(
    i = Environments, j = Environments, x = 0)
}
if (Space == "Ring" || Space == "Line")
  DistanceMatrix <- Matrix::bandSparse(
    Environments, k = c(-1, 1),
    diagonals = list(rep(PerIslandDistance, Environments - 1),
                     rep(PerIslandDistance, Environments - 1))
  )
if (Space == "Ring") {
  DistanceMatrix[Environments, 1] <- PerIslandDistance
  DistanceMatrix[1, Environments] <- PerIslandDistance
}
if (Space == "Grid") {
  # Given matrix(1:4, nrow = 2), trying 1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4.
  # I.e. matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0), nrow = 4)
  # Use divisor closest to but <= square root for number of rows.
}
if (Space == "Full") {
  DistanceMatrix <- matrix(1, nrow = Environments, ncol = Environments)
  diag(DistanceMatrix) <- 0
}

DispersalMatrices <- lapply(
  Pools, function(p, dm, speed) {
    CreateDispersalMatrix(
      EnvironmentDistances = dm,
      SpeciesSpeeds = rep(speed, nrow(p))
    )
  },
  dm = DistanceMatrix, speed = SpeciesSpeeds
)

results <- lapply(
  seq_along(Pools), function(
    i, pools, envs, mats, events, dynams, dispmats,
    elm, arrdens, maxt, bett
  ) {
    print(Sys.time())
    MultipleNumericalAssembly_Dispersal(
      Pool = pools[[i]], NumEnvironments = envs,
      InteractionMatrices = mats[[i]],
      Events = events[[i]],
      PerCapitaDynamics = dynams[[i]],
      DispersalMatrix = dispmats[[i]],
      EliminationThreshold = elm,
      ArrivalDensity = arrdens,
      MaximumTimeStep = maxt,
      BetweenEventSteps = bett,
      Verbose = FALSE
    )
  },
  pools = Pools,
  envs = Environments,
  mats = InteractionMatrices,
  events = Events,
  dynams = PerCapitaDynamics,
  dispmats = DispersalMatrices,
  elm = EliminationThreshold,
  arrdens = ArrivalDensity,
  maxt = MaximumTimeStep,
  bett = BetweenEventSteps
)
print(Sys.time())

`%>%` <- magrittr::`%>%`

log0 <- function(x, base = exp(1)) {
  xneq0 <- x != 0 & !is.na(x)
  x[xneq0] <- log(x[xneq0], base = base)
  return(x)
}

Calculate_Diversity <- function(result) {
  Diversity <- lapply(
    1:result$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      richness <- rowSums(env != 0)
      abundSum <- rowSums(env)
      #NOTE: THIS CAN YIELD NAN'S (0/0).
      # THIS IS NOT NECESSARILY A PROBLEM.
      # IT MIGHT BE WORTH IT JUST TO USE 0 OR
      # TO CATCH IT EXPLICITLY AND REPLACE WITH NAN.
      entropy <- env / abundSum
      entropy <- - apply(
        entropy, MARGIN = 1,
        FUN = function(x) {
          sum(x * log0(x))
        })
      species <- apply(
        env, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      evenness <- entropy / log(richness)
      data.frame(Time = time,
                 Richness = richness,
                 Entropy = entropy,
                 Evenness = evenness,
                 Species = species,
                 Environment = i,
                 stringsAsFactors = FALSE)
    },
    abund = result$Abundance,
    numSpecies = (ncol(result$Abundance) - 1) / result$NumEnvironments
  )


  Diversity <- dplyr::bind_rows(Diversity)
  Diversity_alpha <- Diversity
  # Diversity_alpha <- Diversity_alpha %>% dplyr::mutate(
  #   Evenness = Entropy / log(Richness)
  # )

  # Modify to do the gamma bits right here.
  Diversity_gamma <- Diversity %>% dplyr::group_by(
    Time
  ) %>% dplyr::summarise(
    Mean = mean(Richness),
    SpeciesTotal = toString(sort(unique(unlist(strsplit(paste(
      Species, collapse = ", "), split = ", ", fixed = TRUE))))),
    Gamma = unlist(lapply(strsplit(
      SpeciesTotal, split = ", ", fixed = TRUE), length))
  ) %>% tidyr::pivot_longer(
    cols = c(Mean, Gamma),
    names_to = "Aggregation",
    values_to = "Richness"
  )

  # Combine the two types of results
  Diversity_alpha <- Diversity_alpha %>% dplyr::select(
    -Species
  ) %>% tidyr::pivot_longer(
    cols = c(Richness, Entropy, Evenness),
    names_to = "Measurement",
    values_to = "Value"
  ) %>% dplyr::mutate(
    Environment = as.character(Environment)
  )

  Diversity_gamma <- Diversity_gamma %>% dplyr::select(
    -SpeciesTotal
  ) %>% dplyr::rename(
    Environment = Aggregation,
    Value = Richness
  ) %>% dplyr::mutate(
    Measurement = "Richness"
  )

  Diversity_beta <- Diversity_alpha %>% dplyr::filter(
    Measurement == "Richness"
  ) %>% dplyr::select(
    -Measurement
  ) %>% dplyr::left_join(
    y = Diversity_gamma %>% dplyr::filter(
      Measurement == "Richness", Environment == "Gamma"
    ) %>% dplyr::select(
      -Measurement, -Environment
    ),
    by = "Time",
    suffix = c("_Alpha", "_Gamma")
    # ) %>% dplyr::group_by(
    #   Time
  ) %>% dplyr::mutate(
    BetaSpeciesMissing = Value_Gamma - Value_Alpha,
    BetaSpeciesPercentage = Value_Alpha/Value_Gamma
  ) %>% dplyr::select(
    -Value_Gamma, -Value_Alpha
  ) %>% tidyr::pivot_longer(
    names_to = "Measurement",
    values_to = "Value",
    cols = c(BetaSpeciesMissing, BetaSpeciesPercentage)
    # ) %>% dplyr::ungroup(
  )

  #print(c(colnames(Diversity_alpha), colnames(Diversity_beta), colnames(Diversity_gamma)))
  Diversity <- rbind(
    Diversity_alpha,
    Diversity_beta,
    Diversity_gamma
  )

  return(Diversity)
}

Diversities <- lapply(
  results, Calculate_Diversity
)

# We bin and burn-in, as discussed.
# Look at the ggplot2s to decide burn in. Usually about 100 * 1/CharRate.
# ggplot2::ggplot(
#   diversity1 %>% dplyr::filter(
#     Measurement == "Richness"
#   ), ggplot2::aes(
#     x = Time, y = Value, color = Environment
#   )) + ggplot2::geom_line()

Richnesses <- lapply(
  seq_along(Diversities),
  function(i, div, charrate) {
    div[[i]] %>% dplyr::filter(
      Measurement == "Richness",
      Environment != "Mean",
      Environment != "Gamma",
      Time > 100 / charrate[[i]]
    ) %>% dplyr::mutate(
      Time = floor(Time / 10) * 10
    ) %>% tidyr::pivot_wider(
      names_from = Environment,
      values_from = Value,
      values_fn = median
      # Steps where a species is removed or added seem to
      # be resulting in two entries.
    )
  },
  div = Diversities,
  charrate = CharacteristicRates
)

RichnessesCentred <- lapply(Richnesses, function(x) {
  x %>% dplyr::mutate(
    dplyr::across(`1`:`10`,
                  .fns =  ~ .x - mean(.x, na.rm = TRUE))
  )
})

RichnessesDeAR1ed <- lapply(RichnessesCentred, function(x) {
  x %>% dplyr::mutate(
    dplyr::across(`1`:`10`,
                  .fns =  ~ residuals(arima(.x, order = c(1, 0, 0))))
  )
})

results <- lapply(
  Richnesses,
  function(x) testcorr::rcorr.test(x %>% dplyr::select(`1`:`10`))
  )

MCs <- lapply(Richnesses, function(rich) {
  coda::mcmc.list(
    lapply(1:Environments, function(i, r = rich) {
      coda::mcmc(
        r[[toString(i)]],
        start = r$Time[1],
        end = r$Time[nrow(r)],
        thin = unique(diff(r$Time)))
    })
  )
})

# Not sure there is a real rule of thumb here but a Stata blog claims 1.1.
# https://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/
# It is also pointed out that gelman.plot should be monotonic.
GelRubs <- lapply(MCs, autoburnin = FALSE, coda::gelman.diag)
