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
Space <- match.arg("None", c("None", "Ring", "Line", "Full"))

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

Pool1 <- LawMorton1996_species(
  Basal = Species1[1],
  Consumer = Species1[2],
  Parameters = LMParameters,
  LogBodySize = LMLogBodySize,
  seed = PoolSeed
)

InteractionMatrices1 <- CreateEnvironmentInteractions(
  Pool = Pool1, NumEnvironments = Environments,
  ComputeInteractionMatrix = LawMorton1996_CommunityMat,
  Parameters = LMParameters,
  EnvironmentSeeds = EnvironmentSeed
)

InteractionMatrices1$Mats <- lapply(InteractionMatrices1$Mats,
                                    function(m) m * TimeScaling)

Pool1$ReproductionRate <- - rowMeans(do.call(
  cbind,
  lapply(InteractionMatrices1$Mats,
         function(m) m %*% c(1, 1))
))

stopifnot(
  all(unlist(
    lapply(InteractionMatrices1$Mats,
           function(m) - solve(m) %*% Pool1$ReproductionRate)) > 0)
)

Pool2 <- LawMorton1996_species(
  Basal = Species2[1],
  Consumer = Species2[2],
  Parameters = LMParameters,
  LogBodySize = LMLogBodySize,
  seed = PoolSeed
)

InteractionMatrices2 <- CreateEnvironmentInteractions(
  Pool = Pool2, NumEnvironments = Environments,
  ComputeInteractionMatrix = LawMorton1996_CommunityMat,
  Parameters = LMParameters,
  EnvironmentSeeds = EnvironmentSeed
)

InteractionMatrices2$Mats <- lapply(InteractionMatrices2$Mats,
                                    function(m) m * TimeScaling)

Pool2$ReproductionRate <- - rowMeans(do.call(
  cbind,
  lapply(InteractionMatrices2$Mats,
         function(m) m %*% c(1, 1, 1))
))

stopifnot(
  all(unlist(
    lapply(InteractionMatrices2$Mats,
           function(m) - solve(m) %*% Pool2$ReproductionRate)) > 0)
)

Pool3 <- LawMorton1996_species(
  Basal = Species3[1],
  Consumer = Species3[2],
  Parameters = LMParameters,
  LogBodySize = LMLogBodySize,
  seed = PoolSeed
)

InteractionMatrices3 <- CreateEnvironmentInteractions(
  Pool = Pool3, NumEnvironments = Environments,
  ComputeInteractionMatrix = LawMorton1996_CommunityMat,
  Parameters = LMParameters,
  EnvironmentSeeds = EnvironmentSeed
)

InteractionMatrices3$Mats <- lapply(InteractionMatrices3$Mats,
                                    function(m) m * TimeScaling)

Pool4 <- LawMorton1996_species(
  Basal = Species4[1],
  Consumer = Species4[2],
  Parameters = LMParameters,
  LogBodySize = LMLogBodySize,
  seed = PoolSeed
)

InteractionMatrices4 <- CreateEnvironmentInteractions(
  Pool = Pool4, NumEnvironments = Environments,
  ComputeInteractionMatrix = LawMorton1996_CommunityMat,
  Parameters = LMParameters,
  EnvironmentSeeds = EnvironmentSeed
)

InteractionMatrices4$Mats <- lapply(InteractionMatrices4$Mats,
                                    function(m) m * TimeScaling)

# Pool5 <- LawMorton1996_species(
#   Basal = Species5[1],
#   Consumer = Species5[2],
#   Parameters = LMParameters,
#   LogBodySize = LMLogBodySize,
#   seed = PoolSeed
# )
#
# InteractionMatrices5 <- CreateEnvironmentInteractions(
#   Pool = Pool5, NumEnvironments = Environments,
#   ComputeInteractionMatrix = LawMorton1996_CommunityMat,
#   Parameters = LMParameters,
#   EnvironmentSeeds = EnvironmentSeed
# )
#
# InteractionMatrices5$Mats <- lapply(InteractionMatrices5$Mats,
#                                     function(m) m * TimeScaling)


Pool6 <- LawMorton1996_species(
  Basal = Species6[1],
  Consumer = Species6[2],
  Parameters = LMParameters,
  LogBodySize = LMLogBodySize,
  seed = PoolSeed
)

InteractionMatrices6 <- CreateEnvironmentInteractions(
  Pool = Pool6, NumEnvironments = Environments,
  ComputeInteractionMatrix = LawMorton1996_CommunityMat,
  Parameters = LMParameters,
  EnvironmentSeeds = EnvironmentSeed
)

InteractionMatrices6$Mats <- lapply(InteractionMatrices6$Mats,
                                    function(m) m * TimeScaling)

CharacteristicRate1 <- max(unlist(lapply(
  InteractionMatrices1$Mats, function(m) abs(eigen(m)$values))))
CharacteristicRate2 <- max(unlist(lapply(
  InteractionMatrices2$Mats, function(m) abs(eigen(m)$values))))
CharacteristicRate3 <- max(unlist(lapply(
  InteractionMatrices3$Mats, function(m) abs(eigen(m)$values))))
CharacteristicRate4 <- max(unlist(lapply(
  InteractionMatrices4$Mats, function(m) abs(eigen(m)$values))))
# CharacteristicRate5 <- max(unlist(lapply(
#   InteractionMatrices5$Mats, function(m) abs(eigen(m)$values))))
CharacteristicRate6 <- max(unlist(lapply(
  InteractionMatrices6$Mats, function(m) abs(eigen(m)$values))))
print(c(CharacteristicRate1, CharacteristicRate2,
        CharacteristicRate3, CharacteristicRate4))

Events1 <- CreateAssemblySequence(
  Species = sum(Species1),
  NumEnvironments = Environments,
  ArrivalEvents = EventsEach,
  ArrivalRate = CharacteristicRate1,
  ArrivalFUN = ArrivalFUN_Example,
  ExtinctEvents = EventsEach,
  ExtinctRate = CharacteristicRate1,
  ExtinctFUN = ExtinctFUN_Example,
  HistorySeed = HistorySeed
)

Events2 <- CreateAssemblySequence(
  Species = sum(Species2),
  NumEnvironments = Environments,
  ArrivalEvents = EventsEach,
  ArrivalRate = CharacteristicRate2,
  ArrivalFUN = ArrivalFUN_Example,
  ExtinctEvents = EventsEach,
  ExtinctRate = CharacteristicRate2,
  ExtinctFUN = ExtinctFUN_Example,
  HistorySeed = HistorySeed
)

Events3 <- CreateAssemblySequence(
  Species = sum(Species3),
  NumEnvironments = Environments,
  ArrivalEvents = EventsEach,
  ArrivalRate = CharacteristicRate3,
  ArrivalFUN = ArrivalFUN_Example,
  ExtinctEvents = EventsEach,
  ExtinctRate = CharacteristicRate3,
  ExtinctFUN = ExtinctFUN_Example,
  HistorySeed = HistorySeed
)

Events4 <- CreateAssemblySequence(
  Species = sum(Species4),
  NumEnvironments = Environments,
  ArrivalEvents = EventsEach,
  ArrivalRate = CharacteristicRate4,
  ArrivalFUN = ArrivalFUN_Example,
  ExtinctEvents = EventsEach,
  ExtinctRate = CharacteristicRate4,
  ExtinctFUN = ExtinctFUN_Example,
  HistorySeed = HistorySeed
)

# Events5 <- CreateAssemblySequence(
#   Species = sum(Species5),
#   NumEnvironments = Environments,
#   ArrivalEvents = Environments * ceiling(100 * (log(100 + 5))),
#   ArrivalRate = CharacteristicRate5,
#   ArrivalFUN = ArrivalFUN_Example,
#   ExtinctEvents = Environments * ceiling(100 * (log(100 + 5))),
#   ExtinctRate = CharacteristicRate5,
#   ExtinctFUN = ExtinctFUN_Example,
#   HistorySeed = HistorySeed
# )

Events6 <- CreateAssemblySequence(
  Species = sum(Species6),
  NumEnvironments = Environments,
  ArrivalEvents = EventsEach,
  ArrivalRate = CharacteristicRate6,
  ArrivalFUN = ArrivalFUN_Example,
  ExtinctEvents = EventsEach,
  ExtinctRate = CharacteristicRate6,
  ExtinctFUN = ExtinctFUN_Example,
  HistorySeed = HistorySeed
)

print(table(Events1$Events$Species,
            Events1$Events$Environment))
print(table(Events2$Events$Species,
            Events2$Events$Environment))
print(table(Events3$Events$Species,
            Events3$Events$Environment))
print(table(Events4$Events$Species,
            Events4$Events$Environment))

IntMat1 <- Matrix::bdiag(InteractionMatrices1$Mats)
IntMat2 <- Matrix::bdiag(InteractionMatrices2$Mats)
IntMat3 <- Matrix::bdiag(InteractionMatrices3$Mats)
IntMat4 <- Matrix::bdiag(InteractionMatrices4$Mats)
# IntMat5 <- Matrix::bdiag(InteractionMatrices5$Mats)
IntMat6 <- Matrix::bdiag(InteractionMatrices6$Mats)
PerCapitaDynamics1 <- PerCapitaDynamics_Type1(
  Pool1$ReproductionRate, IntMat1,
  NumEnvironments = Environments
)
PerCapitaDynamics2 <- PerCapitaDynamics_Type1(
  Pool2$ReproductionRate, IntMat2,
  NumEnvironments = Environments
)
PerCapitaDynamics3 <- PerCapitaDynamics_Type1(
  Pool3$ReproductionRate, IntMat3,
  NumEnvironments = Environments
)
PerCapitaDynamics4 <- PerCapitaDynamics_Type1(
  Pool4$ReproductionRate, IntMat4,
  NumEnvironments = Environments
)
# PerCapitaDynamics5 <- PerCapitaDynamics_Type1(
#   Pool5$ReproductionRate, IntMat5,
#   NumEnvironments = Environments
# )
PerCapitaDynamics6 <- PerCapitaDynamics_Type1(
  Pool6$ReproductionRate, IntMat6,
  NumEnvironments = Environments
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

DispersalMatrix1 <- CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool1))
)

DispersalMatrix2 <- CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool2))
)

DispersalMatrix3 <- CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool3))
)

DispersalMatrix4 <- CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool4))
)

# DispersalMatrix5 <- CreateDispersalMatrix(
#   EnvironmentDistances = DistanceMatrix,
#   SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool5))
# )

DispersalMatrix6 <- CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool6))
)


print(Sys.time())
result1 <- MultipleNumericalAssembly_Dispersal(
  Pool = Pool1, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices1,
  Events = Events1,
  PerCapitaDynamics = PerCapitaDynamics1,
  DispersalMatrix = DispersalMatrix1,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = FALSE
)
print(Sys.time())
result2 <- MultipleNumericalAssembly_Dispersal(
  Pool = Pool2, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices2,
  Events = Events2,
  PerCapitaDynamics = PerCapitaDynamics2,
  DispersalMatrix = DispersalMatrix2,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = FALSE
)
print(Sys.time())
result3 <- MultipleNumericalAssembly_Dispersal(
  Pool = Pool3, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices3,
  Events = Events3,
  PerCapitaDynamics = PerCapitaDynamics3,
  DispersalMatrix = DispersalMatrix3,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = FALSE
)
print(Sys.time())
result4 <- MultipleNumericalAssembly_Dispersal(
  Pool = Pool4, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices4,
  Events = Events4,
  PerCapitaDynamics = PerCapitaDynamics4,
  DispersalMatrix = DispersalMatrix4,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = FALSE
)
print(Sys.time())
# result5 <- MultipleNumericalAssembly_Dispersal(
#   Pool = Pool5, NumEnvironments = Environments,
#   InteractionMatrices = InteractionMatrices5,
#   Events = Events5,
#   PerCapitaDynamics = PerCapitaDynamics5,
#   DispersalMatrix = DispersalMatrix5,
#   EliminationThreshold = EliminationThreshold,
#   ArrivalDensity = ArrivalDensity,
#   MaximumTimeStep = MaximumTimeStep,
#   BetweenEventSteps = BetweenEventSteps,
#   Verbose = FALSE
# )
# print(Sys.time())
result6 <- MultipleNumericalAssembly_Dispersal(
  Pool = Pool6, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices6,
  Events = Events6,
  PerCapitaDynamics = PerCapitaDynamics6,
  DispersalMatrix = DispersalMatrix6,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = FALSE
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
      SpeciesTotal, split = ", ", fixed = TRUE), function(x) length(x[x!=""]) ))
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

diversity1 <- Calculate_Diversity(result1)
diversity2 <- Calculate_Diversity(result2)
diversity3 <- Calculate_Diversity(result3)
diversity4 <- Calculate_Diversity(result4)
# diversity5 <- Calculate_Diversity(result5)
diversity6 <- Calculate_Diversity(result6)

# We bin and burn-in, as discussed.
# Look at the ggplot2s to decide burn in.
# ggplot2::ggplot(
#   diversity1 %>% dplyr::filter(
#     Measurement == "Richness"
#   ), ggplot2::aes(
#     x = Time, y = Value, color = Environment
#   )) + ggplot2::geom_line()

richness1 <- diversity1 %>% dplyr::filter(
  Measurement == "Richness",
  Environment != "Mean",
  Environment != "Gamma",
  Time > 5000
) %>% dplyr::mutate(
  Time = floor(Time / 10) * 10
) %>% tidyr::pivot_wider(
  names_from = Environment,
  values_from = Value,
  values_fn = median
  # Steps where a species is removed or added seem to
  # be resulting in two entries.
)
richness2 <- diversity2 %>% dplyr::filter(
  Measurement == "Richness",
  Environment != "Mean",
  Environment != "Gamma",
  Time > 5000
) %>% dplyr::mutate(
  Time = floor(Time / 10) * 10
)  %>% tidyr::pivot_wider(
  names_from = Environment,
  values_from = Value,
  values_fn = median
)
richness3 <- diversity3 %>% dplyr::filter(
  Measurement == "Richness",
  Environment != "Mean",
  Environment != "Gamma",
  Time > 2.5E4
) %>% dplyr::mutate(
  Time = floor(Time / 10) * 10
)  %>% tidyr::pivot_wider(
  names_from = Environment,
  values_from = Value,
  values_fn = median
)
richness4 <- diversity4 %>% dplyr::filter(
  Measurement == "Richness",
  Environment != "Mean",
  Environment != "Gamma",
  Time > 5E4
)  %>% dplyr::mutate(
  Time = floor(Time / 10) * 10
) %>% tidyr::pivot_wider(
  names_from = Environment,
  values_from = Value,
  values_fn = median
)
# richness5 <- diversity5 %>% dplyr::filter(
#   Measurement == "Richness",
#   Environment != "Mean",
#   Environment != "Gamma",
#   Time > 5E4 # TODO: The big question
# )  %>% dplyr::mutate(
#   Time = floor(Time / 10) * 10
# ) %>% tidyr::pivot_wider(
#   names_from = Environment,
#   values_from = Value,
#   values_fn = median
# )
richness6 <- diversity6 %>% dplyr::filter(
  Measurement == "Richness",
  Environment != "Mean",
  Environment != "Gamma",
  Time > 5E4
)  %>% dplyr::mutate(
  Time = floor(Time / 10) * 10
) %>% tidyr::pivot_wider(
  names_from = Environment,
  values_from = Value,
  values_fn = median
)

results <- lapply(
  list(richness1, richness2, richness3, richness4 #,
       # richness5
       ),
  function(x) testcorr::rcorr.test(x %>% dplyr::select(`1`:`10`)))
