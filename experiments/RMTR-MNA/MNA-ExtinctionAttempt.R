library(RMTRCode2)
# source('~/Random-Matrix-Theory/RMTRCode2/R/MultipleCommunityAssembly.R', echo=FALSE)
# source('~/Random-Matrix-Theory/RMTRCode2/R/Tests-MultipleCommunityAssembly.R', echo=FALSE)

# Parameters: ##################################################################
Species <- c(Basal = 34, Consumer = 66)
Environments <- 10
EventsArr <- Environments * ceiling(sum(Species) * (log(sum(Species) + 5)))
EventsExtRate <- 1/10 # Relative to Arrival Event Rate.
EventsExt <- EventsExtRate * EventsArr
# If 1/10 the arrival rate, need 1/10 the events to reach the same length.
# If 10 * the arrival rate, need 10 * the events to reach the same length.

LMParameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
LMLogBodySize <- c(-2, -1, -1, 0)

PerIslandDistance <- 1E9 # of 1 (None,Line), 1E3, 1E6, 1E9 (Line)
# About 6 orders of magnitude diff between high abund and elimination.
SpeciesSpeeds <- 1
Space <- match.arg("Line", c("None", "Ring", "Line", "Full"))

EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 10 # Number of steps to reach next event to smooth.

CalculatePoolAndMatrices <- FALSE
dir <- getSrcDirectory(function(){})

# > runif(1) * 1e8
# [1] 28665115q
PoolSeed <- 28665115
# > runif(1) * 1e8
# [1] 75027622
EnvironmentSeed <- 75027622
# > runif(1) * 1e8
# [1] 27746746
# > runif(1) * 1e8
# [1] 82176895
if (EventsExtRate == 10) {
  HistorySeed <- 27746746
} else if (EventsExtRate == 0.1) {
  HistorySeed <- 82176895
}

# Setup: #######################################################################
if (CalculatePoolAndMatrices) {
  Pool <- LawMorton1996_species(
    Basal = Species[1],
    Consumer = Species[2],
    Parameters = LMParameters,
    LogBodySize = LMLogBodySize,
    seed = PoolSeed
  )

  InteractionMatrices <- CreateEnvironmentInteractions(
    Pool = Pool, NumEnvironments = Environments,
    ComputeInteractionMatrix = LawMorton1996_CommunityMat,
    Parameters = LMParameters,
    EnvironmentSeeds = EnvironmentSeed
  )
  save(Pool, InteractionMatrices,
       file = file.path(dir, paste0(
         "MNA-Ext", EventsExtRate,
         "-PoolMats-Env", Environments,
         ".RData")
       ))

} else {
  load(file = file.path(dir, paste0(
    "MNA-Ext", EventsExtRate,
    "-PoolMats-Env", Environments,
    ".RData")
  ))
}

# Technically, not quite right, but should be good enough.
# In particular,
# max(unlist(lapply(InteractionMatrices$Mats, function(m) {abs(eigen(m)$values)})))
# gives nearly the same value.
# Also note: eigenvalues of block matrices are the eigenvalues of the blocks.
CharacteristicRate <- max(abs(eigen(InteractionMatrices$Mats[[1]])$values))

Events <- CreateAssemblySequence(
  Species = sum(Species),
  NumEnvironments = Environments,
  ArrivalEvents = EventsArr,
  ArrivalRate = CharacteristicRate,
  ArrivalFUN = ArrivalFUN_Example,
  ExtinctEvents = EventsExt,
  ExtinctRate = EventsExtRate * CharacteristicRate,
  ExtinctFUN = ExtinctFUN_Example,
  HistorySeed = HistorySeed
)

print(table(Events$Events$Species,
            Events$Events$Environment))

IntMat <- Matrix::bdiag(InteractionMatrices$Mats)
PerCapitaDynamics <- PerCapitaDynamics_Type1(
  Pool$ReproductionRate, IntMat,
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

DispersalMatrix <- CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool))
)

print(Sys.time())
result <- MultipleNumericalAssembly_Dispersal(
  Pool = Pool, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices,
  Events = Events,
  PerCapitaDynamics = PerCapitaDynamics,
  DispersalMatrix = DispersalMatrix,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = FALSE
)
print(Sys.time())

save(result,
     file = file.path(dir, paste0(
       "MNA-Dist",
       format(PerIslandDistance, scientific = TRUE),
       "-Ext", EventsExtRate, "-Env", Environments,
       "-", Space, ".RData")
     )
)
