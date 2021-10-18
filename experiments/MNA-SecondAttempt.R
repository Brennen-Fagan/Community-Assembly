# Second attempt at generating data.
# Broadly, first attempt was successful, but we ran into a few issues.
# The primary issue is probably the lack of planning regarding file names.
# A secondary issue is that, while the code works, the first thread needs to
# create the pool and matrices for all following threads.
# Another secondary issue is that we are dealing so far with only one
# realisation, but we could have any number of realised histories or systems
# and it is not obvious as to when one should stop generating either.
# This is, of course, before we even consider proper variations.
# Perhaps the most obvious thing to do then would be to do 10 histories for each
# of 10 systems for each parameter combination, in order to be safe.
# Each of the "10 systems for each parameter combination" then needs its pools
# and matrices initialised first before the history is generated.
# Then each history (and pool and matrix system) is then run on each of the
# spatial treatments.

library(RMTRCode2)

# Parameters: ##################################################################
Species <- c(Basal = 34, Consumer = 66)
Environments <- 10
EventsEach <- Environments * ceiling(sum(Species) * (log(sum(Species) + 5)))

LMParameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
LMLogBodySize <- c(-2, -1, -1, 0)

PerIslandDistance <- 1E9
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
# [1] 28665115
PoolSeed <- 28665115
# > runif(1) * 1e8
# [1] 75027622
EnvironmentSeed <- 75027622
# > runif(1) * 1e8
# [1] 64713671
HistorySeed <- 64713671

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
         "MNA-FirstAttempt-PoolMats-Env", Environments, ".RData")))
} else {
  load(file = file.path(dir, paste0(
    "MNA-FirstAttempt-PoolMats-Env", Environments, ".RData")))
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
  ArrivalEvents = EventsEach,
  ArrivalRate = CharacteristicRate,
  ArrivalFUN = ArrivalFUN_Example,
  ExtinctEvents = EventsEach,
  ExtinctRate = CharacteristicRate,
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
       "MNA-FirstAttempt",
       format(PerIslandDistance, scientific = TRUE),
       "-Result-Env", Environments,
       "-", Space, ".RData")
     )
)
