# We are going to look for a few stable systems specifically and see how the
# known-independent runs are correlated when the pool is small (and should
# only have a single known steady-state with no transitions between
# steady-states).
# Since we only want stable systems, we can solve for the appropriate
# reproductive rates in order to maintain a system of the given sizes.
library(RMTRCode2)

Species1 <- c(Basal = 1, Consumer = 1)
Species2 <- c(Basal = 1, Consumer = 2)
Species3 <- c(Basal = 3, Consumer = 0)
Species4 <- c(Basal = 10, Consumer = 0)
Environments <- 10
EventsEach <- Environments * ceiling(10 * (log(10 + 5)))

LMParameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
LMLogBodySize <- c(-2, -1, -1, 0)

PerIslandDistance <- Inf
SpeciesSpeeds <- 1
Space <- match.arg("Ring", c("None", "Ring", "Line", "Full"))

EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 10 # Number of steps to reach next event to smooth.

# > runif(1) * 1e8
# [1] 38427042
PoolSeed <- 38427042
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
