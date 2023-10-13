library(RMTRCode2)
source("Mutualism.R")

# Data_2023-05-02 had almost all species able to live together.
# Afterwards, I increased the competition by a factor of 10 while keeping
# intraspecific:interspecific competition ratio the same.

# Parameters: ##################################################################
Species <- c(Producer = 34, Pollinator = 66) * 2
Environments <- 10
EventsEach <- Environments * ceiling(sum(Species) * (log(sum(Species) + 5)))
EventRateModifiers <- c(1, 1) # Immigration, Extirpation


PerIslandDistance <- 10^1 # 10^5 # Inf # 10^0
SpeciesSpeeds <- 1
Space <- match.arg("Ring", c("None", "Ring", "Line", "Full"))

EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007
ExtinctionProportion <- 1

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 10 # Number of steps to reach next event to smooth.

CalculatePoolAndMatrices <- FALSE
dir <- paste0("Data_", Sys.Date()) # getSrcDirectory(function(){})

if (!dir.exists(dir)) {
  dir.create(dir, showWarnings = FALSE)
}

# > runif(1) * 1e8
# [1] 9927818
PoolSeed <- 9927818
# > runif(1) * 1e8
# [1] 97403448
EnvironmentSeed <- 97403448
# > runif(1) * 1e8
# [1] 94001959
HistorySeed <- 94001959
# > runif(1) * 1e8
# [1] 91068274
SatSeed <- 91068274

# Setup: #######################################################################

## Pools and Interaction Matrices: #############################################
if (CalculatePoolAndMatrices) {
  Pool <- Mutualism_species(
    SpeciesTypes = Species,
    seed = PoolSeed
  )

  Saturation <- Mutualism_saturation(
    Species, seed = SatSeed
  )

  InteractionMatrices <- RMTRCode2::CreateEnvironmentInteractions(
    Pool = Pool, NumEnvironments = Environments,
    ComputeInteractionMatrix = Mutualism_CommunityMat,
    EnvironmentSeeds = EnvironmentSeed
  )
  save(Pool, Saturation, InteractionMatrices,
       file = file.path(dir, paste0(
         "Mutualism-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
} else {
  load(file = file.path(dir, paste0(
    "Mutualism-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
}

## Events: #####################################################################

# Note: eigenvalues of block matrices are the eigenvalues of the blocks.
CharacteristicRate <- max(unlist(lapply(
  InteractionMatrices$Mats, function(m) {abs(eigen(m)$values)}
)))

Events <- RMTRCode2::CreateAssemblySequence(
  Species = sum(Species),
  NumEnvironments = Environments,
  ArrivalEvents = EventsEach * EventRateModifiers[1],
  ArrivalRate = CharacteristicRate * EventRateModifiers[1],
  ArrivalFUN = RMTRCode2::ArrivalFUN_Example2,
  ExtinctEvents = EventsEach * EventRateModifiers[2],
  ExtinctRate = CharacteristicRate * EventRateModifiers[2],
  ExtinctFUN = RMTRCode2::ExtinctFUN_Example2,
  HistorySeed = HistorySeed
)

print(table(Events$Events$Species,
            Events$Events$Environment))

## Dynamics: ###################################################################

IntMat <- Matrix::bdiag(InteractionMatrices$Mats)
PerCapitaDynamics <- PerCapitaDynamics_Mutualistic1(
  Pool$ReproductionRate, IntMat,
  NumEnvironments = Environments,
  SpeciesTypes = Species, Saturations = Saturation
)

### Spatial/Dispersal: #########################################################
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

DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(SpeciesSpeeds, nrow(Pool))
)

## Run: #######################################################################
print(Sys.time())
if (exists("MultipleNumericalAssembly_Dispersal")) {
  theFun <- MultipleNumericalAssembly_Dispersal
} else {
  theFun <- RMTRCode2::MultipleNumericalAssembly_Dispersal
}

result <- theFun(
  Pool = Pool, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices,
  Events = Events,
  PerCapitaDynamics = PerCapitaDynamics,
  DispersalMatrix = DispersalMatrix,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  ExtinctionProportion = ExtinctionProportion,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = TRUE
)
print(Sys.time())

save(result,
     file = file.path(dir, paste0(
       "Mutualism-ExampleExtProp-Result-Env", Environments,
       "-", Space,
       "-", gsub(round(log10(PerIslandDistance)),
                 pattern = "-", replacement = "_", fixed = TRUE),
       "-", EventRateModifiers[1], "-", EventRateModifiers[2],
       "-ExtProp", ExtinctionProportion, ".RData")
     )
)
print(Sys.time())
