library(RMTRCode2)
library(Matrix)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

# Parameters: ##################################################################
Species <- c(Basal = 34, Consumer = 66) * 2
Environments <- 10
EventsEach <- Environments * ceiling(sum(Species) * (log(sum(Species) + 0)))
EventRateModifiers <- c(1, 1) # Immigration, Extirpation

LMParameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
LMLogBodySize <- c(-2, -1, -1, 0)

PerIslandDistance <- 10^c(Inf, 9:-1) # 10^5 # Inf # 10^0
SpeciesSpeeds <- 1
Space <- match.arg("Ring", c("None", "Ring", "Line", "Full"))

EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007
ExtinctionProportion <- 1

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 10 # Number of steps to reach next event to smooth.

CalculatePoolAndMatrices <- TRUE
dir <- paste0("Data_", "2023-09-24")#Sys.Date()) # getSrcDirectory(function(){})

if (!dir.exists(dir)) {
  dir.create(dir, showWarnings = FALSE)
}

# > runif(3) * 1e8
# [1] 11365664 91994571 20423344 # Data_2023-07-06
# [1] 65566924 64305636 14447307 # Data_2023-09-22
# [1] 71113291 29907014 76606233 # Data_2023-09-23
# [1] 53606086 70944574 40035408 # Data_2023-09-24
seeds <- c(
  # 11365664, 91994571, 20423344 # Data_2023-07-06
  # 65566924, 64305636, 14447307 # Data_2023-09-22
  # 71113291, 29907014, 76606233 # Data_2023-09-23
  53606086, 70944574, 40035408 # Data_2023-09-24
)
PoolSeed <- seeds[1]
EnvironmentSeed <- seeds[2]
HistorySeed <- seeds[3]

# Setup: #######################################################################

## Pools and Interaction Matrices: #############################################
if (CalculatePoolAndMatrices) {
  Pool <- RMTRCode2::LawMorton1996_species(
    Basal = Species[1],
    Consumer = Species[2],
    Parameters = LMParameters,
    LogBodySize = LMLogBodySize,
    seed = PoolSeed
  )

  InteractionMatrices <- RMTRCode2::CreateEnvironmentInteractions(
    Pool = Pool, NumEnvironments = Environments,
    ComputeInteractionMatrix = RMTRCode2::LawMorton1996_CommunityMat,
    Parameters = LMParameters,
    EnvironmentSeeds = EnvironmentSeed
  )
  save(Pool, InteractionMatrices,
       file = file.path(dir, paste0(
         "MNA-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
} else {
  load(file = file.path(dir, paste0(
    "MNA-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
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

print(combinations <-
        table(Events$Events$Species,
              Events$Events$Environment))
if(any(combinations == 0)) {warning(
  "Exists a species which doesn't immigrate to an environment."
)}

## Dynamics: ###################################################################

IntMat <- Matrix::bdiag(InteractionMatrices$Mats)
PerCapitaDynamics <- RMTRCode2::PerCapitaDynamics_Type1(
  Pool$ReproductionRate, IntMat,
  NumEnvironments = Environments
)

records <- foreach::foreach(
  dist = iterators::iter(PerIslandDistance),
  # .export = c(
  #   "Pool", "Environments",
  #   "InteractionMatrices",
  #   "Events",
  #   "PerCapitaDynamics",
  #   "EliminationThreshold",
  #   "ArrivalDensity",
  #   "ExtinctionProportion",
  #   "MaximumTimeStep",
  #   "BetweenEventSteps"),
  .export = "reprate",
  .packages = "RMTRCode2"
) %dopar% {
  print(paste(dist, "in"))
  # table(Pool)
  pool <- data.frame(n = 1:sum(Species))#Hack
  # print(paste(dist, "hack"))
  ### Spatial/Dispersal: #########################################################
  if (Space == "None") {
    DistanceMatrix <- Matrix::sparseMatrix(
      i = Environments, j = Environments, x = 0)
  }
  if (Space == "Ring" || Space == "Line")
    DistanceMatrix <- Matrix::bandSparse(
      Environments, k = c(-1, 1),
      diagonals = list(rep(dist, Environments - 1),
                       rep(dist, Environments - 1))
    )
  if (Space == "Ring") {
    DistanceMatrix[Environments, 1] <- dist
    DistanceMatrix[1, Environments] <- dist
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
    SpeciesSpeeds = rep(SpeciesSpeeds, nrow(pool))
  )

  # print(paste(dist, "disp mat"))
  ## Run: #######################################################################
  print(Sys.time())
  if (exists("MultipleNumericalAssembly_Dispersal")) {
    theFun <- MultipleNumericalAssembly_Dispersal
  } else {
    theFun <- RMTRCode2::MultipleNumericalAssembly_Dispersal
  }
  # print(paste(dist, "fun"))
  print(record <- Sys.time())

  result <- theFun(
    pool, NumEnvironments = Environments,
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

  save(result,
       file = file.path(dir, paste0(
         "MNA-ExampleExtProp-Result-Env", Environments,
         "-", Space,
         "-", gsub(round(log10(dist)),
                   pattern = "-", replacement = "_", fixed = TRUE),
         "-", EventRateModifiers[1], "-", EventRateModifiers[2],
         "-ExtProp", ExtinctionProportion, ".RData")
       )
  )
  print(Sys.time())

  return(Sys.time() - record)
}

#parallel::stopCluster(clust)
