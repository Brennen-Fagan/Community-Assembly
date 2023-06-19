library(RMTRCode2)
library(parallel)
library(doParallel)
library(iterators)
library(foreach)

clust <- parallel::makeCluster(3, outfile = "")
doParallel::registerDoParallel(clust)

# Parameters: ##################################################################
Species <- c(Producer = 34, Pollinator = 66) * 2
Environments <- 10
EventsEach <- Environments * ceiling(sum(Species) * (log(sum(Species)) + 0))
EventRateModifiers <- c(1, 1) # Immigration, Extirpation


PerIslandDistance <- 10^c(Inf, 9:-1) # 10^5 # Inf # 10^0
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
# [1] 12032489
PoolSeed <- 12032489
# > runif(1) * 1e8
# [1] 28665115
EnvironmentSeed <- 28665115
# > runif(1) * 1e8
# [1] 75027622
HistorySeed <- 75027622
# > runif(1) * 1e8
# [1] 64713671
SatSeed <- 64713671

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
    ComputeInteractionMatrix = Mutualism_CommunityMat_ByBlock,
    EnvironmentSeeds = EnvironmentSeed,
    MinimumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
                                c(-2/2.5, 2/2.5, 2/2.5, -2/2.5)),
    MaximumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
                                c(-1/2.5, 3/2.5, 3/2.5, -1/2.5))
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

print(combinations <-
        table(Events$Events$Species,
            Events$Events$Environment))
if(any(combinations == 0)) {warning(
  "Exists a species which doesn't immigrate to an environment."
)}

## Dynamics: ###################################################################

IntMat <- Matrix::bdiag(InteractionMatrices$Mats)
reprate <- Pool$ReproductionRate
PerCapitaDynamics <- PerCapitaDynamics_Mutualistic1(
  reprate, IntMat,
  NumEnvironments = Environments,
  SpeciesTypes = Species, Saturations = Saturation
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
  print(record <- Sys.time())

  save(result,
       file = file.path(dir, paste0(
         "Mutualism-ExampleExtProp-Result-Env", Environments,
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

parallel::stopCluster(clust)
