library(RMTRCode2)
library(parallel)
library(doParallel)
library(foreach)

clust <- parallel::makeCluster(14)
doParallel::registerDoParallel(clust)

# Parameters: ##################################################################
Species <- c(Basal = 34, Consumer = 66)
Environments <- 10
EventsBase <- Environments * ceiling(sum(Species) * (log(sum(Species) + 5)))

LMParameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
LMLogBodySize <- c(-2, -1, -1, 0)

SpeciesSpeeds <- 1
Space <- match.arg("Line", c("None", "Ring", "Line", "Full"))

EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 10 # Number of steps to reach next event to smooth.

CalculatePoolAndMatrices <- TRUE
dir <- getSrcDirectory(function(){})

# > runif(1) * 1e8
# [1] 28665115q
PoolSeed <- 28665115
# > runif(1) * 1e8
# [1] 75027622
EnvironmentSeed <- 75027622

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
         "MNA-2",
         "-PoolMats-Env", Environments,
         ".RData")
       ))
  
} else {
  load(file = file.path(dir, paste0(
    "MNA-2", 
    "-PoolMats-Env", Environments,
    ".RData")
  ))
}

# Technically, not quite right, but should be good enough.
# In particular,
# max(unlist(lapply(InteractionMatrices$Mats, function(m) {abs(eigen(m)$values)})))
# gives nearly the same value.
# Also note: eigenvalues of block matrices are the eigenvalues of the blocks.
CharacteristicRate <- max(unlist(lapply(
  InteractionMatrices$Mats, function(m) {abs(eigen(m)$values)}
)))

IntMat <- Matrix::bdiag(InteractionMatrices$Mats)

PerCapitaDynamics <- PerCapitaDynamics_Type1(
  Pool$ReproductionRate, IntMat,
  NumEnvironments = Environments
)

parametersDF <- expand.grid(
  EventsArrRate = 10^(-1:1),
  EventsExtRate = 10^(-1:1),
  PerIslandDistance = c(10^(3*(0:3)), Inf)
)

# > runif(1)*1E8
# [1] 48175667
if (exists(".Random.seed")) {
  oldSeed <- .Random.seed
}
set.seed(48175667)
parametersDF$EventsArr <- parametersDF$EventsArrRate * EventsBase
parametersDF$EventsExt <- parametersDF$EventsExtRate * EventsBase
parametersDF$HistorySeed <- runif(nrow(parametersDF)) * 1E8
if (exists("oldSeed")) {
  set.seed(oldSeed)
} else {
  set.seed()
}

parallel::clusterExport(
  clust, varlist = c(
    "Pool", "Species", "Environments", "CharacteristicRate",
    "SpeciesSpeeds", "PerCapitaDynamics", "IntMat",
    "EliminationThreshold",
    "ArrivalDensity",
    "MaximumTimeStep",
    "BetweenEventSteps"
  )
  )

success <- foreach::foreach(
  params = iterators::iter(parametersDF, by = "row")
) %dopar% {
  tryCatch({
    
  Events <- with(params, CreateAssemblySequence(
    Species = sum(Species),
    NumEnvironments = Environments,
    ArrivalEvents = EventsArr,
    ArrivalRate = EventsArrRate * CharacteristicRate,
    ArrivalFUN = RMTRCode2::ArrivalFUN_Example2,
    ExtinctEvents = EventsExt,
    ExtinctRate = EventsExtRate * CharacteristicRate,
    ExtinctFUN = RMTRCode2::ExtinctFUN_Example2,
    HistorySeed = HistorySeed
  ))
  
  if (Space == "None") {
    DistanceMatrix <- Matrix::sparseMatrix(
      i = Environments, j = Environments, x = 0)
  }
  if (Space == "Ring" || Space == "Line")
    DistanceMatrix <- with(params, Matrix::bandSparse(
      Environments, k = c(-1, 1),
      diagonals = list(rep(PerIslandDistance, Environments - 1),
                       rep(PerIslandDistance, Environments - 1))
    ))
  if (Space == "Ring") {
    DistanceMatrix[Environments, 1] <- params$PerIslandDistance
    DistanceMatrix[1, Environments] <- params$PerIslandDistance
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
  
  result <- with(params, MultipleNumericalAssembly_Dispersal(
    Pool = Pool, NumEnvironments = Environments, 
    CharacteristicRate = CharacteristicRate, 
    Events = Events,
    PerCapitaDynamics = PerCapitaDynamics,
    DispersalMatrix = DispersalMatrix,
    EliminationThreshold = EliminationThreshold,
    ArrivalDensity = ArrivalDensity,
    MaximumTimeStep = MaximumTimeStep,
    BetweenEventSteps = BetweenEventSteps,
    Verbose = FALSE
  ))
  
  with(params, save(result,
       file = file.path(dir, paste0(
         "MNA-2-Dist",
         format(PerIslandDistance, scientific = TRUE),
         "-", "Arr", EventsArrRate,
         "-", "Ext", EventsExtRate,
         "-Env", Environments,
         "-", Space,
         ".RData")
       )
  ))
  
  return(TRUE)
  }, error = function(e) {
    return(e)
  })
}

parametersDF$success <- success

print(parametersDF)

parallel::stopCluster(clust)