# NOTE: Major differences from T&F 2010:
# Not worrying about varying connectivity in this script. With 200 species
# T&F suggest around 0.05 connectivity, but that's for a final network.
# Between 20 and 50 species occur with a 0.20 connectance.

# Libraries: ##################################################################
#print("Loading libraries")
#librarypath <- file.path(".", "Rlibs")
#if (!dir.exists(librarypath)) {
#  dir.create(librarypath, showWarnings = FALSE)
#}
#.libPaths(c(librarypath, .libPaths()))
#
#allLibraryPaths <- .libPaths()
#
#packages <- c(
#  "Matrix",     # Common Format
#  "parallel",   # Base parallel dependency
#  "doParallel", # For compatible cluster
#  "foreach",    # For parallel evaluation
#  "iterators"   # For parallel evaluation
#)
#
## Does not work because it tries to use the system-wide libraries... Oops.
## update.packages(repos = 'https://cloud.r-project.org',
##                 oldPkgs = packages, ask = FALSE)
#
#for (package in packages) {
#  if (!require(package, character.only = TRUE)) {
#    install.packages(package, lib = librarypath,
#                     repos = 'https://cloud.r-project.org',
#                     dependencies = TRUE)
#  }
#  library(package, character.only = TRUE)
#}
#
#if (!require("RMTRCode2", character.only = TRUE)) {
#  install.packages(
#    "RMTRCode2_0.3.tar.gz", lib = librarypath,
#    repos = NULL, type = "source"
#  )
#}
#library(RMTRCode2)#, lib.loc = librarypath) # lib.loc shouldn't be necessary.
#
#cargs <- as.numeric(commandArgs(trailingOnly = TRUE))
#
library(RMTRCode2)
library(Matrix)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

clust <- parallel::makeCluster(1, outfile = "")
doParallel::registerDoParallel(clust)
parallel::clusterCall(clust, function() {
librarypath <- file.path(".", "Rlibs")
.libPaths(c(librarypath, .libPaths()))
library(RMTRCode2)
library(Matrix)
}) # Implicitly needed.


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
dir <- paste0("Data_", "2023-08-01")#Sys.Date()) # "2023-06-30")
# See also 07-24, which we are going to compare with to make sure equivalent.

if (!dir.exists(dir)) {
  dir.create(dir, showWarnings = FALSE)
}

# > runif(4) * 1e8
# [1] 54668004 73065730 63168946  4069184
PoolSeed <- 54668004
EnvironmentSeed <- 73065730
HistorySeed <- 63168946
SatSeed <- 4069184

# Setup: #######################################################################

## Pools and Interaction Matrices: #############################################
if (CalculatePoolAndMatrices) {
  Pool <- Mutualism_species(
    SpeciesTypes = Species,
    MinimumRepRates = c(0.3, -0.2),
    MaximumRepRates = c(0.5, 0),
    seed = PoolSeed
  )

  Saturation <- Mutualism_saturation(
    rep(1, sum(Species)), seed = SatSeed
  )

  InteractionMatrices <- RMTRCode2::CreateEnvironmentInteractions(
    Pool = Pool, NumEnvironments = Environments,
    ComputeInteractionMatrix = Mutualism_CommunityMat_ByBlock,
    EnvironmentSeeds = EnvironmentSeed,
    MinimumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
                                c(0, -3, 2, 0)),
    MaximumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
                                c(0, -2, 3, 0))
  )

  # To induce the proper connectance, you can Hadamard product an undirected
  # network with each interaction matrix. Due to the block structure, we need
  # to use the same block structure in the adj. matrix. We also fix diagonals.
  adjmat <- CreateConnectanceByBlock(
    Species,
    Connectances = matrix(c(0, 0.2, NA, 0),
                          byrow = TRUE, nrow = 2),
    Directed = FALSE
  )

  InteractionMatrices$Mats <- lapply(InteractionMatrices$Mats, function(mat) {
    mat <- mat * adjmat
    diag(mat) <- runif(nrow(mat), min = -2, max = -1)
    return(mat)
  })

  # T&F use the connectance of the saturations as the primary connectance...
  diagSat <- diag(Saturation)
  Saturation <- Saturation * adjmat
  diag(Saturation) <- diagSat

  # Additionally, Saturation needs to be species x species x patch:
  Saturation <- Matrix::bdiag(rep(list(Saturation), Environments))

  save(Pool, Saturation, InteractionMatrices,
       file = file.path(dir, paste0(
         "PrPrTF2010-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
} else {
  load(file = file.path(dir, paste0(
    "PrPrTF2010-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
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
PerCapitaDynamics <- PerCapitaDynamics_Mutualistic4(
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
         "PrPrTF2010-ExampleExtProp-Result-Env", Environments,
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
