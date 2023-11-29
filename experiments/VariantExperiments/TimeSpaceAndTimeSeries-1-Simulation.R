# Introduction: ###############################################################
# As a sequel to TimeSpaceAndTimeSeries-1-Bootstrap.R and -1-Intervention.R,
# we are now introducing a simulation based intervention.
# Please see the previous files for some design choices, although we aim to
# improve design at each stage.

# Create a Master seed, which we'll use to generate simulation seeds.
# runif(1) * 1e8
simulations <- 100
simulationsDictionaryChoice <- 1
cores <- 3

simulationsDictionary <- data.frame(
  SourceSimulations = c( 
    # Sets of files, as a string, that can be split with ", ".
    # NOTE: Files should be generated with same parameters for comparability.
    # NOTE: Use file.path if you need files in directories.
    paste0(c(), sep = ", ")
  ),
  SourceSimulations2 = c(
    # Sets of files, as a string, that can be split with ", ".
    # NOTE: Files here should be second halves of files in SourceSimulations.
    # NOTE: Use NA if a corresponding second half does not exist.
    # NOTE: Use file.path if you need files in directories.
    paste0(c(), sep = ", ")
  ),
  SourcePools = c(
    # Sets of files, as a string, that can be split with ", ".
    # NOTE: Should correspond one-to-one with the SourceSimulations.
    paste0(c(), sep = ", ")
  ),
  Seeds = c(# runif(1) * 1e8
    69148611
  )
)[simulationsDictionaryChoice, ]

calculationsPlotLong <- FALSE
logarithmicTimeScale <- TRUE

# Libraries: ##################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(patchwork)

library(RMTRCode2)

library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

directory <- '.' # VariousExperiments folder.
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

# Files: ######################################################################
# Convert the simulationsDictionary into proper file pairs.
.sSeed <- simulationsDictionary$Seeds
.sFile1 <- strsplit(simulationsDictionary$SourceSimulations, split = ", ")[[1]]
.sFile2 <- strsplit(simulationsDictionary$SourceSimulations2, split = ", ")[[1]]
.sPool <- strsplit(simulationsDictionary$SourcePools, split = ", ")[[1]]

stopifnot(length(.sFile1) == length(.sFile2),
          length(.sFile1) == length(.sPool))

.s <- foreach::foreach(
  id = seq_along(.sFile1),
  file1 = iterators::iter(.sFile1),
  file2 = iterators::iter(.sFile2),
  fpool = iterators::iter(.sFile2)
) %op% {
  simulation <- loadSimulation(file1, if(file2 != "NA") file2)
  
  # Make sure we clear out any extinct abundances:
  simulation$Abundance[, -1] <- ifelse(
    simulation$Abundance[, -1] <
      simulation$Parameters$EliminationThreshold,
    0, simulation$Abundance[, -1]
  )
  
  load(fpool)
  
  if (exists("pools") && !exists("Pool")) {
    IDNumbers <- sub(basename(file1), pattern = ".RData", replacement = "")
    IDNumbers <- strsplit(IDNumbers, split = "-", fixed = TRUE)[[1]]
    IDNumbers <- IDNumbers[(which(IDNumbers == "Prepared") + 1):length(IDNumbers)]
    IDNumbers <- as.numeric(IDNumbers)
    Pool <- pools[[cases$Parameters[IDNumbers[2]]]][[cases$System[IDNumbers[2]]]]
    # Why 2? 1 == File / Main Case, 2 == row of cases (derived from row of CSV)
    # 3 == The seeds used, 4 == which part of the simulation was saved.
    # Note 3 and 4 are optional if all simulations of a row were saved together.
  }
  
  return(list(
    ID = id,
    File1 = file1,
    File2 = file2,
    FPool = fpool,
    NEnvs = simulation$NumEnvironments,
    NSpec = nrow(Pool),
    Simulation = simulation,
    Pool = Pool
  ))
}

# Setup and Framing: ##########################################################