# Introduction: ###############################################################
# As a sequel to TimeSpaceAndTimeSeries-1-Bootstrap.R and -1-Intervention.R,
# we are now introducing a simulation based intervention.
# Please see the previous files for some design choices, although we aim to
# improve design at each stage.

# Create a Master seed, which we'll use to generate simulation seeds.
# runif(1) * 1e8
simulations <- 100
simulationsDictionaryChoice <- 1
interventionChoice <-
  1. # The change does not occur, but the scientists don't know that!
# 2. # The change occurs, is discontinuous, and local, e.g. forest fire.
# 3. # The change is slow and local, e.g. gradual deforestation.
# 4. # The change occurs on all patches, e.g. climate change.
# 5. # Two changes, one following 2, the other 4.
# 6. # Two changes, one following 3, the other 4.

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

  # Interesting problem: we need everything to be on the same time scales.
  # We already use the reaction times as our best proxy/control of the time
  # scale, so we'll reformat all of the times using those.
  simulation$Abundance[, 1] <-
    simulation$Abundance[, 1] / simulation$ReactionTime
  simulation$Events$Times <-
    simulation$Events$Times / simulation$ReactionTime

  load(fpool)

  # Some files we might be interested in have multiple pools bundled.
  # This code uses the file names to extract the correct one.
  if (exists("pools") && !exists("Pool")) {
    IDNumbers <- sub(basename(file1), pattern = ".RData", replacement = "")
    IDNumbers <- strsplit(IDNumbers, split = "-", fixed = TRUE)[[1]]
    IDNumbers <- IDNumbers[
      (which(IDNumbers == "Prepared") + 1):length(IDNumbers)
      ]
    IDNumbers <- as.numeric(IDNumbers)
    Pool <- pools[[cases$Parameters[IDNumbers[2]]]][[
      cases$System[IDNumbers[2]]
      ]]
    # Why 2? 1 == File / Main Case, 2 == row of cases (derived from CSV row)
    # 3 == The seeds used, 4 == which part of the simulation was saved.
    # Note 3 & 4 are optional if all simulations of a row were saved together.
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
# (Supposing omniscience == know all information about the present, !future:)
# Two researchers, an omniscient recorder, and a dimension hopping omniscient
# recorder come upon a set of sites after having been told that the
# local government has designated them for some form of land use change
# (e.g. fertilizer, w/e). The first researcher arrives with some time
# before the change is implemented and will compare the sites before and after
# the land use change. The second researcher is too late and needs to compare
# instead with surrounding sites that did not undergo the land use change.
# Meanwhile the omniscient recorder knows the actual changes (no sampling), and
# the dimension hopper *also* knows what would have happened if the land use
# change did not, in fact, occur.
#
# Furthermore, suppose that there are 6 universes, each with a separate change:
#   1. The change does not occur, but the scientists don't know that!
#   2. The change occurs, is discontinuous, and local, e.g. forest fire.
#   3. The change is slow and local, e.g. gradual deforestation.
#   4. The change occurs on all patches, e.g. climate change.
#   5. Two changes, one following 2, the other 4.
#   6. Two changes, one following 3, the other 4.
#
# How do the four observers' results differ between each other in each universe?
# Consider the impact of sites monitored (incl. quantity),
#                        monitoring effectiveness (e.g. sampling intensity),
#                        and length of monitoring.
#
# For comparability, interventions are all considered to start at the same time,
# but the samples we compare may be midway through or after interventions, esp.
# when comparing slow to fast or immediate interventions.

# Choose a sigmoidal interpolation between the matrices.
# Problem: we want no deviation at 0 or t = timespan, but not too much in mid.
# Crit. Value = first numb. in tanh. 4 => 2% dev at edge, ~75% done in mid. 50%.
# Require: Output is a function.
# Require: Output: t = timespan / 2 => out = (mat1 + mat2) / 2.
interpolateMatrices <- function(matrix1, matrix2, timespan) {
  stopifnot(dim(matrix1) == dim(matrix2))
  force(matrix1);force(matrix2);force(timespan)
  function(t, ...) {
    matrix1 + (matrix2 - matrix1) * (tanh(4 * (t/timespan - 0.5)) + 1) / 2
  }
}
switchMatrices <- function(matrix1, matrix2, time) {
  stopifnot(dim(matrix1) == dim(matrix2))
  force(matrix1);force(matrix2);force(timespan)
  function(t, ...) {
    if (t < time) {
      matrix1
    } else if (t > time) {
      matrix2
    } else {
      (matrix1 + matrix2)/2
    }
  }
}

# If we've coded things correctly, we should be able to re-run our original
# RMTRCode2::MultipleNumericalAssembly_Dispersal(...) simulation.
# Recall that this function
#   1) can start from arbitrary initial populations and times,
#   2) does not need the interaction matrices themselves,
#   3) only needs the pool size to not change, the rest of the pool can,
#   4) can take the fixed events list, and
#   5) can take any PerCapitaDynamics that accepts (t, y, parms) arguments.
# Downsides:
#   1) we'll have to fix a characteristic rate. Pick the most frequent sampling.
#   2) The DispersalMatrix will need to remain fixed.

# In order to charge ahead then, we need to
#   1) pick a base file and history,
#   2) reformat it to meet the MNAD criteria,
#   3) create a PerCapitaDynamics function to pass to MNAD.
# The last of these can have some preparation ahead of time.


