# Introduction: ###############################################################
# As a sequel to TimeSpaceAndTimeSeries-1a-Bootstrap.R and 
# TimeSpaceAndTimeSeries-2a-Intervention.R, we are now introducing an actual 
# (file-substitution) intervention that is not limited to a single file pair.
# Please see the previous files for some design choices.

# Create a Master seed, which we'll use to generate seeds for each individual
# bootstrap.
# > runif(1) * 1e8
# [1] 97870743
bootstraps <- 100
bootstrapSeed <-
  # 97870743 # used for 2023-09-25, 5 -> Hik6_2023-03-04, 56-6
  21139057 # used for 2023-09-25, 5 -> Inf

directory <- '.' # VariousExperiments folder.

calculationsPlotLong <- FALSE

# Parameters/Rules: ###########################################################
logarithmicTimeScale <- TRUE # Sampling on log spacings or linear spacings.
relabelPoolIntervention <- FALSE # Second file's pool is only new species?
timeInterventionRandom <- FALSE # Intervention moves to new file at new time?
patchInterventionRandom <- FALSE # Intervention moves to new file at new patch?

# Libraries: ##################################################################

# Data handling:
library(dplyr)
library(tidyr)

# Core personal functions
library(RMTRCode2)

# Parallelisation
library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

# Additional functions used exclusively for this analysis.
# (Should be folded into a new package or into a new function file.)
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))

# Parallelisation: ############################################################
clust <- parallel::makeCluster(3, outfile = "")
doParallel::registerDoParallel(clust)

# Files: ######################################################################
# One problem with this approach is that we either need to load up all results
# and pools in the main thread and then copy them from their to the workers, or
# we need to load the appropriate ones directly to the workers (and deal with
# reloading them huge amounts of times).
# Part of this depends on implementation (thinking deep versus shallow copy).