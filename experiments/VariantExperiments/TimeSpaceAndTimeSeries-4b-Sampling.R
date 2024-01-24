# Introduction: ###############################################################
# Follows TimeSpaceAndTimeSeries-4a-Simulation.R
# Takes in the Simulations, which have documented:
#  abundances (matrix)
#  intervention times and locations (data.frame)
#  intervention types and magnitudes (list?)
#  event time series (data.frame)
#  sampling random seeds (to maintain a high degree of reproducibility)
# and performs sampling immediately prior to, during, and after intervention.
# We then present the usual set of summary statistics in order to verify
# functionality.
# Together, this file and TimeSpaceAndTimeSeries-1a-Simulation.R should produce
# output similar to TimeSpaceAndTimeSeries-1-Bootstraps and *Intervention.R

# Note that this file should be agnostic to the intervention ideally, outside
# of the information suggested above.

# Parameters: #################################################################
pipeline <- FALSE

if (!exists("directory")) {
  directory <- '.'
}

if (!pipeline) {
  datfolder <-
  # "TSTS_Simulations_1-1-1_2024-01-16"
  # "TSTS_Simulations_1-2-1_2024-01-10"
  # "TSTS_Simulations_1-2-2_2024-01-10"
  # "TSTS_Simulations_1-2-3_2024-01-10"
  # "TSTS_Simulations_1-3-1_2024-01-12"
  # "TSTS_Simulations_1-3-2_2024-01-12"
  # "TSTS_Simulations_1-3-3_2024-01-15"
  # "TSTS_Simulations_2-1-6_2024-01-19"
  # "TSTS_Simulations_2-2-6_2024-01-19"
   "TSTS_Simulations_2-3-6_2024-01-19"
  # "TSTS_Simulations_3-1-7_2024-01-22"
  # "TSTS_Simulations_3-2-7_2024-01-22"
  # "TSTS_Simulations_3-3-7_2024-01-22"
} else {
  # pipeline mode implies that everything should already be defined!!!
  stopifnot(exists("datfolder"))
}

### Quantities: ###############################################################
# When:
samplingQuantity <- 100 # Not guaranteed!
samplingTimeScaleLogarithmic <- TRUE

# How:
samplingFailureRate <- 0.1
samplingPerAbundance <- 1/100
# calculationsPlotLong <- FALSE

# Runs:
samplingRuns <- 20

# Parallelization
cores <- 1

### Automatically load seed matched to each datfolder: ########################
samplingSeed <- switch(
  datfolder,
  "TSTS_Simulations_1-1-1_2024-01-16" = 15351143,
  "TSTS_Simulations_1-2-1_2024-01-10" = 12230846,
  "TSTS_Simulations_1-2-2_2024-01-10" = 17289358,
  "TSTS_Simulations_1-2-3_2024-01-10" = 77022241,
  "TSTS_Simulations_1-3-1_2024-01-12" = 14042213,
  "TSTS_Simulations_1-3-2_2024-01-12" = 40157837,
  "TSTS_Simulations_1-3-3_2024-01-15" = 74027828,
  "TSTS_Simulations_2-1-6_2024-01-19" = 17425645,
  "TSTS_Simulations_2-2-6_2024-01-19" = 51355175,
  "TSTS_Simulations_2-3-6_2024-01-19" = 19450925,
  "TSTS_Simulations_3-1-7_2024-01-22" = 10694873,
  "TSTS_Simulations_3-2-7_2024-01-22" =  3765671,
  "TSTS_Simulations_3-3-7_2024-01-22" = 52689997#,
)

### Fix format change during development: #####################################
defaultTimeSpan <- function(iTS) {
  # Up through 22/01/2024, I had been using either
  # timespan = 20 (1-1-1 through 1-3-3) or
  # timespan = 200 (2-1-6 through 3-3-7).
  # Change of design to include timespan when intervention doesn't have one
  # in order to allow sampling to depend on it means defaults are needed.
  if (is.null(iTS)) {
    return(switch(
      datfolder,
      "TSTS_Simulations_1-1-1_2024-01-16" =  20,
      "TSTS_Simulations_1-2-1_2024-01-10" =  20,
      "TSTS_Simulations_1-2-2_2024-01-10" =  20,
      "TSTS_Simulations_1-2-3_2024-01-10" =  20,
      "TSTS_Simulations_1-3-1_2024-01-12" =  20,
      "TSTS_Simulations_1-3-2_2024-01-12" =  20,
      "TSTS_Simulations_1-3-3_2024-01-15" =  20,
      "TSTS_Simulations_2-1-6_2024-01-19" = 200,
      "TSTS_Simulations_2-2-6_2024-01-19" = 200,
      "TSTS_Simulations_2-3-6_2024-01-19" = 200,
      "TSTS_Simulations_3-1-7_2024-01-22" = 200,
      "TSTS_Simulations_3-2-7_2024-01-22" = 200,
      "TSTS_Simulations_3-3-7_2024-01-22" = 200,
      100 # "Default"
    ))
  } else {
    return(iTS)
  }
}

### Folder Implied Properties: ################################################
datfolder_properties <- strsplit(datfolder, split = "_")
stopifnot(length(datfolder_properties) == 1,
          datfolder_properties[[1]][1] == "TSTS",
          datfolder_properties[[1]][2] == "Simulations")

filename <- file.path(
  datfolder,
  paste0("TSTS_Sampling_", datfolder_properties[[1]][3], ".RData")
)

if(file.exists(filename)) {
  warning(paste(filename, "already exists."))
}

# Load Data: ##################################################################
if (!pipeline) {
  results <- lapply(dir(datfolder, full.names = TRUE,
                        pattern = "Simulation"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1) # Should only be one thing in each file.
    return(get(names))
  })
} else {
  stopifnot(exists("results"))
}

if (as.Date(datfolder_properties[[1]][4], format = "%Y-%m-%d") <=
    as.Date("2024-01-22", format = "%Y-%m-%d")) {
  results <- lapply(results, function(x) {
    # Results collects interventionSimulations.
    x$Ellipsis$Intervention$TimeSpan <-
      defaultTimeSpan(x$Ellipsis$Intervention$TimeSpan)

    if (!"Time" %in% names(x$Ellipsis$Intervention)) {
      x$Ellipsis$Intervention$Time <- x$Ellipsis$TimeIntervention
    }

    return(x)
  })
}

# Libraries: ##################################################################
library(dplyr)

library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Interventions.R"))

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

### Sampling Regime: ##########################################################
createSamplingTimes <- function(timeSpan) {
  # The times (with t = 0 == the intervention time) at which we should sample.
  # Note that this means we will need to add these values to the intervention
  # time (and probably will want to add a manual sampling at the actual start
  # of the intervention run (i.e. pre-intervention)).
  if(samplingTimeScaleLogarithmic) {
    # This version runs from 0 to twice the time span to make sure we have
    # good after intervention coverarge. We take half the timespan as a point
    # of symmetry and scale logarithmically backwards (to 0) and forwards.
    samplingTimes <- sort(unique(c(
      0,
      timeSpan/2 - (
        exp(seq(from = log(1), # We remove the additional 1 at the end.
                to = log(timeSpan/2 - 0 + 1), # We pad 1 for the removal.
                length.out = floor(samplingQuantity / 3))) - 1
        ),
      timeSpan/2 + (
        exp(seq(from = log(1),
                to = log(2*timeSpan - timeSpan/2 + 1),
                length.out = ceiling(samplingQuantity * 2 / 3))) - 1
        ),
      2*timeSpan
    )))
  } else {
    samplingTimes <- seq(from = 0,
                         length.out = samplingQuantity,
                         to = timeSpan * 2)
  }
}

# Perform Sampling: ###########################################################
samples <- foreach::foreach(
  id = 1:samplingRuns,
  r = iterators::iter(results, recycle = TRUE), # File and History (not random!).
  .options.RNG = samplingSeed,
  #.combine = "rbind", # Each result is a list of different objects.
  .packages = c("dplyr")
) %dorng% {

  # Note all on the characteristic time scale.
  samplingTimes <- unique(c(# Guard
    r$Abundance[1, 1], # Pre-intervention, then moment of intervention forward.
    createSamplingTimes(r$Ellipsis$Intervention$TimeSpan) +
    r$Ellipsis$Intervention$Time
  ))

  # In comparison to previous runs, we aren't just looking at a single time
  # in theory versus a single baseline (either previous in time or adjacent in
  # space).
  # Instead, we aren't sure which times we might need, so we need all times
  # and all spaces (with possible exception of the before intervention on
  # non-intervention patches.)
  samplingDataFrame <- data.frame(expand.grid(
    Time = samplingTimes, # True time on characterstic scale.
    Patch = 1:r$NumEnvironments,
    SamplingRun = id
  )) %>% dplyr::arrange(Time) %>% dplyr::mutate(
    PatchType = ifelse(Patch %in% r$Ellipsis$Intervention$Patches,
                       "Experiment", "Control"),
    TimeBase = Time - r$Ellipsis$Intervention$Time, # Time from intervention.
    ParentRun = r$Ellipsis$FullID
  )

  samplingResults <- sampleFromResults2(
    resultAbundance = r$Abundance, # With Time Column
    sampling = samplingDataFrame,
    control = c(1:10)[!1:10 %in% r$Ellipsis$Intervention$Patches],
    intervention = r$Ellipsis$Intervention$Time,
    nSpecies = (ncol(r$Abundance) - 1) / r$NumEnvironments,
    samplingPerAbundance = samplingPerAbundance,
    samplingFailureRate = samplingFailureRate,
    PoolTypes = r$Ellipsis$OriginalRun$PoolTypes
  )
}

save(samples, file = filename)

# Cleanup: ####################################################################
if (exists("clust")) {
  parallel::stopCluster(clust)
}
