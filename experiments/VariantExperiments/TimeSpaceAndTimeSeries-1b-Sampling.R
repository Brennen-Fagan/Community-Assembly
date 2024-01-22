# Introduction: ###############################################################
# Follows TimeSpaceAndTimeSeries-1a-Simulation.R
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
    "TSTS_Simulations_1-1-1_2024-01-16"
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
      "TSTS_Simulations_3-3-7_2024-01-22" = 200#,
    ))
  } else {
    return(iTS)
  }
}

# Load Data: ##################################################################
if (!pipeline) {
  results <- lapply(dir(datfolder, full.names = TRUE), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1) # Should only be one thing in each file.
    return(get(names))
  })
} else {
  stopifnot(exists("results"))
}

# Libraries: ##################################################################
library(dplyr)
library(ggplot2)
library(patchwork)

source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Interventions.R"))

### Sampling Regime: ##########################################################
# The times (with t = 0 == the intervention time) at which we should sample.
if(samplingTimeScaleLogarithmic) {
  # This version is symmetric on the log scale, centred on 1 time unit,
  # and ends at the time gap. Number of sampling times not guaranteed.
  # The centre is chosen for its relevance to the characteristic time scale.
  samplingTimes <- c(0, unique(exp(c(
    seq(from = log(1),
        to = -log(samplingMaxTime),
        length.out = floor(samplingQuantity/2)),
    seq(from = log(1),
        to = log(samplingMaxTime),
        length.out = ceiling(samplingQuantity/2))
  ))))
} else {
  samplingTimes <- seq(from = 0,
                       by = samplingMaxTime/samplingQuantity,
                       to = samplingMaxTime)
}
