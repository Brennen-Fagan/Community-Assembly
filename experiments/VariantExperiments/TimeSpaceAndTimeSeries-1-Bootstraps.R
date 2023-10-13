# Introduction: ###############################################################
# As a sequel to TimeSpaceAndTimeSeries.R, can we make the calculation more
# robust by bootstrapping our choices?
# Future Ideas:
#   We can incorporate the file chosen into the bootstrap.

# Things worth adding still:
#    beta diversity at both abundance and species levels (True and Observed).
#    biodiversity intactness index
# And then potentially wrapping it in loops/replicates to get cred. intervals.

# Create a Master seed, which we'll use to generate seeds for each individual
# bootstrap.
# > runif(1) * 1e8
# [1] 36641133
bootstraps <- 100
bootstrapSeed <-
  # 36641133 # Used for 2023-07-06
  # 56901098 # Used for 2023-09-23
  # 83418616 # Used for 2023-09-24
  12854863 # Used for 2023-09-25

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

library(betapart)
library(car)
library(lme4)

clust <- parallel::makeCluster(3)#, outfile = "")
doParallel::registerDoParallel(clust)

# Files: ######################################################################
directory <- '.' # Mutualism folder.
file_result <- "MNA-ExampleExtProp-Result-Env10-Ring-5-1-1-ExtProp1.RData"
file_pool <- "MNA-ExampleOutcome-PoolMats-Env10.RData"
load(file.path(directory, "Data_2023-09-25", file_result))
load(file.path(directory, "Data_2023-09-25", file_pool))
# Make sure we clear out any extinct abundances:
result$Abundance[, -1] <- ifelse(
  result$Abundance[, -1] < result$Parameters$EliminationThreshold,
  0, result$Abundance[, -1]
)

nPatches <- 10 # Fixed in between simulations.
nSpecies <- (ncol(result$Abundance) - 1) / nPatches

stopifnot(nSpecies == floor(nSpecies),
          nSpecies == ceiling(nSpecies),
          nSpecies > 0) # Check if positive integer

source("TimeSpaceAndTimeSeries-Functions.R")

# Setup and Framing: ##########################################################
# Two researchers come upon a set of sites after having been told that the
# local government has designated them for some form of land use change
# (e.g. fertilizer, w/e). They know that they have some amount of time
# before the change is implemented and when the change will be implemented.
# The first researcher decides to examine the sites after the land use change.
# The second researcher decides to do so through the land use change.
# Neither researcher realises that the government doesn't perform the action
# (e.g. lack of funds, w/e).
# How do their results differ between each other and truth?
# Consider the impact of sites monitored (incl. quantity),
#                        monitoring effectiveness (e.g. sampling intensity),
#                        and length of monitoring.

# Time after arrival that intervention takes place.
interventionGap <- 10 * result$ReactionTime

# What should our model of sampling be?
# The typical one sounds like it is be in a place for a period of time,
# Edit: In\^es has found in her experience that space for time at best has
# two temporal sampling points (calculate diff between time points and compare
# between control and experiment/disturbed). At worst, it only compares the
# difference in number of species for a single time point between C & E/D.
# Susan and Jon agreed with In\^es above.
# Set samplingGap to match interventionGap to make for only one time sample.
samplingGap <- result$ReactionTime;

# notate everything that you perceive
samplingFailureRate <- 0.1
# and report back.
# If we assume that sampling is neutral (obviously not) with respect to our
# species then it would presumably be based on abundance.
# So a Poisson dist. based on the total abundance and time period, and then
# assign an identity to each event drawn from the Poisson dist (multinomial).
# For the sake of argument:
samplingPerAbundance <- 1/100
# We'll take a sampling period to be instantaneous.

lastTimeSampleable <- result$Events$Times[length(result$Events$Times)]
lastTimeSampleable <- lastTimeSampleable - samplingGap * 2
firstTimeSampleable <- 1000 # to make sure we're past the simulation burnin.

# The Bootstraps: #############################################################
bootstrapSamples <- foreach::foreach(
  bootstrapID = 1:bootstraps,
  .options.RNG = bootstrapSeed,
  .combine = "rbind",
  .packages = c("dplyr")
) %dorng% {
  # Pick a random patch as control, rest as experiment.
  # is adding 1:5 (or w/e) okay? Yes, we're assuming contiguous patches.
  control <- ((sample.int(nPatches, 1) + 1:(nPatches / 2) ) %% nPatches) + 1
  # print(paste(bootstrapID, ":", toString(control)))
  experiment <- c(1:nPatches)[!c(1:nPatches) %in% control]
  # print(paste(bootstrapID, ":", toString(experiment)))

  # Pick a random start time.
  burnin <- runif(
    n = 1, min = firstTimeSampleable, max = lastTimeSampleable
  ) # When the researchers can arrive, presume steady-state.
  intervention <- burnin + interventionGap # When land use change "occurs".
  # print(paste(bootstrapID, ":", toString(intervention)))

  # Time Series:
  # The premise is that the Time Series researcher is only looking at the
  # plots where the experiment / land use change is taking place, but they
  # instead study the sites before the land use change occurs.
  # They'd be interested in trying to estimate the regional pool / diversity,
  # as well as local diversity and inter-patch diversity.
  sampling_TimeSeries <- data.frame(expand.grid(
    Time = unique(c(seq(from = intervention - samplingGap,
                        by = -samplingGap,
                        to = burnin),
                    seq(from = intervention + samplingGap,
                        by = samplingGap,
                        to = intervention + (intervention - burnin)))),
    Patch = experiment,
    Type = "Time series"
  )) %>% dplyr::arrange(Time)

  # Space for Time:
  # The premise is that the Space for Time researcher is looking between the
  # control and experiment / land use change plots to try to understand how
  # diversity has changed.
  # They'd also be interested in estimating the regional pool / diversity,
  # as well as local diversity and inter-patch diversity.
  sampling_SpaceForTime <- data.frame(expand.grid(
    Time = seq(from = intervention + samplingGap,
               by = samplingGap,
               to = intervention + (intervention - burnin)),
    Patch = c(control, experiment),
    Type = "Space for time"
  ))

  sampling_SpaceForTime <- sampleFromResults(
    sampling = sampling_SpaceForTime,
    resultAbundance = result$Abundance,
    control = control,
    intervention = intervention,
    nSpecies = nSpecies,
    samplingPerAbundance = samplingPerAbundance,
    samplingFailureRate = samplingFailureRate
  )

  sampling_TimeSeries <- sampleFromResults(
    sampling = sampling_TimeSeries,
    resultAbundance = result$Abundance,
    control = control,
    intervention = intervention,
    nSpecies = nSpecies,
    samplingPerAbundance = samplingPerAbundance,
    samplingFailureRate = samplingFailureRate
  )

  sampling_SpaceForTime <- computeSpeciesInControl(sampling_SpaceForTime)
  sampling_TimeSeries <- computeSpeciesInControl(sampling_TimeSeries)
  return(
    dplyr::bind_rows(
      sampling_SpaceForTime,
      sampling_TimeSeries
    ) %>% dplyr::mutate(
      Bootstrap = bootstrapID
    )
  )
}

bootstrapSamples <- bootstrapSamples %>% dplyr::group_by(
  Type, Control, Bootstrap, Patch
) %>% dplyr::arrange(
  Time
) %>% dplyr::mutate(
  TimeSinceStart = Time - min(Time)
) %>% dplyr::ungroup()

parallel::stopCluster(clust)
