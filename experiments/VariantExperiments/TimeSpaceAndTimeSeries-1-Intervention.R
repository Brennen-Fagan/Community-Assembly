# Introduction: ###############################################################
# As a sequel to TimeSpaceAndTimeSeries-1-Bootstrap.R,
# we are now introducing an actual (file-substitution) intervention.
# Please see the previous file for some design choices.

# Create a Master seed, which we'll use to generate seeds for each individual
# bootstrap.
# > runif(1) * 1e8
# [1] 97870743
bootstraps <- 100
bootstrapSeed <-
  # 97870743 # used for 2023-09-25, 5 -> Hik6_2023-03-04, 56-6
  21139057 # used for 2023-09-25, 5 -> Inf

calculationsPlotLong <- FALSE
logarithmicTimeScale <- TRUE
relabelPoolIntervention <- FALSE
timeInterventionRandom <- FALSE
patchInterventionRandom <- FALSE
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

library(vegan) # vegdist

library(betapart)
library(car)
library(lme4)

directory <- '.' # VariousExperiments folder.
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))

clust <- parallel::makeCluster(3, outfile = "")
doParallel::registerDoParallel(clust)

# Files: ######################################################################

# Base
fileBaseFolder <- "Data_2023-09-25"
fileBaseResult <- "MNA-ExampleExtProp-Result-Env10-Ring-5-1-1-ExtProp1.RData"
fileBaseResult2 <- NULL
fileBasePool <- "MNA-ExampleOutcome-PoolMats-Env10.RData"

# Intervention
fileInterventionFolder <- "Data_2023-09-25"
fileInterventionResult <- "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData"
fileInterventionResult2 <- NULL
fileInterventionPool <- "MNA-ExampleOutcome-PoolMats-Env10.RData"

resultBase <- loadSimulation(
  file.path(directory, fileBaseFolder, fileBaseResult),
  if(!is.null(fileBaseResult2))
    file.path(directory, fileBaseFolder, fileBaseResult2)
)
resultIntervention <- loadSimulation(
  file.path(directory, fileInterventionFolder, fileInterventionResult),
  if(!is.null(fileInterventionResult2))
    file.path(directory, fileInterventionFolder, fileInterventionResult2)
)

load(file.path(directory, fileBaseFolder, fileBasePool))
if (exists("pools") && !exists("Pool")) {
  IDNumbers <- sub(fileBaseResult, pattern = ".RData", replacement = "")
  IDNumbers <- strsplit(IDNumbers, split = "-", fixed = TRUE)[[1]]
  IDNumbers <- IDNumbers[(which(IDNumbers == "Prepared") + 1):length(IDNumbers)]
  IDNumbers <- as.numeric(IDNumbers)
  Pool <- pools[[cases$Parameters[IDNumbers[2]]]][[cases$System[IDNumbers[2]]]]
  # Why 2? 1 == File / Main Case, 2 == row of cases (derived from row of CSV)
  # 3 == The seeds used, 4 == which part of the simulation was saved.
  # Note 3 and 4 are optional if all simulations of a row were saved together.
}
PoolBase <- Pool; rm(Pool)

load(file.path(directory, fileInterventionFolder, fileInterventionPool))
if (exists("pools") && !exists("Pool")) {
  IDNumbers <- sub(fileInterventionResult, pattern = ".RData", replacement = "")
  IDNumbers <- strsplit(IDNumbers, split = "-", fixed = TRUE)[[1]]
  IDNumbers <- IDNumbers[(which(IDNumbers == "Prepared") + 1):length(IDNumbers)]
  IDNumbers <- as.numeric(IDNumbers)
  Pool <- pools[[cases$Parameters[IDNumbers[2]]]][[cases$System[IDNumbers[2]]]]
  # Why 2? 1 == File / Main Case, 2 == row of cases (derived from row of CSV)
  # 3 == The seeds used, 4 == which part of the simulation was saved.
  # Note 3 and 4 are optional if all simulations of a row were saved together.
}
PoolIntervention <- Pool; rm(Pool)

# Make sure we clear out any extinct abundances:
resultBase$Abundance[, -1] <- ifelse(
  resultBase$Abundance[, -1] <
    resultBase$Parameters$EliminationThreshold,
  0, resultBase$Abundance[, -1]
)
resultIntervention$Abundance[, -1] <- ifelse(
  resultIntervention$Abundance[, -1] <
    resultIntervention$Parameters$EliminationThreshold,
  0, resultIntervention$Abundance[, -1]
)

nPatches <- resultBase$NumEnvironments
# Why not Intervention? Because the Base patches will be 1-to-1 replaced.

if (relabelPoolIntervention) {
  # Pool Problems:
  # NOTE: We use "which" in the sampling scheme.
  #       We SHOULD insert the Pool$ID[which(...)] instead, but
  #       all of our pools label species 1:nrow(Pool) anyways.
  # Brute force labeling and relabeling.
  if (!"ID" %in% colnames(PoolIntervention)) {
    PoolIntervention$ID <- 1:nrow(PoolIntervention)
  }
  flagAdditions <- 0 # Used to update Events later.
  while (any(PoolIntervention$ID %in% PoolBase$ID)) {
    flagAdditions <- flagAdditions + 1
    PoolIntervention$ID <- length(PoolBase$ID) + PoolIntervention$ID
  }

  # Abundance Problems:
  # Since we perform our sampling based on column position, we need to
  # give the different species different column positions.
  # => Inserting empty (rep(0)) columns.
  # Easiest method is probably to
  #   separate matrices into environments,
  #   cbind to each environment the appropriate new columns,
  #   then recombine the environments as before.

  # Note BASE first, then INTERVENTION
  resultBase$Abundance <- cbind(
    resultBase$Abundance[, 1], # Time Column
    do.call(cbind, lapply(
      1:resultBase$NumEnvironments,
      function(i, ab, nspec) {
        # Time Col. + Past Env. Cols + This Env.'s Cols
        temp <- ab[, 1 + (i - 1) * nspec + (1:nspec)]

        return(cbind(
          temp, matrix(0, nrow = nrow(temp), ncol = nrow(PoolIntervention))
        ))
      },
      ab = resultBase$Abundance,
      nspec = nrow(PoolBase)
    ))
  )

  # Note reverse ordering!
  resultIntervention$Abundance <- cbind(
    resultIntervention$Abundance[, 1], # Time Column
    do.call(cbind, lapply(
      1:resultIntervention$NumEnvironments,
      function(i, ab, nspec) {
        # Time Col. + Past Env. Cols + This Env.'s Cols
        temp <- ab[, 1 + (i - 1) * nspec + (1:nspec)]

        return(cbind(
          matrix(0, nrow = nrow(temp), ncol = nrow(PoolBase)), temp
        ))
      },
      ab = resultIntervention$Abundance,
      nspec = nrow(PoolIntervention)
    ))
  )

  # Event Problems:
  resultIntervention$Events$Species <-
    resultIntervention$Events$Species + flagAdditions * length(PoolBase$ID)

  nSpecies <- nrow(PoolBase) + nrow(PoolIntervention)
} else if (
  length(setdiff(PoolBase$ID, PoolIntervention$ID)) > 0 ||
  length(setdiff(PoolIntervention$ID, PoolBase$ID)) > 0 ) {
  stop(
    "relabelPoolIntervention == FALSE but different Pool IDs not implemented"
  )
  # Pool Problems: We're not relabeling, but one pool is obv. different.
  # Nothing to do here.

  # Abundance Problems: We are, however, missing some columns.
  # The cols to add are non-trivial!
  # We will have to trust that the different IDs AND
  # the orders they appear in are correct!

  MissingFromIntervention <- setdiff(PoolBase$ID, PoolIntervention$ID)
  MissingFromBase <- setdiff(PoolBase$ID, PoolIntervention$ID)

  if (length(MissingFromIntervention)) {
    # Need to add columns, in the right places!, to Intervention Abundance.

  }

  # Event Problems: Differences in species in pool isn't a problem for events.
  # Nothing to do here.


  nSpecies <- length(union(PoolBase, PoolIntervention))
} else {
  nSpecies <- nrow(PoolBase) # Should be the same as PoolIntervention.
}

stopifnot(nSpecies == floor(nSpecies),
          nSpecies == ceiling(nSpecies),
          nSpecies > 0) # Check if positive integer

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

# Interesting problem: we need everything to be on the same time scales.
# We're already using the reaction times as our best proxy/control of the time
# scale, so we'll reformat all of the times using those.
resultBase$Abundance[, 1] <-
  resultBase$Abundance[, 1] / resultBase$ReactionTime
resultBase$Events$Times <-
  resultBase$Events$Times / resultBase$ReactionTime
resultIntervention$Abundance[, 1] <-
  resultIntervention$Abundance[, 1] / resultIntervention$ReactionTime
resultIntervention$Events$Times <-
  resultIntervention$Events$Times / resultIntervention$ReactionTime

# Time after arrival that intervention takes place.
# Since same time scale, opting for the same value here.
timeGap <- 100

# What should our model of sampling be?
# The typical one sounds like it is be in a place for a period of time,
# Edit: In\^es has found in her experience that space for time at best has
# two temporal sampling points (calculate diff between time points and compare
# between control and experiment/disturbed). At worst, it only compares the
# difference in number of species for a single time point between C & E/D.
# Susan and Jon agreed with In\^es above.
# Set samplingGap to match interventionGap to make for only one time sample.
samplingGap <- 1

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

lastTimeSampleable <- resultBase$Events$Times[
  length(resultBase$Events$Times)
  ]
lastTimeSampleable <-
  lastTimeSampleable - timeGap * 2

lastTimeSampleableAlternate <- resultIntervention$Events$Times[
  length(resultIntervention$Events$Times)
  ]
lastTimeSampleableAlternate <-
  lastTimeSampleableAlternate - timeGap * 2

# to make sure we're past the simulation burnin.
firstTimeSampleableBase <- 125
firstTimeSampleableIntervention <- 125

# Enforce the first sampleable row for the intervention.
interventionFirstRowSampleable <-
  which.max(resultIntervention$Abundance[, 1] > firstTimeSampleableIntervention)

if(logarithmicTimeScale) {
  # This version is symmetric on the log scale, centred on the sampling gap,
  # and ends at the time gap. Number of sampling times not guaranteed.
  timediffs <- unique(exp(c(
    seq(from = -log(timeGap),
        to = log(samplingGap),
        length.out = floor(timeGap/samplingGap/2)),
    seq(from = log(samplingGap),
        to = log(timeGap),
        length.out = ceiling(timeGap/samplingGap/2))
  )))
  timediffs <- timediffs[timediffs > min(c(
    diff(resultBase$Abundance[, 1]),
    diff(resultIntervention$Abundance[, 1])
    ))]
} else {
  timediffs <- seq(from = samplingGap,
                   by = samplingGap,
                   to = timeGap)
}

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
    n = 1, min = firstTimeSampleableBase, max = lastTimeSampleable
  ) # When the researchers can arrive, presume steady-state.
  timeSwitch <- burnin + timeGap # When land use change "occurs".
  # print(paste(bootstrapID, ":", toString(intervention)))

  timeAlternate <- if (timeInterventionRandom) {
    runif(
      n = 1, min = firstTimeSampleableIntervention,
      max = lastTimeSampleableAlternate
    )
  } else {
    timeSwitch
  }

  patchAlternate <- if(patchInterventionRandom) {
    sample(resultIntervention$NumEnvironments)
  } else {
    1:resultBase$NumEnvironments
  }

  # Time Series:
  # The premise is that the Time Series researcher is only looking at the
  # plots where the experiment / land use change is taking place, but they
  # instead study the sites before the land use change occurs.
  # They'd be interested in trying to estimate the regional pool / diversity,
  # as well as local diversity and inter-patch diversity.
  sampling_TimeSeries <- data.frame(expand.grid(
    TimeBase = unique(c(timeSwitch - timediffs,
                        timeSwitch + timediffs)),
    Patch = experiment,
    Type = "Time series"
  )) %>% dplyr::arrange(TimeBase) %>% dplyr::mutate(
    TimeIntervention = TimeBase - timeSwitch + timeAlternate,
    PatchIntervention = patchAlternate[Patch]
  )

  # Space for Time:
  # The premise is that the Space for Time researcher is looking between the
  # control and experiment / land use change plots to try to understand how
  # diversity has changed.
  # They'd also be interested in estimating the regional pool / diversity,
  # as well as local diversity and inter-patch diversity.
  sampling_SpaceForTime <- data.frame(expand.grid(
    TimeBase = timeSwitch + timediffs,
    Patch = c(control, experiment),
    Type = "Space for time"
  )) %>% dplyr::arrange(TimeBase) %>% dplyr::mutate(
    TimeIntervention = TimeBase - timeSwitch + timeAlternate,
    PatchIntervention = patchAlternate[Patch]
  )

  # print(paste(bootstrapID, "made it."))

  sampling_SpaceForTime <- sampleFromResultsIntervention(
    sampling = sampling_SpaceForTime,
    baseAbundance = resultBase$Abundance,
    interventionAbundance = resultIntervention$Abundance,
    control = control,
    interventionTime = timeSwitch,
    nSpecies = nSpecies,
    samplingPerAbundance = samplingPerAbundance,
    samplingFailureRate = samplingFailureRate,
    Pool = rbind(PoolBase, PoolIntervention)
  )

  # print(paste(bootstrapID, "made it."))

  sampling_TimeSeries <- sampleFromResultsIntervention(
    sampling = sampling_TimeSeries,
    baseAbundance = resultBase$Abundance,
    interventionAbundance = resultIntervention$Abundance,
    control = control,
    interventionTime = timeSwitch,
    nSpecies = nSpecies,
    samplingPerAbundance = samplingPerAbundance,
    samplingFailureRate = samplingFailureRate,
    Pool = rbind(PoolBase, PoolIntervention)
  )

  # print(paste(bootstrapID, "made it."))

  sampling_SpaceForTime <- computeSpeciesInControl(
    sampling_SpaceForTime, Time = "TimeBase")
  sampling_TimeSeries <- computeSpeciesInControl(
    sampling_TimeSeries, Time = "TimeBase")

  print(paste(bootstrapID, "made it."))

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
  TimeBase
) %>% dplyr::mutate(
  TimeSinceStart = TimeBase - min(TimeBase)
) %>% dplyr::ungroup()

parallel::stopCluster(clust)

# Plotting: ###################################################################
### Plot 0: Sense Checking: ###################################################
