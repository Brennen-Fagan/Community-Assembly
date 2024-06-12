# TimeSpaceAndTimeSeries-6e-Sampling.R
# Returning now to the problem of how sampling methods/patterns change
# perceptions of biodiversity.
# See TSTS-4b-Sampling.R for the previous implementation.

# Parameters: #################################################################
datfolder <-
    "TSTS_Simulations_18-1_6-6_2024-05-23"

### Quantities: ###############################################################
# When:
samplingQuantity <- 10 # 100 # Not guaranteed! Number of samples taken by a run
samplingTimeScaleLogarithmic <- TRUE # Logarithmically spread from intervention

# How:
samplingFailureRate <- 0.1 # Change an individual fails to be recorded
samplingPerAbundance <- 1/100 # Converts abundance to rate of encounter
# calculationsPlotLong <- FALSE

# Runs:
samplingRunsPerFile <- 20 # A human readable number. Number of runs.

# Parallelization
cores <- 1


# Files: ######################################################################
# Move all interventions to the front so that we can record the intervention
# times for use with the simulation files, if there are any.
files <- dir(datfolder, full.names = TRUE, pattern = "Simulation")
anySimulations <- length(files)
files <-  c(dir(datfolder, full.names = TRUE, pattern = "Intervention"),
            files)

pool <- load(dir(datfolder, full.names = TRUE, pattern = "Pool"))
pool <- get(pool)

### Automatically load seed matched to each datfolder: ########################
# > runif(1)*1E8
samplingSeed <- switch(
  datfolder,
  "TSTS_Simulations_18-1_6-6_2024-05-23" = 66624262
)

samplingSeedsForFiles <- withRandom(runif(length(files))*1E8, samplingSeed)

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
      10 # "Default"
    ))
  } else {
    return(iTS)
  }
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

# Sampling: ###################################################################
interventionTimes <- c()

samples <- foreach::foreach(
  file = iterators::iter(files),
  seed = iterators::iter(samplingSeedsForFiles)
) %do% {
  fileproperties <- strsplit(file, split = "_", fixed = TRUE)
  fileproperties[[1]][[2]] <- "Sampling"
  filename <- paste0(fileproperties[[1]], collapse = "_")
  if (file.exists(filename)) {
    warning(paste("File", filename, "exists."))
    return(NULL)
  }

  # Load:
  results <- load(file)
  stopifnot(length(names) == 1)
  results <- get(results)

  # Characteristic Time Scale
  if (results$Ellipsis$Timescale == "Simulation") {
    results$Ellipsis$Timescale <- "Characteristic"
    results$Abundance[, 1] <- results$Abundance[, 1] / results$ReactionTime
    results$Events$Times <- results$Events$Times / results$ReactionTime
    if ("TimeIntervention" %in% names(results$Ellipsis$Affinity)) {
      results$Ellipsis$Affinity$TimeIntervention <-
        results$Ellipsis$Affinity$TimeIntervention / results$ReactionTime
    }
  }

  if ("TimeIntervention" %in% names(results$Ellipsis$Affinity)) {
    interventionTimes <-
      c(interventionTimes, results$Ellipsis$Affinity$TimeIntervention)
  } else if (length(unique(interventionTimes) -> uIT) == 1) {
    results$Ellipsis$Affinity$TimeIntervention <- uIT
  } else {return(NULL)}

  # Sampling Time Span
  results$Ellipsis$TimeSpan <-
    defaultTimeSpan(results$Ellipsis$TimeSpan)

  # Make sure to remove spurious abundance if present.
  results$Abundance[, -1] <- ifelse(
    results$Abundance[, -1] > results$Parameters$EliminationThreshold,
    results$Abundance[, -1], 0
  )

  # We want a row pre-intervention, but not too far away from it.
  # In the initial TSTS format, this could reliably be row 1, since the
  # format was to take an existing file, cut it before intervention and make
  # a new file starting at the cut point.
  # In the data format from MNA-Image-ExampleOutcome-Create.R, it's a full
  # file, and so we need to be a bit more careful about how we look around.
  interventionRow <-
    which.max(results$Abundance[, 1] >=
                results$Ellipsis$Affinity$TimeIntervention )

  if (interventionRow - 10 < 0) {
    interventionRow <- 1
  } else {
    interventionRow <- interventionRow - 10
  }

  samplingTimes <- unique(c(# Guard
    # r$Abundance[1, 1], # Pre-intervention, then moment of intervention forward.
    results$Abundance[interventionRow, 1],
    createSamplingTimes(results$Ellipsis$TimeSpan) +
      results$Ellipsis$Affinity$TimeIntervention
  ))

  samplingDataFrame <- data.frame(expand.grid(
    Time = samplingTimes, # True time on characterstic scale.
    Patch = 1:results$NumEnvironments#,
    # SamplingRun = id
  )) %>% dplyr::arrange(Time) %>% dplyr::mutate(
    PatchType = ifelse(Patch %in% r$Ellipsis$Affinity$PatchInterventions,
                       "Experiment", "Control"),
    # Time from intervention.
    TimeBase = Time - results$Ellipsis$Affinity$TimeIntervention,
    ParentRun = results$Ellipsis$ID,
    PatchAffinity = ifelse(
      PatchType == "Experiment",
      results$Ellipsis$Affinity$PatchAffinitiesIntervention,
      result$Ellipsis$Affinity$PatchAffinitiesOld)
  )

  samplingSeedsForRuns <-
    withRandom(runif(samplingRunsPerFile) * 1E8,
               samplingSeedsForFiles)

  samplingResults <- foreach::foreach(
    id = 1:samplingRunsPerFile,
    sd = iterators::iter(samplingSeedsForRuns),
    .packages = "dplyr"
  ) %dopar% {
    withRandom(sampleFromResults2(
      resultAbundance = results$Abundance,
      sampling = samplingDataFrame %>% dplyr::mutate(SamplingRun = id),
      control = c(1:results$NumEnvironments)[
        ! 1:results$NumEnvironments %in%
          results$Ellipsis$Affinity$PatchInterventions
        ],
      intervention = results$Ellipsis$Affinity$TimeIntervention,
      nSpecies = (ncol(results$Abundance) - 1) / results$NumEnvironments,
      samplingPerAbundance = samplingPerAbundance,
      samplingFailureRate = samplingFailureRate,
      PoolTypes = table(pool$Affinity)
    ) %>% dplyr::mutate(
      ParentRun = results$Ellipsis$ID,
      Seed = sd
      ),
    sd)
  }

  save(samplingResults, file = filename)

  return(samplingResults)
}

# Cleanup: ####################################################################
if (exists("clust")) {
  parallel::stopCluster(clust)
}
