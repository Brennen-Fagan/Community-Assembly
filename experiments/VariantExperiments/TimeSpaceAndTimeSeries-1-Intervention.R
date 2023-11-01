# Introduction: ###############################################################
# As a sequel to TimeSpaceAndTimeSeries-1-Bootstrap.R,
# we are now introducing an actual (file-substitution) intervention.
# Please see the previous file for some design choices.

# Create a Master seed, which we'll use to generate seeds for each individual
# bootstrap.
# > runif(1) * 1e8
# [1] 97870743
bootstraps <- 100
bootstrapSeed <-97870743

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

clust <- parallel::makeCluster(3)#, outfile = "")
doParallel::registerDoParallel(clust)

# Files: ######################################################################

# Base
fileBaseFolder <- "Data_Hik6_2023-03-04"
fileBaseResult <- "MNA-Hik6A-Cases-Prepared-1-56-6-1.RData"
fileBaseResult2 <- "MNA-Hik6A-Cases-Prepared-1-56-6-2.RData"
fileBasePool <- "MNA-Hik6A-Cases-Prepared.RData"

# Intervention
fileInterventionFolder <- "Data_Hik6_2023-03-04"
fileInterventionResult <- "MNA-Hik6A-Cases-Prepared-1-56-6-1.RData"
fileInterventionResult2 <- "MNA-Hik6A-Cases-Prepared-1-56-6-2.RData"
fileInterventionPool <- "MNA-Hik6A-Cases-Prepared.RData"

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
  resultBase$Abundance <- dplyr::bind_cols(
    resultBase$Abundance[, 1], # Time Column
    lapply(
      1:resultBase$NumEnvironments,
      function(i, ab, nspec) {
        # Time Col. + Past Env. Cols + This Env.'s Cols
        temp <- ab[, 1 + (i - 1) * nspec + (1:nspec)]

        return(dplyr::bind_cols(
          temp, matrix(0, nrow = nrow(temp), ncol = nrow(PoolIntervention))
        ))
      },
      ab = resultBase$Abundance,
      nspec = nrow(PoolBase)
    )
  )

  # Note reverse ordering!
  resultIntervention$Abundance <- dplyr::bind_cols(
    resultIntervention$Abundance[, 1], # Time Column
    lapply(
      1:resultIntervention$NumEnvironments,
      function(i, ab, nspec) {
        # Time Col. + Past Env. Cols + This Env.'s Cols
        temp <- ab[, 1 + (i - 1) * nspec + (1:nspec)]

        return(dplyr::bind_cols(
          matrix(0, nrow = nrow(temp), ncol = nrow(PoolBase)), temp
        ))
      },
      ab = resultIntervention$Abundance,
      nspec = nrow(PoolIntervention)
    )
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

# Time after arrival that intervention takes place.
interventionGap <- 100 * result$ReactionTime

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

lastTimeSampleable <- resultBase$Events$Times[
  length(resultBase$Events$Times)
  ]
lastTimeSampleable <- lastTimeSampleable - interventionGap * 2

lastTimeSampleableAlternate <- resultIntervention$Events$Times[
  length(resultIntervention$Events$Times)
  ]
lastTimeSampleableAlternate <- lastTimeSampleableAlternate - interventionGap * 2

firstTimeSampleable <- 1000 # to make sure we're past the simulation burnin.

# Enforce the first sampleable row for the intervention.
interventionFirstRowSampleable <-
  which.max(resultIntervention$Abundance[, 1] > firstTimeSampleable)

if(logarithmicTimeScale) {
  #timediffs <- exp(seq(from = -log(interventionGap/samplingGap),
  #                 to = log(interventionGap/samplingGap),
  #                 length.out = round(interventionGap/samplingGap)))
  #timediffs <- exp(seq(from = -log(interventionGap),
  #                     to = log(interventionGap),
  #                     length.out = round(interventionGap/samplingGap)))
  # This version is symmetric on the log scale, centred on the sampling gap,
  # and ends at the intervention gap. Number of sampling times not guaranteed.
  timediffs <- unique(exp(c(
    seq(from = -log(interventionGap),
        to = log(samplingGap),
        length.out = floor(interventionGap/samplingGap/2)),
    seq(from = log(samplingGap),
        to = log(interventionGap),
        length.out = ceiling(interventionGap/samplingGap/2))
  )))
  timediffs <- timediffs[timediffs > min(diff(result$Abundance[, 1]))]
} else {
  timediffs <- seq(from = samplingGap,
                   by = samplingGap,
                   to = interventionGap)
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
    n = 1, min = firstTimeSampleable, max = lastTimeSampleable
  ) # When the researchers can arrive, presume steady-state.
  timeSwitch <- burnin + interventionGap # When land use change "occurs".
  # print(paste(bootstrapID, ":", toString(intervention)))

  timeAlternate <- if (timeInterventionRandom) {
    runif(
      n = 1, min = firstTimeSampleable, max = lastTimeSampleableAlternate
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
  )) %>% dplyr::arrange(Time) %>% dplyr::mutate(
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
    Time = interventionTime + timediffs,
    Patch = c(control, experiment),
    Type = "Space for time"
  )) %>% dplyr::arrange(Time) %>% dplyr::mutate(
    TimeIntervention = TimeBase - timeSwitch + timeAlternate,
    PatchIntervention = patchAlternate[Patch]
  )

  sampling_SpaceForTime <- sampleFromResultsIntervention(
    sampling = sampling_SpaceForTime,
    baseAbundance = resultBase$Abundance,
    interventionAbundance = resultIntervention$Abundance,
    control = control,
    interventionTime = timeSwitch,
    nSpecies = nSpecies,
    samplingPerAbundance = samplingPerAbundance,
    samplingFailureRate = samplingFailureRate
  )

  sampling_TimeSeries <- sampleFromResultsIntervention(
    sampling = sampling_TimeSeries,
    baseAbundance = resultBase$Abundance,
    interventionAbundance = resultIntervention$Abundance,
    control = control,
    interventionTime = interventionTime,
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

# Plotting: ###################################################################
### Plot 0: Sense Checking: ###################################################
# Richness Comparison
plot_0_Richness <- ggplot2::ggplot(
  bootstrapSamples %>% dplyr::mutate(
    TrueRichness = unlist(lapply(
      strsplit(SamplingNonZeroSpecies, split = ", ", fixed = TRUE), length
    ))
  ),
  ggplot2::aes(x = Time/result$ReactionTime, y = TrueRichness, color = Patch, group = Patch)
) + ggplot2::geom_line(
) + ggplot2::geom_line(
  ggplot2::aes(y = SamplingAlpha * 3),
  alpha = 0.25,
  linetype = 2
) + ggplot2::scale_y_continuous(
  sec.axis = ggplot2::sec_axis(~ . / 3, name = "Richness Observed")
) + ggplot2::facet_wrap(
  . ~ Patch
) + ggplot2::labs(
  x = "Time (Units: Characteristic Times)",
  y = "True Richness"
)

# But correlation remains low? (but significant) (For file 5-1-1)
# with(bootstrapSamples %>% dplyr::mutate(
#   TrueRichness = unlist(lapply(
#     strsplit(SamplingNonZeroSpecies, split = ", ", fixed = TRUE), length
#   ))
# ), cor.test(TrueRichness, SamplingAlpha))

plot_0_Abundance <- ggplot2::ggplot(
  bootstrapSamples,
  ggplot2::aes(x = Time/result$ReactionTime, y = SamplingAbundance,
               color = Patch, group = Patch)
) + ggplot2::geom_line(
) + ggplot2::geom_line(
  ggplot2::aes(y = SamplingObserved * 100),
  alpha = 0.25,
  linetype = 2
) + ggplot2::scale_y_continuous(
  name = "True Total Abundance",
  sec.axis = ggplot2::sec_axis(~ . / 100, name = "Abundance Observed")
) + ggplot2::facet_wrap(. ~ Patch
) + ggplot2::labs(
  x = "Time (Units: Characteristic Times)"
)

# Abundance is (unsurprisingly) much more correlated. (For file 5-1-1)
# with(bootstrapSamples, cor.test(SamplingAbundance, SamplingObserved))

if (calculationsPlotLong) {
  plot_0_AbundanceSpecies <- lapply(
    1:result$NumEnvironments,
    function(i) RMTRCode2::LawMorton1996_PlotAbundance(
      result$Abundance[
        result$Abundance[, 1] > min(bootstrapSamples$Time) &
          result$Abundance[, 1] < max(bootstrapSamples$Time)
        , c(1, (i-1) * 200 + 1:200 + 1)], guides = FALSE
    ) + ggplot2::scale_y_log10()
  ) %>% patchwork::wrap_plots()
}

# bootstrapSamplesOccupancy <- bootstrapSamples %>% dplyr::group_by(
#   Bootstrap, Time, Type, Control, TimeActualRow, TimeActual, TimeSinceStart
# ) %>% dplyr::group_map(
#   function(x, y) {
#     # Per Bootstrap, Time, Patch Set (Control/Experiment)
#     # We are reproducing our Figure 3a in principle.
#     # Note that that would mean No. of Bootstrap plots!
#     # It doesn't make sense to average these on a per species basis.
#     # It might make sense, however, to average on the dominant niche
#     # (according to the pool) though.
#     # How do we minimise the number of plots? We can reduce to 100 + 2
#     # by having Fig 3a equivalent have nowhere present as NA, and present
#     # as -5 : 5 (all Control : all Experiment) and the 2 as richness plots.
#
#     # Second thoughts: maybe we only need the true occupancy.
#
#     TruePresent <- unique(unlist(
#       strsplit(x$SamplingNonZeroSpecies, split = ", ")
#       ))
#
#   }
# )

plot_0_OccupancyTrue <- (
  speciesAbundances <- RMTRCode2::Calculate_Species(
    result
  )
) %>% dplyr::group_by(
  Time, Species
) %>% dplyr::summarise(
  Count = dplyr::n()
) %>% dplyr::left_join(
  Pool %>% dplyr::arrange(Size) %>% dplyr::mutate(SizeID = 1:nrow(Pool)),
  by = c("Species" = "ID")
) %>% ggplot2::ggplot(
  ggplot2::aes(x = Time/result$ReactionTime, y = SizeID, color = Count)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1, limits = c(1, 10)
) + ggplot2::geom_hline(
  yintercept = nrow(Pool %>% dplyr::filter(Type == Type[1])) + 0.5,
  color = "red"
) + ggplot2::labs(
  x = "Time (Units: Characteristic Times)",
  y = "Species by Size"
) + ggplot2::theme_bw()


(speciesAbundanceStats <- speciesAbundances %>% dplyr::group_by(
  Time, Species
) %>% dplyr::summarise(
  Count = dplyr::n(),
  AvgAbundance = mean(Abundance),
  AvgLogAbundance = mean(log(Abundance))
))

plot_0_AbundanceAvg <- speciesAbundanceStats %>% ggplot2::ggplot(
  ggplot2::aes(x = Time/result$ReactionTime, y = AvgAbundance,
               color = as.character(Species), group = Species)
) + ggplot2::geom_line(
) + ggplot2::labs(
  x = "Time (Units: Characteristic Times)",
  y = "Average Abundance"
) + ggplot2::facet_wrap(
  . ~ Count
) + ggplot2::theme_bw() + ggplot2::scale_y_log10(
) + ggplot2::guides(color = "none")

plot_0_AbundanceGeoAvg <- speciesAbundanceStats %>% ggplot2::ggplot(
  ggplot2::aes(x = Time/result$ReactionTime, y = exp(AvgLogAbundance),
               color = as.character(Species), group = Species)
) + ggplot2::geom_line(
) + ggplot2::labs(
  x = "Time (Units: Characteristic Times)",
  y = "Average Abundance (Geometric)"
) + ggplot2::facet_wrap(
  . ~ Count
) + ggplot2::theme_bw() + ggplot2::scale_y_log10(
) + ggplot2::guides(color = "none")


# Still need to do Observed Species Abundances, consumer richness,
# basal richness, and to add correlations to the plots.
#
# plot_0_Consumers <-

plot_0_JaccardTrue <- bootstrapSamples %>% dplyr::filter(
  Type == "Space for time", TimeSinceStart == 0
) %>% dplyr::group_by(
  Bootstrap, Time
) %>% dplyr::group_modify(
  .f = function(.x, .y) {
    stopifnot(length(.x$Patch) == length(unique(.x$Patch)))

    uniqueSpecies <- sort(unique(unlist(strsplit(
      .x$SamplingNonZeroSpecies,
      split = ", "))))

    comdatmat <- matrix(0,
                        nrow = length(.x$Patch),
                        ncol = length(uniqueSpecies))
    colnames(comdatmat) <- uniqueSpecies
    rownames(comdatmat) <- .x$Patch

    for(i in seq_along(.x$Patch)) {
      # print(strsplit(.x$SamplingNonZeroSpecies[i], split = ", ")[[1]])
      # print(strsplit(.x$SamplingNonZeroAbundances[i], split = ", ")[[1]])

      comdatmat[.x$Patch[i],
                strsplit(.x$SamplingNonZeroSpecies[i], split = ", ")[[1]]] <-
        as.numeric(strsplit(.x$SamplingNonZeroAbundances[i], split = ", ")[[1]])
    }

    Jacs <- vegan::vegdist(method = "jaccard", x = comdatmat > 0)

    data.frame(
      Jaccard = Jacs[1:length(Jacs)],
      Patch1 = rep(attr(Jacs, "Labels"), (length(attr(Jacs, "Labels"))-1):0),
      Patch2 = attr(Jacs, "Labels")[
        sequence(from = seq_along(attr(Jacs, "Labels"))[-1],
                 nvec = (length(attr(Jacs, "Labels")) - 1):1)
        ]
    )
  }
) %>% dplyr::group_by(
  Bootstrap, Time
) %>% dplyr::summarise(
  ymin = quantile(Jaccard, probs = 0.95),
  ymax = quantile(Jaccard, probs = 0.05),
  Jaccard = quantile(Jaccard, probs = 0.50)
) %>% ggplot2::ggplot(
  ggplot2::aes(x = Time/result$ReactionTime,
               y = Jaccard, ymin = ymin, ymax = ymax)
) + ggplot2::geom_ribbon(
  alpha = 0.25
) + ggplot2::geom_line(
) + ggplot2::labs(
  x = "Time (Units: Characteristic Times)",
  y = "True Jaccard Dissimilarity (Presence/Absence)",
  caption = "Inner 90% Interval and Median, Evaluated Only at (False) Disturbance Time"
) + ggplot2::theme_bw() + ggplot2::geom_rug(sides = "b")
