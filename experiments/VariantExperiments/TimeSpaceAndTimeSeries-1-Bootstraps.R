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
  # 12854863 # Used for 2023-09-25
  # 99701559 # Used for Hik6_2023-03-04, 56-6
  11460190 # Used for 2023-11-13, Medium, not HIK6

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

library(vegan) # vegdist

library(betapart)
library(car)
library(lme4)

clust <- parallel::makeCluster(3)#, outfile = "")
doParallel::registerDoParallel(clust)

# Files: ######################################################################
directory <- '.' # Mutualism folder.
file_result <-
  #"MNA-Hik6A-Cases-Prepared-1-56-6-1.RData"
  "MNA-ExampleExtProp-Result-Env10-Ring-5-1-1-ExtProp1.RData"
file_result2 <- NULL # "MNA-Hik6A-Cases-Prepared-1-56-6-2.RData"
file_pool <-
  #"MNA-Hik6A-Cases-Prepared.RData"
  "MNA-ExampleOutcome-PoolMats-Env10.RData"
load(file.path(directory, "Data_2023-11-13", file_result))
load(file.path(directory, "Data_2023-11-13", file_pool))

if (!is.null(file_result2)) {
  results1 <- results
  load(file.path(directory, "Data_Hik6_2023-03-04", file_result2))

  EntriesToCheck <- !names(results) %in% c("Events", "Abundance")

  stopifnot(isTRUE(all.equal(results[EntriesToCheck],
                             results1[EntriesToCheck])))

  results1$Abundance <- results1$Abundance[
    results1$Abundance[, 1] < min(results$Abundance[, 1]),
  ]

  results$Events <- rbind(results1$Events, results$Events)
  results$Events <- results$Events %>% dplyr::distinct()
  results$Abundance <- rbind(results1$Abundance, results$Abundance)

  result <- results # Note the "s", there was a change in formatting!

  rm(results1)
  rm(results)
}

if (exists("pools") && !exists("Pool")) {
  IDNumbers <- sub(file_result, pattern = ".RData", replacement = "")
  IDNumbers <- strsplit(IDNumbers, split = "-", fixed = TRUE)[[1]]
  IDNumbers <- IDNumbers[(which(IDNumbers == "Prepared") + 1):length(IDNumbers)]
  IDNumbers <- as.numeric(IDNumbers)
  Pool <- pools[[cases$Parameters[IDNumbers[2]]]][[cases$System[IDNumbers[2]]]]
  # Why 2? 1 == File / Main Case, 2 == row of cases (derived from row of CSV)
  # 3 == The seeds used, 4 == which part of the simulation was saved.
  # Note 3 and 4 are optional if all simulations of a row were saved together.
}

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

source("TimeSpaceAndTimeSeries-0-Functions.R")

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

lastTimeSampleable <- result$Events$Times[length(result$Events$Times)]
lastTimeSampleable <- lastTimeSampleable - interventionGap * 2
firstTimeSampleable <- 1000 # to make sure we're past the simulation burnin.

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
  intervention <- burnin + interventionGap # When land use change "occurs".
  # print(paste(bootstrapID, ":", toString(intervention)))

  # Time Series:
  # The premise is that the Time Series researcher is only looking at the
  # plots where the experiment / land use change is taking place, but they
  # instead study the sites before the land use change occurs.
  # They'd be interested in trying to estimate the regional pool / diversity,
  # as well as local diversity and inter-patch diversity.
  sampling_TimeSeries <- data.frame(expand.grid(
    Time = unique(c(intervention - timediffs,
                    intervention + timediffs)),
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
    Time = intervention + timediffs,
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

# Plotting: ###################################################################
### Plot 0: Sense Checking: ###################################################
# Richness Comparison
plot_0_Richness <- ggplot2::ggplot(
  bootstrapSamples %>% dplyr::mutate(
    TrueRichness = unlist(lapply(
      strsplit(SamplingNonZeroSpecies, split = ", ", fixed = TRUE), length
    ))
  ),
  ggplot2::aes(x = Time/result$ReactionTime, y = TrueRichness,
               color = Patch, group = Patch)
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
) + ggplot2::coord_cartesian(
  ylim = c(0, 30)
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
  .f = SpeciesStringsToBeta
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
) + ggplot2::theme_bw() + ggplot2::geom_rug(
  sides = "b"
  ) + ggplot2::coord_cartesian(
  ylim = c(0, 1)
)
