# Setup: ######################################################################
# Plot of the effects of intervention overall as a set of heatmaps.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

figure5 <- list(
  pref = "100% 0", #"Uniform(0, 1)"
  heatmapTimes = c(10, 10000)
)

# Main Plots: #################################################################
### Plot 5: ###################################################################
##### Data: ###################################################################

# Interventions store the time right before intervention, then the
# time of the intervention itself. Retrieve this second time.
# Note it is per-simulation (timescale set by PoolPatch effectively).
# (If done correctly, the set-up as of 23/01/2026 means there will be
# round numbers from 1:10, evens from 12:20, and then by 3s 20:50.)
figure5$interventionTimes <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == figure5$pref,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset), # Won't matter, so less data
  InterventionInitial != InterventionFinal
) |> tidytable::select(
  PoolPatch, PoolPatchSeed, Time
) |> tidytable::group_by(
  PoolPatch, PoolPatchSeed
) |> tidytable::summarise(
  Time = min(
    Time[round(Time, 6)!=round(min(Time), 6)]
  ), # Not min(Time) but the next time.
  .groups = "drop"
)

figure5$dataBase <- tidytable::bind_rows(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure5$pref,
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0" 
    # Need all subsets
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure5$pref,
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance" 
    # Need all subsets
  )
) |> tidytable::left_join(
  figure5$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ), 
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::mutate(
  Metric = factor(Metric, levels = c("Alpha Hill:0", "Alpha Abundance"),
                  labels = c("Richness", "Abundance"), ordered = TRUE),
  Time = Time - InterventionTime
) |> tidytable::filter(
  Time > -1000, Time < 15000,
  # Avoid singletons.
  abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0 
)

# Why to the level of summary? Because the PlotMeanAndInner function
# isn't built to handle the multiple resolutions that we have in the
# actual data, which makes it harder to portray the data accurately.
figure5$dataOverallSummary <- figure5$dataBase |> tidytable::filter(
  Metric %in% c("Richness", "Abundance"),
  is.na(Subset) # Not overall values
) |> tidytable::mutate(
  Time = tidytable::case_when( # Create groupings for times.
    Time < -50 ~ round(Time, -2),
    Time < 0 ~ -25, # In the last bin before regime change.
    Time <= 50 ~ round(Time, 0),
    Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
    Time < 16350 ~ round(Time, -2),
    TRUE ~ Time
  )
) |> tidytable::filter(
  Time %in% figure5$heatmapTimes
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, Time
) |> tidytable::summarise(
  Average = mean(Value)
) 

# Ratios need to be handled slightly differently due to consumer/basal
# resulting in row changes.
figure5$dataBCSummary <- figure5$dataBase |> tidytable::filter(
  Metric %in% c("Richness", "Abundance"),
  !is.na(Subset) # Not overall values
) |> tidytable::separate_wider_delim(
  delim = "_", cols = Subset, names = c("Guild", "AffinityBins")
) |> unifyAffinityBins( # if many preference types.
) |> tidytable::group_by(
  # Aggregate Over the AffinityBins.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Guild, Time
) |> tidytable::summarise(
  Value = sum(Value)
) |> tidytable::pivot_wider(
  names_from = Guild, values_from = Value
) |> tidytable::mutate(
  Time = tidytable::case_when( # Create groupings for times.
    Time < -50 ~ round(Time, -2),
    Time < 0 ~ -25, # In the last bin before regime change.
    Time <= 50 ~ round(Time, 0),
    Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
    Time < 16350 ~ round(Time, -2),
    TRUE ~ Time
  ),
  Subset = NA,
  Value = Consumer/Basal
) |> tidytable::filter(
  Time %in% figure5$heatmapTimes
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, Time
) |> tidytable::summarise(
  Average = mean(Value)
) 


##### b: ######################################################################
# HEATMAPS: B Richness, B Abundance, 
# B Richness guild Ratio, B Abundance guild Ratio,
# B Richness time Difference, B Abundance time Difference
figure4$plotBR <- ggplot2::ggplot(
  figure4$dataOverallSummary |> tidytable::filter(
    Time %in% figure4$heatmapTimes
  )
)
figure4$plotBA
figure4$plotBRR
figure4$plotBAR
figure4$plotBRD
figure4$plotBAD