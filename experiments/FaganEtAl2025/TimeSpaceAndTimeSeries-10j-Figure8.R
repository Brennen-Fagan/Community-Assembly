# Setup: ######################################################################
# Plot of the effects of intervention overall as a set of heatmaps.
# Quite a bit of reproduction of effort from Figure 4, but for a different
# set of data.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

library(colorspace) # for diverging color scales with midpoint control.

figure8 <- list(
  pref = "Uniform(0, 1)",
  heatmapTimes = c(10, 10000),
  emphasise = c("(0.5)", "(0.5)->(0)", "(0.5)->(1)")
)

# Main Plots: #################################################################
### Plot 8: ###################################################################
##### Data: ###################################################################

# Interventions store the time right before intervention, then the
# time of the intervention itself. Retrieve this second time.
# Note it is per-simulation (timescale set by PoolPatch effectively).
# (If done correctly, the set-up as of 23/01/2026 means there will be
# round numbers from 1:10, evens from 12:20, and then by 3s 20:50.)
figure8$interventionTimes <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == figure8$pref,
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

figure8$dataBase <- tidytable::bind_rows(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure8$pref,
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0"
    # Need all subsets
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure8$pref,
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance"
    # Need all subsets
  )
) |> tidytable::left_join(
  figure8$interventionTimes |> tidytable::rename(
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
figure8$dataOverallSummary <- figure8$dataBase |> tidytable::filter(
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
  Time %in% c(0, figure8$heatmapTimes)
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
figure8$dataBCSummary <- figure8$dataBase |> tidytable::filter(
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
  Time %in% c(0, figure8$heatmapTimes)
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

##### FUNCTION: ###############################################################
# Different Scales mean we have to separate out the data, so we define a
# function to perform the plotting repeatedly/consistently.
plotTextHeatmap <- function(data, legendName, legendtrans = "identity") {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = InterventionInitial,
      y = InterventionFinal,
      fill = Average
    )
  ) + ggplot2::geom_tile(
    width = 1, height = 1, color = NA
  ) + ggplot2::geom_tile(
    data = function(x) x |> tidytable::filter(Emphasis),
    fill = NA, color = "black", linewidth = 1
  ) + ggplot2::geom_text(
    ggplot2::aes(label = Average)
  ) + ggplot2::facet_grid(
    Metric ~ Time
  ) + ggplot2::scale_fill_viridis_c(
    transform = legendtrans, begin = 0.1
  ) + ggplot2::theme_minimal(
  ) + ggplot2::labs(
    fill = legendName,
    x = "Initial Habitat Type",
    y = "Final Habitat Type"
  )
}

##### KEY: ####################################################################
# HEATMAPS: Richness, Abundance,
# Richness guild Ratio, Abundance guild Ratio,
# Richness time Difference, Abundance time Difference

figure8$plotR <- plotTextHeatmap(
  figure8$dataOverallSummary |> tidytable::filter(
    Metric == "Richness", Time != 0
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    Average = signif(Average, digits = 2)
  ),
  "Richness"
)
figure8$plotRR <- plotTextHeatmap(
  figure8$dataBCSummary |> tidytable::filter(
    Metric == "Richness", Time != 0
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    Average = signif(Average, digits = 2)
  ),
  "Richness\nRatio"
)
figure8$plotA <- plotTextHeatmap(
  figure8$dataOverallSummary |> tidytable::filter(
    Metric == "Abundance", Time != 0
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    Average = signif(Average, digits = 2)
  ),
  "Abundance", "log10"
)
figure8$plotAR <- plotTextHeatmap(
  figure8$dataBCSummary |> tidytable::filter(
    Metric == "Abundance", Time != 0
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    Average = signif(Average, digits = 2)
  ),
  "Abundance\nRatio", "log10"
)


figure8$plotRD <- plotTextHeatmap(
  figure8$dataOverallSummary |> tidytable::filter(
    Metric == "Richness"
  ) |> dplyr::mutate(
    Time = paste("Time", Time)
  ) |> tidytable::pivot_wider(
    names_from = Time, values_from = Average
  ) |> tidytable::rename(
    "Initial" = `Time 0`
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    tidytable::across(
      tidytable::starts_with("Time"),
      function(x, init) signif(x - init, digits = 2),
      init = Initial
    )
  ) |> tidytable::pivot_longer(
    names_to = "Time", values_to = "Average", # Works because division by same
    # number in both fractions, so we can rearrange numerators for equiv.
    cols = tidytable::starts_with("Time")
  ) |> tidytable::mutate(
    # Remove the "Time" tag for consistency
    Time = as.numeric(substring(Time, first = 5))
  ),
  "Richness\nDifference"
) + colorspace::scale_fill_continuous_diverging(
  palette = "Green-Brown", mid = 0
)

figure8$plotAD <- plotTextHeatmap(
  figure8$dataOverallSummary |> tidytable::filter(
    Metric == "Abundance"
  ) |> dplyr::mutate(
    Time = paste("Time", Time)
  ) |> tidytable::pivot_wider(
    names_from = Time, values_from = Average
  ) |> tidytable::rename(
    "Initial" = `Time 0`
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    tidytable::across(
      tidytable::starts_with("Time"),
      function(x, init) signif(x - init, digits = 2),
      init = Initial
    )
  ) |> tidytable::pivot_longer(
    names_to = "Time", values_to = "Average", # Works because division by same
    # number in both fractions, so we can rearrange numerators for equiv.
    cols = tidytable::starts_with("Time")
  ) |> tidytable::mutate(
    # Remove the "Time" tag for consistency
    Time = as.numeric(substring(Time, first = 5))
  ),
  "Abundance\nDifference"
) + colorspace::scale_fill_continuous_diverging(
  palette = "Green-Brown", mid = 0
)


##### Combine: ################################################################
figure8$plotRichness <- ggpubr::ggarrange(
  plotlist = list(
    figure8$plotR,
    figure8$plotRR
  ), nrow = 2
)
figure8$plotSupplement <- ggpubr::ggarrange(
  plotlist = list(
    figure8$plotA,
    figure8$plotAR
  ), nrow = 2
)

ggplot2::ggsave(
  plot = figure8$plotRichness,
  filename = file.path(dirImages, "figure8_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotRichness,
  filename = file.path(dirImages, "figure8_Prototype1.png"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotSupplement,
  filename = file.path(dirImages, "figure8S_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotR,
  filename = file.path(dirImages, "figure8R_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotRR,
  filename = file.path(dirImages, "figure8RR_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotA,
  filename = file.path(dirImages, "figure8A_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotAR,
  filename = file.path(dirImages, "figure8AR_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotAD,
  filename = file.path(dirImages, "figure8AD_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(
  plot = figure8$plotRD,
  filename = file.path(dirImages, "figure8RD_Prototype1.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2)
