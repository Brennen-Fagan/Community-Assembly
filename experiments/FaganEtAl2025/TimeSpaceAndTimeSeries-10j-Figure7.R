# Setup: ######################################################################
# Plot of the effects of intervention overall as a set of heatmaps.
# Quite a bit of reproduction of effort from Figure 4, but for a different
# set of data.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")
source(file.path("R", "plotTextHeatmap.R"))

library(colorspace) # for diverging color scales with midpoint control.
library(scales) # conversion to percentages

figure7 <- list(
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
figure7$interventionTimes <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == figure7$pref,
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

figure7$dataBase <- tidytable::bind_rows(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure7$pref,
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0"
    # Need all subsets
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure7$pref,
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance"
    # Need all subsets
  )
) |> tidytable::left_join(
  figure7$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ),
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::mutate(
  Metric = factor(Metric, levels = c("Alpha Hill:0", "Alpha Abundance"),
                  labels = c("Richness", "Abundance"), ordered = TRUE),
  Time = round(Time - InterventionTime, 6) # remove false differences
) |> tidytable::filter(
  Time > -1000, Time <= 16000,
  # Avoid singletons.
  abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0
)

# Why to the level of summary? Because the PlotMeanAndInner function
# isn't built to handle the multiple resolutions that we have in the
# actual data, which makes it harder to portray the data accurately.
figure7$dataLogScale <- figure7$dataBase |> tidytable::filter(
  Metric %in% c("Richness", "Abundance")#,
  # is.na(Subset) # Not overall values
) |> tidytable::separate_wider_delim(
  Subset, names = c("Subset", "Affinity"), delim = "_"
) |> tidytable::group_by(
  # Need to aggregate to trophic levels before we reconcile times.
  Intervention, InterventionInitial, InterventionFinal, Metric, Subset,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = sum(Value), .groups = "drop"
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
  Time %in% c(0, figure7$heatmapTimes)
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric, Subset,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Within simulation proportions
  Intervention, InterventionInitial, InterventionFinal, Metric, Subset,
  PoolPatchSeed, SpeciesPreferences
) |> tidytable::mutate(
  Value = log10(Value/Value[Time == 0]) # Both x3, x/3 same magnitude.
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric, Subset,
  SpeciesPreferences, Time
) |> tidytable::summarise(
  Average = 10^mean(Value),
  StdDev = sd(Value),
  CI025 = 10^quantile(Value, p = 0.025, na.rm = TRUE), # Some subsets
  CI975 = 10^quantile(Value, p = 0.975, na.rm = TRUE), # have x/0.
  .groups = "drop"
) |> dplyr::mutate( # Change labelling, dplyr for conversion (can't in dt)
  Time = factor(
    Time, levels = c(0, range(figure7$heatmapTimes)),
    labels = c("Time 0", paste0(
      c("Short (t = ", "Long (t = "),
      range(figure7$heatmapTimes), ")"
    ))
  )
)

figure7$dataRawValues <- figure7$dataBase |> tidytable::filter(
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
  Time %in% c(0, figure7$heatmapTimes)
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
) |> dplyr::mutate( # Change labelling, dplyr for conversion (can't in dt)
  Time = factor(
    Time, levels = c(0, range(figure7$heatmapTimes)),
    labels = c("Time 0", paste0(
      c("Short (t = ", "Long (t = "),
      range(figure7$heatmapTimes), ")"
    ))
  )
)

##### KEY: ####################################################################

figure7$plot <- plotTextHeatmap(
  figure7$dataLogScale |> tidytable::filter(
    Time != "Time 0", is.na(Subset)
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure7$emphasise,
    Average = Average
  ),
  ""
) + colorspace::scale_fill_discrete_sequential(
  palette = "Hawaii", drop = FALSE
)

figure7$plotB <- plotTextHeatmap(
  figure7$dataLogScale |> tidytable::filter(
    Time != "Time 0", Subset == "Basal"
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure7$emphasise,
    Average = Average
  ),
  "Basal"
) + colorspace::scale_fill_discrete_sequential(
  palette = "Hawaii", drop = FALSE
) # See end comment if you think it looks weird!

figure7$plotC <- plotTextHeatmap(
  figure7$dataLogScale |> tidytable::filter(
    Time != "Time 0", Subset == "Consumer"
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure7$emphasise,
    Average = Average
  ),
  "Consumer"
) + colorspace::scale_fill_discrete_sequential(
  palette = "Hawaii", drop = FALSE
)

figure7$plotRaw <- ggpubr::ggarrange(
  plotTextHeatmap(
    figure7$dataRawValues |> tidytable::filter(
      Time != "Time 0", Metric == "Richness"
    ) |> tidytable::mutate(
      Emphasis = Intervention %in% figure7$emphasise
    ),
    legendName = "Raw Values",
    scalestrans = scales::label_number(accuracy = 0.1)
  ) + ggplot2::scale_fill_viridis_d(
    drop = FALSE
    # transform = "identity", begin = 0.1,
    # labels = scales::label_number(accuracy = 0.1)
  ),
  plotTextHeatmap(
    figure7$dataRawValues |> tidytable::filter(
      Time != "Time 0", Metric == "Abundance"
    ) |> tidytable::mutate(
      Emphasis = Intervention %in% figure7$emphasise
    ),
    legendName = "Raw Values", legendtrans = "log10",
    scalestrans = scales::label_number(accuracy = 1)
  ) + ggplot2::scale_fill_viridis_d(
    drop = FALSE
    # transform = "log10", begin = 0.1,
    # labels = scales::label_number(accuracy = 0.1)
  ),
  nrow = 2
)

ggplot2::ggsave(
  plot = figure7$plot,
  filename = file.path(dirImages, "Figure7_Prototype2.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plot,
  filename = file.path(dirImages, "Figure7_Prototype2.png"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plotB,
  filename = file.path(dirImages, "Figure7B_Prototype2.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plotC,
  filename = file.path(dirImages, "Figure7C_Prototype2.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plotRaw,
  filename = file.path(dirImages, "Figure7Raw_Prototype2.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)


### Does it look weird? ######################################################
# If you're here, you probably saw that figure 8 abundance looks oddly high. And
# this goes for the basal and consumer versions to varying extents. If you
# compare the raw values between 10 and 10000, the 10 looks a lot lower as well.
# So, are the values generating this data particularly weird? (I.e., when can we
# reject the null hypothesis that the data before and after the intervention
# point are drawn from the same distribution, esp. when there is a null
# intervention?)

# Drag the data back from before summarising:
figure7$suppTesting <- figure7$dataBase |> tidytable::filter(
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
  Time %in% c(0, figure7$heatmapTimes)
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> dplyr::group_by(
  # compare across seeds, but within each metric and intervention
  Metric, Intervention
) |> dplyr::group_map(
  function(data, key) {
    cbind(
      key,
      data.frame(
        vs = c(10, 10000),
        p.value = c( # Test against null two samples drawn from same dist.
          ks.test(data |> filter(Time == 0) |> pull(Value),
                  data |> filter(Time == 10) |> pull(Value),
                  simulate.p.value = TRUE)$p.value,
          ks.test(data |> filter(Time == 0) |> pull(Value),
                  data |> filter(Time == 10000) |> pull(Value),
                  simulate.p.value = TRUE)$p.value)
      )
    )
  }
) |> tidytable::bind_rows(
) |> tidytable::mutate(
  # Adjust for the mass testing situation.
  p.adj = p.adjust(p.value, method = "fdr")
)

# An alternative way of looking at it: how many are increasing/decreasing?
figure7$suppSummary <- figure7$dataBase |> tidytable::filter(
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
  Time %in% c(0, figure7$heatmapTimes)
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::pivot_wider(
  names_from = Time, values_from = Value
) |> tidytable::mutate(
  Ratio10 = `10`/`0`,
  Ratio10000 = `10000`/`0`
) |> tidytable::group_by(
  Metric, Intervention
) |> tidytable::summarise(
  Increase10 = mean(Ratio10 > 1),
  Neutral10 = mean(Ratio10 == 1),
  Decrease10 = mean(Ratio10 < 1),
  Increase10000 = mean(Ratio10000 > 1),
  Neutral10000 = mean(Ratio10000 == 1),
  Decrease10000 = mean(Ratio10000 < 1)
)
