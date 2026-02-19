# Setup: ######################################################################
# Plot of the effects of intervention overall as a set of heatmaps.
# Quite a bit of reproduction of effort from Figure 4, but for a different
# set of data.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

library(colorspace) # for diverging color scales with midpoint control.
library(scales) # conversion to percentages

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
  Time = round(Time - InterventionTime, 6) # remove false differences
) |> tidytable::filter(
  Time > -1000, Time <= 16000,
  # Avoid singletons.
  abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0
)

# Why to the level of summary? Because the PlotMeanAndInner function
# isn't built to handle the multiple resolutions that we have in the
# actual data, which makes it harder to portray the data accurately.
figure8$dataLogScale <- figure8$dataBase |> tidytable::filter(
  Metric %in% c("Richness", "Abundance")#,
  # is.na(Subset) # Not overall values
) |> tidytable::separate_wider_delim(
  Subset, names = c("Subset", "Affinity"), delim = "_"
) |> tidytable::group_by(
  # Need to aggregate over trophic levels before we reconcile times.
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
  Time %in% c(0, figure8$heatmapTimes)
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
    Time, levels = c(0, range(figure8$heatmapTimes)),
    labels = c("Time 0", paste0(
      c("Short (t = ", "Long (t = "),
      range(figure8$heatmapTimes), ")"
    ))
  )
)

figure8$dataRawValues <- figure8$dataBase |> tidytable::filter(
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
) |> dplyr::mutate( # Change labelling, dplyr for conversion (can't in dt)
  Time = factor(
    Time, levels = c(0, range(figure8$heatmapTimes)),
    labels = c("Time 0", paste0(
      c("Short (t = ", "Long (t = "),
      range(figure8$heatmapTimes), ")"
    ))
  )
)

##### FUNCTION: ###############################################################
# Different Scales mean we have to separate out the data, so we define a
# function to perform the plotting repeatedly/consistently.
plotTextHeatmap <- function(
  data, legendName,
  legendtrans = "identity",
  scalestrans = scales::label_percent(accuracy = 0.1)
) {
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
    ggplot2::aes(label = scalestrans(Average))
  ) + ggplot2::facet_grid(
    Metric ~ Time
  ) + ggplot2::scale_fill_viridis_c(
    transform = legendtrans, begin = 0.1, labels = scalestrans
  ) + ggplot2::theme_minimal(
  ) + ggplot2::labs(
    fill = legendName,
    x = "Initial Habitat Type",
    y = "Final Habitat Type"
  ) + ggplot2::coord_fixed(
  )
}

##### KEY: ####################################################################

figure8$plot <- plotTextHeatmap(
  figure8$dataLogScale |> tidytable::filter(
    Time != "Time 0", is.na(Subset)
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    Average = Average
  ),
  ""
) + colorspace::scale_fill_continuous_diverging(
  palette = "Blue-Yellow 3", mid = log10(1), transform = "log10",
  labels = scales::label_percent(accuracy = 0.1)
  # Old Version of the package, github.com/tidyverse/ggplot2/issues/3198
)

figure8$plotB <- plotTextHeatmap(
  figure8$dataLogScale |> tidytable::filter(
    Time != "Time 0", Subset == "Basal"
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    Average = Average
  ),
  "Basal"
) + colorspace::scale_fill_continuous_diverging(
  palette = "Blue-Yellow 3", mid = log10(1), transform = "log10",
  labels = scales::label_percent(accuracy = 0.1)
  # Old Version of the package, github.com/tidyverse/ggplot2/issues/3198
) # See end comment if you think it looks weird!

figure8$plotC <- plotTextHeatmap(
  figure8$dataLogScale |> tidytable::filter(
    Time != "Time 0", Subset == "Consumer"
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure8$emphasise,
    Average = Average
  ),
  "Consumer"
) + colorspace::scale_fill_continuous_diverging(
  palette = "Blue-Yellow 3", mid = log10(1), transform = "log10",
  labels = scales::label_percent(accuracy = 0.1)
  # Old Version of the package, github.com/tidyverse/ggplot2/issues/3198
)

figure8$plotRaw <- ggpubr::ggarrange(
  plotTextHeatmap(
    figure8$dataRawValues |> tidytable::filter(
      Time != "Time 0", Metric == "Richness"
    ) |> tidytable::mutate(
      Emphasis = Intervention %in% figure8$emphasise
    ),
    legendName = "Raw Values",
    scalestrans = scales::label_number(accuracy = 0.1)
  ),
  plotTextHeatmap(
    figure8$dataRawValues |> tidytable::filter(
      Time != "Time 0", Metric == "Abundance"
    ) |> tidytable::mutate(
      Emphasis = Intervention %in% figure8$emphasise
    ),
    legendName = "Raw Values", legendtrans = "log10",
    scalestrans = scales::label_number(accuracy = 1)
  ),
  nrow = 2
)

ggplot2::ggsave(
  plot = figure8$plot,
  filename = file.path(dirImages, "figure8_Prototype2.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure8$plot,
  filename = file.path(dirImages, "figure8_Prototype2.png"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure8$plotB,
  filename = file.path(dirImages, "figure8B_Prototype2.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure8$plotC,
  filename = file.path(dirImages, "figure8C_Prototype2.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure8$plotRaw,
  filename = file.path(dirImages, "figure8Raw_Prototype2.pdf"),
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
figure8$suppTesting <- figure8$dataBase |> tidytable::filter(
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
