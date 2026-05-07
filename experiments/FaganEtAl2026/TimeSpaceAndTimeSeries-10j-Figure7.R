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

figure7$timeDiameter <- 10^(log10(max(figure7$heatmapTimes))-1)

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
  Time >= -figure7$timeDiameter, Time <= 16000,
  # Avoid singletons.
  abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time <= 0
)

# Why to the level of summary? Because the PlotMeanAndInner function
# isn't built to handle the multiple resolutions that we have in the
# actual data, which makes it harder to portray the data accurately.

# We want to do two analogous but not identical comparisons.
# One is the shorter drop-off between two points that have autocorrelation.
# The other is long-run behaviours.
# The differing goals, and differing amount of data we can use, suggests we
# should look at two separate but similar analyses.
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
  # Time = tidytable::case_when( # Create groupings for times, avoid doublecount
  #   Time < -50 ~ round(Time, -2),
  #   Time < 0 ~ -25, # In the last bin before regime change.
  #   Time <= 50 ~ round(Time, 0),
  #   Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
  #   Time < 16350 ~ round(Time, -2),
  #   TRUE ~ Time
  # )
  Time = tidytable::case_when(
    round(Time) == 0 ~ 0, # Base Case for reference point analysis
    Time < 0 ~ -1, # Base Case for reference time segment analysis
    round(Time) == min(figure7$heatmapTimes) ~
      min(figure7$heatmapTimes), # Contrast point
    max(figure7$heatmapTimes) - figure7$timeDiameter/2 < Time &
      max(figure7$heatmapTimes) + figure7$timeDiameter/2 > Time ~
      max(figure7$heatmapTimes), # Contrast time segment
    TRUE ~ Time # To Discard.
  ),
  Analysis = ifelse(Time %in% c(0, min(figure7$heatmapTimes)),
                    "Point", "Segment")
) |> tidytable::filter( # Discard To-Discards.
  Time %in% c(-1, 0, figure7$heatmapTimes)
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric, Subset,
  PoolPatchSeed, SpeciesPreferences, Time, Analysis
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
)

# Because only non-intervention systems have a past (-> -1), we need to provide
# the intervention systems with their corresponding past before continuing.
figure7$dataLogScale <- tidytable::bind_rows(
  figure7$dataLogScale,
  figure7$dataLogScale |> tidytable::filter(
    Time == -1
  ) |> tidytable::select(
    -Intervention, -InterventionFinal
  ) |> tidytable::full_join(
    tidytable::expand_grid(
      InterventionInitial = unique(figure7$dataLogScale$InterventionInitial),
      InterventionFinal = unique(figure7$dataLogScale$InterventionFinal)
    )
  ) |> tidytable::filter(
    InterventionInitial != InterventionFinal
  ) |> tidytable::mutate(
    Intervention = paste(InterventionInitial, InterventionFinal, sep = "->")
  )
) |> tidytable::group_by(
  # Within simulation proportions
  Intervention, InterventionInitial, InterventionFinal, Metric, Subset,
  PoolPatchSeed, SpeciesPreferences, Analysis
) |> tidytable::mutate(
  Value = log10(Value/Value[Time == min(Time)]) # Both x3, x/3 same magnitude.
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric, Subset,
  SpeciesPreferences, Time
) |> tidytable::summarise(      # Been told we want it 0-centred.
  Average = 10^mean(Value) - 1, # Abusing that (new-old)/old = (new/old) - 1.
  StdDev = sd(Value),           # full: 10^(sum(log10(new/old))/length(new))-1.
  CI025 = 10^quantile(Value, p = 0.025, na.rm = TRUE), # Some subsets
  CI975 = 10^quantile(Value, p = 0.975, na.rm = TRUE), # have x/0.
  .groups = "drop"
) |> dplyr::mutate( # Change labelling, dplyr for conversion (can't in dt)
  Time = factor(
    Time, levels = c(-1, 0, range(figure7$heatmapTimes)),
    labels = c("Time Past", "Time 0", paste0(
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
    Time != "Time 0", Time != "Time Past", is.na(Subset)
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure7$emphasise,
    Average = Average
  ),
  ""
) + colorspace::scale_fill_discrete_qualitative(
  palette = "Set 2", drop = FALSE
)

figure7$plotB <- plotTextHeatmap(
  figure7$dataLogScale |> tidytable::filter(
    Time != "Time 0", Time != "Time Past", Subset == "Basal"
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure7$emphasise,
    Average = Average
  ),
  "Basal"
) + colorspace::scale_fill_discrete_qualitative(
  palette = "Set 2", drop = FALSE
)

figure7$plotC <- plotTextHeatmap(
  figure7$dataLogScale |> tidytable::filter(
    Time != "Time 0", Time != "Time Past", Subset == "Consumer"
  ) |> tidytable::mutate(
    Emphasis = Intervention %in% figure7$emphasise,
    Average = Average
  ),
  "Consumer"
) + colorspace::scale_fill_discrete_qualitative(
  palette = "Set 2", drop = FALSE
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
  filename = file.path(dirImages, "Figure7_Prototype3.pdf"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plot,
  filename = file.path(dirImages, "Figure7_Prototype3.png"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plotB,
  filename = file.path(dirImages, "FigureS8_Basal.png"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plotC,
  filename = file.path(dirImages, "FigureS8_Consumer.png"),
  units = "cm", width = 6.5*4, height = 6.5*3)
ggplot2::ggsave(
  plot = figure7$plotRaw,
  filename = file.path(dirImages, "Figure7Raw_Prototype3.png"),
  units = "cm", width = 6.5*4, height = 6.5*3)

