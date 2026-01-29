# Setup: ######################################################################
# Plot of the effects of intervention from (0.5) on 100% (0) species case.
# This means a) summarised richness time series from time of intervention,
# b) pseudo-heatmap of short and long time scales (statistic intensities, not
# counts), c) Consumer-Basal Richness ratio and d) Abundance ratio through
# short time scales, e) supplement with short to long-term transition
# for richness and abundance on log(1+Time) scales, and f) repeat all of the
# above for from (0) and from (1).

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

# This is better as an environment, but that's more opaque.
figure4 <- list(
  CI = 0.75,
  pref = "100% 0", #"Uniform(0, 1)"
  luinitl = "(0.5)", # Land Use INITiaL
  lufinal = c("(0)", "(0.5)", "(1)") # Land Use FINAL
  # lufinal = c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)") # Land Use FINAL
)

figure4$prefstring <- switch(
  figure4$pref,
  "100% 0" = "", # Base Case
  "50% 0, 50% 1" = "_5050",
  "Uniform(0, 1)" = "_Unif"
)
figure4$lustring <- paste0(switch(
  figure4$luinitl,
  "(0)" = "0",
  "(0.5)" = "", # Base Case
  "(1)" = "1"
), "to", length(figure4$lufinal))

# Main Plots: #################################################################
### Plot 4: ###################################################################
##### Data: ###################################################################

# Interventions store the time right before intervention, then the
# time of the intervention itself. Retrieve this second time.
# Note it is per-simulation (timescale set by PoolPatch effectively).
# (If done correctly, the set-up as of 23/01/2026 means there will be
# round numbers from 1:10, evens from 12:20, and then by 3s 20:50.)
figure4$interventionTimes <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == figure4$pref,
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

figure4$dataBase <- tidytable::bind_rows(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure4$pref,
    NicheDistance == defaultNicheDistance,
    InterventionInitial == figure4$luinitl,
    InterventionFinal %in% figure4$lufinal,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0"
    # Need all subsets
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure4$pref,
    NicheDistance == defaultNicheDistance,
    InterventionInitial == figure4$luinitl,
    InterventionFinal %in% figure4$lufinal,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance"
    # Need all subsets
  )
) |> tidytable::left_join(
  figure4$interventionTimes |> tidytable::rename(
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
figure4$dataOverallSummary <- figure4$dataBase |> tidytable::filter(
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
  Lower = quantile(Value, p = (1 - figure4$CI) + (1 - figure4$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure4$CI + (1 - figure4$CI)/2, na.rm = TRUE)
)

# Ratios need to be handled slightly differently due to consumer/basal
# resulting in row changes.
figure4$dataBCSummary <- figure4$dataBase |> tidytable::filter(
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
  Lower = quantile(Value, p = (1 - figure4$CI) + (1 - figure4$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure4$CI + (1 - figure4$CI)/2, na.rm = TRUE)
)

# Same idea as the overall case, but split by guild.
figure4$dataBCSupplement <- figure4$dataBase |> tidytable::filter(
  Metric %in% c("Richness", "Abundance"),
  !is.na(Subset) # Not overall values
) |> tidytable::separate_wider_delim(
  delim = "_", cols = Subset, names = c("Guild", "AffinityBins")
) |> unifyAffinityBins( # if many preference types.
) |> tidytable::mutate(
  Time = tidytable::case_when( # Create groupings for times.
    Time < -50 ~ round(Time, -2),
    Time < 0 ~ -25, # In the last bin before regime change.
    Time <= 50 ~ round(Time, 0),
    Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
    Time < 16350 ~ round(Time, -2),
    TRUE ~ Time
  ),
  Subset = NA
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Guild, AffinityBins, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, Guild, AffinityBins, Time
) |> tidytable::summarise(
  Lower = quantile(Value, p = (1 - figure4$CI) + (1 - figure4$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure4$CI + (1 - figure4$CI)/2, na.rm = TRUE)
)

# As in dataBCSummary, but broken up by AffinityBins
figure4$dataBCSupplement2 <- figure4$dataBase |> tidytable::filter(
  Metric %in% c("Richness", "Abundance"),
  !is.na(Subset) # Not overall values
) |> tidytable::separate_wider_delim(
  delim = "_", cols = Subset, names = c("Guild", "AffinityBins")
) |> unifyAffinityBins( # if many preference types.
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
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, AffinityBins, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, AffinityBins, Time
) |> tidytable::summarise(
  Lower = quantile(Value, p = (1 - figure4$CI) + (1 - figure4$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure4$CI + (1 - figure4$CI)/2, na.rm = TRUE)
)

##### a: ######################################################################
# Summarised time-series plot for overall richness.
figure4$plotA <- ggplot2::ggplot(
  figure4$dataOverallSummary |> tidytable::filter(
    Metric == "Richness"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax),
  expand = FALSE
)

##### b: ######################################################################
# Summarised time-series plot for short time scale richness ratio.
# Note this should be Consumer / Basal (because Consumers often -> 0).

figure4$plotB <- ggplot2::ggplot(
  figure4$dataBCSummary |> tidytable::filter(
    Metric == "Richness",
    Time > -200, Time < 100
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::coord_cartesian(
  xlim = c(0, 50), expand = FALSE
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness (Cons./Basal)"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
)

##### c: ######################################################################
# Summarised time-series plot for short time scale abundance ratio.
# Note this should be Consumer / Basal (because Consumers often -> 0).

figure4$plotC <- ggplot2::ggplot(
  figure4$dataBCSummary |> tidytable::filter(
    Metric == "Abundance",
    Time > -200, Time < 100
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::coord_cartesian(
  xlim = c(0, 50), expand = FALSE
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Abundance (Cons./Basal)"
) + ggplot2::scale_y_log10(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
)

##### SUPPLEMENT: #############################################################
# Short-long term transition for richness and abundance on log(1+Time) scale
# expecting to hit the ratios as well as the base values.

figure4$SupplementRichness <- ggplot2::ggplot(
  figure4$dataBCSupplement |> tidytable::filter(
    Metric == "Richness"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  factor(Guild, levels = c("Consumer", "Basal"), ordered = TRUE) ~
    AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
)

figure4$SupplementAbundance <- ggplot2::ggplot(
  figure4$dataBCSupplement |> tidytable::filter(
    Metric == "Abundance"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Abundance"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  factor(Guild, levels = c("Consumer", "Basal"), ordered = TRUE) ~
    AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
) + ggplot2::scale_y_log10(
)

figure4$SupplementRichnessRatio <- ggplot2::ggplot(
  figure4$dataBCSupplement2 |> tidytable::filter(
    Metric == "Richness"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness (Cons./Basal)"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  . ~
    AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
)

figure4$SupplementAbundanceRatio <- ggplot2::ggplot(
  figure4$dataBCSupplement2 |> tidytable::filter(
    Metric == "Abundance"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Abundance (Cons./Basal)"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  . ~ AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
) + ggplot2::scale_y_log10(
)

##### Combine: ################################################################
figure4$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure4$plotA,
    ggpubr::ggarrange(plotlist = list(
      figure4$plotB + ggplot2::theme(legend.position = "none"),
      figure4$plotC + ggplot2::theme(legend.position = "none")
    ), ncol = 1)
  ), nrow = 1, common.legend = TRUE
)

figure4$suffix <- paste0("_", figure4$prefstring, "_", figure4$lustring)
figure4$prefix <- "Figure4"
figure4$iter <- "_Prototype5"
figure4$ids <- c("", "", "A", "B", "C", "SR", "SA", "SRR", "SAR")
figure4$ext <- c(".png", rep(".pdf", 8))

for (fnum in 1:(2+3+4)) {
  with(figure4, ggplot2::ggsave(
    plot = switch(fnum,
                  plot, plot,
                  plotA, plotB, plotC,
                  SupplementRichness, SupplementAbundance,
                  SupplementRichnessRatio, SupplementAbundanceRatio),
    filename = file.path(dirImages,
                         paste0(prefix, ids[fnum], iter, suffix, ext[fnum])),
    units = "cm", width = 6.5*3, height = 6.5*2)
  )
}
