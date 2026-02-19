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

library(scales)

# This is better as an environment, but that's more opaque.
figure7 <- list(
  CI = 0.75,
  pref = "Uniform(0, 1)",
  luinitl = "(0.5)", # Land Use INITiaL
  lufinal = c("(0)", "(0.5)", "(1)") # Land Use FINAL
  # lufinal = c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)") # Land Use FINAL
)

figure7$prefstring <- switch(
  figure7$pref,
  "100% 0" = "", # Base Case
  "50% 0, 50% 1" = "_5050",
  "Uniform(0, 1)" = "_Unif"
)
figure7$lustring <- paste0(switch(
  figure7$luinitl,
  "(0)" = "0",
  "(0.5)" = "", # Base Case
  "(1)" = "1"
), "to", length(figure7$lufinal))

# Main Plots: #################################################################
### Plot 7: ###################################################################
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
    InterventionInitial == figure7$luinitl,
    InterventionFinal %in% figure7$lufinal,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0"
    # Need all subsets
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure7$pref,
    NicheDistance == defaultNicheDistance,
    InterventionInitial == figure7$luinitl,
    InterventionFinal %in% figure7$lufinal,
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
  Time <= 16000, # Need the start for the inset.
  # Avoid singletons.
  abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0
)

# Why to the level of summary? Because the PlotMeanAndInner function
# isn't built to handle the multiple resolutions that we have in the
# actual data, which makes it harder to portray the data accurately.
figure7$dataOverallSummary <- figure7$dataBase |> tidytable::filter(
  Time > -1000,
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
  Lower = quantile(Value, p = (1 - figure7$CI) - (1 - figure7$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure7$CI + (1 - figure7$CI)/2, na.rm = TRUE)
)

# Same idea as the overall case, but split by guild.
figure7$dataBCSupplement <- figure7$dataBase |> tidytable::filter(
  Time > -1000,
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
  Lower = quantile(Value, p = (1 - figure7$CI) - (1 - figure7$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure7$CI + (1 - figure7$CI)/2, na.rm = TRUE)
)

##### A: ################################################################
figure7$plotA <- plotMeanAndInner(
  figure7$dataBase |> tidytable::filter(
    Metric == "Richness",
    is.na(Subset) # Not overall values
  ) |> tidytable::mutate(
    Time = Time + InterventionTime
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_vline(
  xintercept = min(figure7$interventionTimes$Time),
  color = "black", linetype = "dashed"
) + ggplot2::labs(
  x = "Time",
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = "none", # already covered by the main plot.
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax),
  xlim = c(0, 31000),
  expand = FALSE
)

##### B: ######################################################################
figure7$plotB <- ggplot2::ggplot(
  figure7$dataOverallSummary |> tidytable::filter(
    Metric == "Richness", Time >= 0
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
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
) + ggplot2::scale_x_continuous(
  trans = "log1p",
  breaks = 10^(0:4)
)

##### C+D: ####################################################################
figure7$plotCD <- ggplot2::ggplot(
  figure7$dataBCSupplement |> tidytable::filter(
    Metric == "Richness", Time >= 0
  ) |> dplyr::mutate(
    Guild = factor(Guild, ordered = TRUE, levels = c("Consumer", "Basal"))
    # Top-Bottom ordering
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
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
) + ggplot2::scale_x_continuous(
  trans = "log1p",
  breaks = 10^(0:4)
) + ggplot2::facet_grid(
  Guild~.
)

##### E: ######################################################################
figure7$plotE <- ggplot2::ggplot(
  figure7$dataOverallSummary |> tidytable::filter(
    Metric == "Abundance", Time >= 0
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
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
  ylim = c(1e-2, 1e5),
  expand = FALSE
) + ggplot2::scale_y_log10(
  label = scales::label_log()
) + ggplot2::scale_x_continuous(
  trans = "log1p",
  breaks = 10^(0:4)
)

##### F+G: ####################################################################
figure7$plotFG <- ggplot2::ggplot(
  figure7$dataBCSupplement |> tidytable::filter(
    Metric == "Abundance", Time >= 0
  ) |> dplyr::mutate(
    Guild = factor(Guild, ordered = TRUE, levels = c("Consumer", "Basal"))
    # Top-Bottom ordering
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
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
  ylim = c(1e-2, 1e5),
  expand = FALSE
) + ggplot2::scale_x_continuous(
  trans = "log1p",
  breaks = 10^(0:4)
) + ggplot2::scale_y_log10(
  label = scales::label_log()
) + ggplot2::facet_grid(
  Guild~.
)

##### Combine: ################################################################
figure7$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure7$plotA + ggplot2::scale_y_continuous(
      breaks = c(0, 10, 20, 30, 40),
      labels = c("0", "10", "20", "30", "")
    ),
    ggpubr::ggarrange(plotlist = list(
      figure7$plotB + ggplot2::theme(legend.position = "none"),
      figure7$plotCD + ggplot2::theme(legend.position = "none",
                                      axis.title.y = ggplot2::element_blank())
    ), nrow = 1),
    ggpubr::ggarrange(plotlist = list(
      figure7$plotE + ggplot2::theme(legend.position = "none",
                                     axis.text.y = ggplot2::element_text(
                                       margin = ggplot2::margin(0, 0, 0, -5)
                                     )),
      figure7$plotFG + ggplot2::theme(legend.position = "none",
                                      axis.title.y = ggplot2::element_blank(),
                                      axis.text.y = ggplot2::element_text(
                                        margin = ggplot2::margin(0, 0, 0, -5)
                                      ))
    ), nrow = 1)
  ), nrow = 3, common.legend = TRUE
)

figure7$suffix <- paste0("_", figure7$prefstring, "_", figure7$lustring)
figure7$prefix <- "figure7"
figure7$iter <- "_Prototype2"
figure7$ids <- c("", "", "A", "B", "CD", "E", "FG")
figure7$ext <- c(".png", rep(".pdf", 6))

for (fnum in 1:7) {
  with(figure7, ggplot2::ggsave(
    plot = switch(fnum,
                  plot, plot,
                  plotA, plotB, plotCD, plotE, plotFG),
    filename = file.path(dirImages,
                         paste0(prefix, ids[fnum], iter, suffix, ext[fnum])),
    units = "cm", width = 6.5*3, height = 6.5*2.1)
  )
}
