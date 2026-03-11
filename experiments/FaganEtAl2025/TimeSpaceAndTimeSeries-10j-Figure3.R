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
library(patchwork)

# This is better as an environment, but that's more opaque.
figure3 <- list(
  CI = 0.75,
  graph = list(
    seed = "2", # "11", "17", "2"!,
    time = 25000
  ),
  pref = "100% 0", #"Uniform(0, 1)"
  luinitl = "(0.5)", # Land Use INITiaL
  lufinal = c("(0)", "(0.5)", "(1)") # Land Use FINAL
  # lufinal = c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)") # Land Use FINAL
)

figure3$prefstring <- switch(
  figure3$pref,
  "100% 0" = "", # Base Case
  "50% 0, 50% 1" = "_5050",
  "Uniform(0, 1)" = "_Unif"
)
figure3$lustring <- paste0(switch(
  figure3$luinitl,
  "(0)" = "0",
  "(0.5)" = "", # Base Case
  "(1)" = "1"
), "to", length(figure3$lufinal))

figure3$graph$specification <- diversitiesRichness |> tidytable::select(c(
  # Which network:
  "Time", "Environment1",
  # Which File (Base):
  "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
  "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
  "Dispersal", "NicheDistance",
  # Oops, there was a collision causing human readable to replace machine.
  # Will be replaced SpeciesAffinity#2 will -> SpeciesPreferences.
  "SpeciesAffinity",
  "SpeciesAffinitySeed", "PatchAffinity", "PatchAffinitySeed",
  # Which File (Intervention):
  "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
  "InterventionTimeSeed", "InterventionDispersal", "InterventionNicheDistance",
  # Ease of Use
  "SpeciesPreferences",
  "Intervention", "InterventionInitial", "InterventionFinal"
)) |> tidytable::filter(
  SpeciesPreferences == figure3$pref,
  NicheDistance == defaultNicheDistance,
  InterventionInitial %in% figure3$luinitl,
  InterventionFinal %in% figure3$lufinal,
  PoolPatchSeed %in% figure3$graph$seed,
  Time == figure3$graph$time
) |> tidytable::distinct(
)

figure3$graph$networks <- generateNetworks(figure3$graph$specification,
                                           Date = "2025-07-30", split = FALSE)

# Main Plots: #################################################################
### Plot 3: ###################################################################
##### Data: ###################################################################

# Interventions store the time right before intervention, then the
# time of the intervention itself. Retrieve this second time.
# Note it is per-simulation (timescale set by PoolPatch effectively).
# (If done correctly, the set-up as of 23/01/2026 means there will be
# round numbers from 1:10, evens from 12:20, and then by 3s 20:50.)
figure3$interventionTimes <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == figure3$pref,
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

figure3$dataBase <- tidytable::bind_rows(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure3$pref,
    NicheDistance == defaultNicheDistance,
    InterventionInitial == figure3$luinitl,
    InterventionFinal %in% figure3$lufinal,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0"
    # Need all subsets
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure3$pref,
    NicheDistance == defaultNicheDistance,
    InterventionInitial == figure3$luinitl,
    InterventionFinal %in% figure3$lufinal,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance"
    # Need all subsets
  )
) |> tidytable::left_join(
  figure3$interventionTimes |> tidytable::rename(
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
figure3$dataOverallSummary <- figure3$dataBase |> tidytable::filter(
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
  Lower = quantile(Value, p = (1 - figure3$CI) - (1 - figure3$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure3$CI + (1 - figure3$CI)/2, na.rm = TRUE)
)

# Same idea as the overall case, but split by guild.
figure3$dataBCSupplement <- figure3$dataBase |> tidytable::filter(
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
  Lower = quantile(Value, p = (1 - figure3$CI) - (1 - figure3$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure3$CI + (1 - figure3$CI)/2, na.rm = TRUE)
)

figure3$indices <- figure3$graph$networks$Index |> tidytable::filter(
  SpeciesPreferences == figure3$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c(figure3$luinitl,
                      paste(figure3$luinitl, figure3$lufinal, sep = "->")),
  PoolPatchSeed %in% figure3$graph$seed
) |> tidytable::arrange(
  Intervention
)

##### A: ################################################################
figure3$plotA <- plotMeanAndInner(
  figure3$dataBase |> tidytable::filter(
    Metric == "Richness",
    is.na(Subset) # Not overall values
  ) |> tidytable::mutate(
    Time = Time + InterventionTime
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_vline(
  xintercept = min(figure3$interventionTimes$Time),
  color = "black", linetype = "dashed"
) + ggplot2::labs(
  x = "Time",
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",#,
  # color = "none",
  # fill = "none"
  fill = ggplot2::guide_legend(title = "Habitat Type",
                               override.aes = list(alpha = 1))
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax),
  xlim = c(0, 31000),
  expand = FALSE
)

##### Networks: ###############################################################
# Example networks from different scenarios of the same simulation, showing
# effects of the current habitat type through time on network shape.
# Previously, these were independent panels, but I'm switching to a facets.
figure3$plotNetworks <- figure3$graph$networks$Plot + ggplot2::facet_grid(
  # Reverse order
  factor(Intervention,
         levels = c("(0.5)->(1)", "(0.5)", "(0.5)->(0)"), ordered = T) ~ .
) + ggplot2::theme(
  axis.title.x = ggplot2::element_blank(),
  panel.border = ggplot2::element_rect(color = "black", fill = NA),
  legend.position = "none",
  axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = -30),
                                      size = 12),
  axis.text.x = ggplot2::element_blank()
) + ggplot2::coord_cartesian(xlim = c(-0.25, 1))

##### B: ######################################################################
figure3$plotB <- ggplot2::ggplot(
  figure3$dataOverallSummary |> tidytable::filter(
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
) + ggplot2::facet_grid( # For consistent plot sizing.
  "" ~ .
)

##### C+D: ####################################################################
figure3$plotCD <- ggplot2::ggplot(
  figure3$dataBCSupplement |> tidytable::filter(
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
figure3$plotE <- ggplot2::ggplot(
  figure3$dataOverallSummary |> tidytable::filter(
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
) + ggplot2::facet_grid(
  "" ~ .
)

##### F+G: ####################################################################
figure3$plotFG <- ggplot2::ggplot(
  figure3$dataBCSupplement |> tidytable::filter(
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
# figure3$plot <- ggpubr::ggarrange(
#   plotlist = list(
#     figure3$plotA + ggplot2::scale_y_continuous(
#       breaks = c(0, 10, 20, 30, 40),
#       labels = c("0", "10", "20", "30", "")
#     ),
#     ggpubr::ggarrange(plotlist = list(
#       figure3$plotB + ggplot2::theme(legend.position = "none"),
#       figure3$plotCD + ggplot2::theme(legend.position = "none",
#                                       axis.title.y = ggplot2::element_blank())
#     ), nrow = 1),
#     ggpubr::ggarrange(plotlist = list(
#       figure3$plotE + ggplot2::theme(legend.position = "none",
#                                      axis.text.y = ggplot2::element_text(
#                                        margin = ggplot2::margin(0, 0, 0, -5)
#                                      )),
#       figure3$plotFG + ggplot2::theme(legend.position = "none",
#                                       axis.title.y = ggplot2::element_blank(),
#                                       axis.text.y = ggplot2::element_text(
#                                         margin = ggplot2::margin(0, 0, 0, -5)
#                                       ))
#     ), nrow = 1)
#   ), nrow = 3, common.legend = TRUE
# )

figure3$plot <-
  figure3$plotA + ggplot2::scale_y_continuous(
    breaks = c(0, 10, 20, 30, 40),
    labels = c("0", "10", "20", "30", "")
  ) +
  figure3$plotB + ggplot2::theme(
    legend.position = "none"
  ) +
  figure3$plotCD + ggplot2::theme(
    legend.position = "none",
    axis.title.y = ggplot2::element_blank(),
    panel.spacing = ggplot2::unit(0.8, "lines")
  ) +
  figure3$plotE + ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(
      margin = ggplot2::margin(0, 0, 0, -5)
    )
  ) +
  figure3$plotFG + ggplot2::theme(
    legend.position = "none",
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(
      margin = ggplot2::margin(0, 0, 0, -5)
    ),
    panel.spacing = ggplot2::unit(0.8, "lines")
  ) +
  patchwork::plot_layout(design = "
    #111111#
    #111111#
    #111111#
    22222333
    22222333
    44444555
    44444555
  ")


figure3$suffix <- paste0("_", figure3$prefstring, "_", figure3$lustring)
figure3$prefix <- "Figure3"
figure3$iter <- "_Prototype6"
figure3$ids <- c("", "", "A", "Networks", "B", "CD", "E", "FG")
figure3$ext <- c(".png", rep(".pdf", 7))

for (fnum in 1:8) {
  with(figure3, ggplot2::ggsave(
    plot = switch(fnum,
                  plot, plot,
                  plotA, plotNetworks, plotB, plotCD, plotE, plotFG),
    filename = file.path(dirImages,
                         paste0(prefix, ids[fnum], iter, suffix, ext[fnum])),
    units = "cm", width = 6.5*3, height = 6.5*2.8)
  )
}
