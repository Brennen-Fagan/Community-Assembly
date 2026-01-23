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
figure4$lustring <- switch(
  figure4$luinitl,
  "(0)" = "0",
  "(0.5)" = "", # Base Case
  "(1)" = "1"
)

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

figure4$data <- tidytable::bind_rows(
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
  Time > -1000, Time < 15000
)

##### a: ######################################################################
# Summarised time-series plot for overall richness.
figure4$plotA <- plotMeanAndInner(
    figure4$data |> tidytable::filter(
      Metric == "Richness", is.na(Subset)
    ), CIs = 0.75, facets = as.formula(. ~ .),
    digits = -2 # Due to temporal resolution of the data, above -2 we can't 
    # see the trends (not aggregated enough) but at -2, we can't see the
    # connection (the downward jump, because too aggregated).
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(-1000, 15000), ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::theme(
  legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95),
  axis.text.x = ggplot2::element_text(hjust = 1)
)

##### b: ######################################################################
# HEATMAPS

##### c: ######################################################################
# Summarised time-series plot for short time scale richness ratio.
# Note this should be Consumer / Basal (because Consumers sometimes -> 0).

figure4$plotC <- plotMeanAndInner(
  figure4$data |> tidytable::filter(
    !is.na(Subset), # Not overall values
    Metric == "Richness",
    Time > -0.01, Time < 51, Time == round(Time) # Avoid singletons.
  ) |> tidytable::separate_wider_delim(
    delim = "_", cols = Subset, names = c("Trophic", "Pref")
  ) |> tidytable::pivot_wider(
    names_from = Trophic, values_from = Value
  ) |> tidytable::group_by(
    Intervention, InterventionInitial, InterventionFinal, 
    SpeciesPreferences, Pref
  ) |> tidytable::mutate(
    Value = Consumer/Basal,
    Subset = NA
  ), CIs = 0.75, facets = as.formula(. ~ Pref),
  digits = 0
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness (Cons./Basal)"
) + ggplot2::guides(
  linetype = "none"#,
  # color = ggplot2::guide_legend(ncol = 5),
  # fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::theme(
  # legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95),
  axis.text.x = ggplot2::element_text(hjust = 1)
)

##### d: ######################################################################
# Summarised time-series plot for short time scale abundance ratio.
# Note this should be Consumer / Basal (because Consumers sometimes -> 0).

figure4$plotD <- plotMeanAndInner(
  figure4$data |> tidytable::filter(
    !is.na(Subset), # Not overall values
    Metric == "Abundance",
    Time > -0.01, Time < 51, Time == round(Time) # Avoid singletons.
  ) |> tidytable::separate_wider_delim(
    delim = "_", cols = Subset, names = c("Trophic", "Pref")
  ) |> tidytable::pivot_wider(
    names_from = Trophic, values_from = Value
  ) |> tidytable::group_by(
    Intervention, InterventionInitial, InterventionFinal, 
    SpeciesPreferences, Pref
  ) |> tidytable::mutate(
    Value = Consumer/Basal,
    Subset = NA
  ), CIs = 0.75, facets = as.formula(. ~ Pref),
  digits = 0
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Abundance (Cons./Basal)"
) + ggplot2::guides(
  linetype = "none"#,
  # color = ggplot2::guide_legend(ncol = 5),
  # fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::theme(
  # legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95),
  axis.text.x = ggplot2::element_text(hjust = 1)
) + ggplot2::scale_y_log10(
)

##### SUPPLEMENT: #############################################################
# Short-long term transition for richness and abundance on log(1+Time) scale
# expecting to hit the ratios as well as the base values.

##### Combine: ################################################################
figure4$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure4$plotA,
    ggpubr::ggarrange(plotlist = list(
      figure4$plotC + ggplot2::theme(legend.position = "none"),
      figure4$plotD + ggplot2::theme(legend.position = "none")
    ), ncol = 1)
  ), nrow = 1, common.legend = TRUE
)

#4prototype5