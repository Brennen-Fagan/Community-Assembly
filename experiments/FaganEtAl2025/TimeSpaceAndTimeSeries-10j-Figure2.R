# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use,
# when species preferences are 100% 0.
# Also functinally an overview plot of network structure.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")
source(file.path("R", "flattenDiversity.R")) # Req'd by below
source(file.path("R", "generateNetworks.R")) # To create inset graphs.

# This is better as an environment, but that's more opaque.
figure2 <- list(
  graph = list(
    seed = "2", # "11", "17", "2"!,
    time = 25000
  ),
  pref = "100% 0"#"Uniform(0, 1)"
)

figure2$graph$specification <- diversitiesRichness |> tidytable::select(c(
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
  "SpeciesPreferences", "Intervention"
)) |> tidytable::filter(
  SpeciesPreferences == figure2$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  PoolPatchSeed %in% figure2$graph$seed,
  Time == figure2$graph$time
) |> tidytable::distinct(
)

figure2$graph$networks <- generateNetworks(figure2$graph$specification,
                                           Date = "2025-07-30", split = FALSE)

# Main Plots: #################################################################
### Plot 2: ###################################################################
##### Data: ###################################################################
figure2$data <- dplyr::full_join(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure2$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) |> tidytable::pivot_wider(
    names_from = Metric, values_from = Value
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure2$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance",
    is.na(Subset)
  ) |> tidytable::pivot_wider(
    names_from = Metric, values_from = Value
  ),
  by = dplyr::join_by(
    Time, Environment1, Environment2, Subset, PoolPatch, PoolPatchSeed,
    Interactions, InteractionsSeed, Events, EventsSeed,
    InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
    SpeciesAffinity, SpeciesAffinitySeed, PatchAffinity, PatchAffinitySeed,
    InterventionPatchType, InterventionPatchSeed, InterventionTimeType,
    InterventionTimeSeed, InterventionDispersal, InterventionNicheDistance,
    Intervention, SpeciesPreferences, InterventionInitial, InterventionFinal,
    DispersalParam, Start, Stop)
)

figure2$indices <- figure2$graph$networks$Index |> tidytable::filter(
  SpeciesPreferences == figure2$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds
) |> tidytable::arrange(
  Intervention
)

##### a: ######################################################################
# Richness through time across simulations, showing stability and separation.
figure2$plotA <- plotMeanAndInner(
  rbind(
    figure2$data |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    figure2$data |> tidytable::filter(
      PoolPatchSeed == figure2$graph$seed,
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - figure2$graph$time) == min(abs(Time - figure2$graph$time))
    ) |> tidytable::mutate(
      `Alpha Hill:0` = -100 # coord_cartesian will eliminate these points.
    )
  ) |> tidytable::rename(
    Value = `Alpha Hill:0`
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = function(x) {x |> tidytable::filter(
    PoolPatchSeed == figure2$graph$seed,
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    abs(Time - figure2$graph$time) == min(abs(Time - figure2$graph$time))
  )},
  mapping = ggplot2::aes(fill = Intervention),
  shape = 21,
  color = "black"
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5, title = "Habitat Type"),
  fill = ggplot2::guide_legend(ncol = 5, title = "Habitat Type")
) + ggplot2::coord_cartesian(
  xlim = c(0, 30000), ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::theme(
  legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95),
  axis.text.x = ggplot2::element_text(hjust = 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
# Example networks from different scenarios of the same simulation, showing
# effects of the current habitat type through time on network shape.
# Previously, these were independent panels, but I'm switching to a facets.
figure2$plotB <- figure2$graph$networks$Plot + ggplot2::facet_grid(
  # Reverse order
  factor(Intervention, levels = c("(1)", "(0.5)", "(0)"), ordered = T) ~ .
  ) + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", fill = NA)
  )

##### c: ######################################################################
# Richness varies with land-use type for our fixed land-use preference (0).
figure2$plotC <- ggplot2::ggplot(
  figure2$data |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) |> tidytable::summarise(
    Value = mean(`Alpha Hill:0`)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 1,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
# ) + ggplot2::geom_line( # Tracks the mean across sims and habitat types.
#   data = ~ summarise(group_by(.x, Intervention), Value = mean(Value)),
#   color = "black", group = 1
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::annotate(
  "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted"),
  size = 3
)

##### d: ######################################################################
# Abundance has a complex relationship with land-use type for fixed preference.
figure2$plotD <- ggplot2::ggplot(
  figure2$data |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) |> tidytable::summarise(
    Value = mean(`Alpha Abundance`)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 1,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
# ) + ggplot2::geom_line( # Tracks the mean across sims and habitat types.
#   data = ~ summarise(group_by(.x, Intervention), Value = mean(Value)),
#   color = "black", group = 1
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  y = "Avg. Total Abundance",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
# ) + ggplot2::annotate(
#   "text", x = c(1.5, 4.5), y = 3e3, label = c("Well\nAdapted", "Poorly\nAdapted"),
#   size = 3
) + ggplot2::scale_y_log10(
)

##### e: ######################################################################
# Richness and abundance co-vary for our scenarios.
figure2$plotE <- ggplot2::ggplot(
  figure2$data |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by( # Reduce to per run (x44 sims for param combns)
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) |> tidytable::summarise(
    `Alpha Hill:0` = mean(`Alpha Hill:0`),
    `Alpha Abundance` = mean(`Alpha Abundance`)
  ),
  ggplot2::aes(
    x = `Alpha Abundance`,
    y = `Alpha Hill:0`,
    color = Intervention
  )
) + ggplot2::geom_point(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  x = "Avg. Total Abundance",
  y = "Avg. Richness"
) + ggplot2::guides(
  color = "none",
  fill = "none"
# ) + ggplot2::scale_x_log10(
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax)
)


##### Combine: ################################################################
figure2$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure2$plotA,
    figure2$plotB + ggplot2::theme(legend.position = "none"),
    ggpubr::ggarrange(plotlist = list(
      figure2$plotC,
      figure2$plotD,
      figure2$plotE
    ), ncol = 1)
  ), nrow = 1, widths = c(2/4, 1/4, 1/4), common.legend = TRUE
)

ggplot2::ggsave(plot = figure2$plot, filename = file.path(dirImages, "Figure2_Prototype9.pdf"),
                units = "cm", width = 6.5*4, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plot, filename = file.path(dirImages, "Figure2_Prototype9.png"),
                units = "cm", width = 6.5*4, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotA, filename = file.path(dirImages, "Figure2A_Prototype9.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*3)
ggplot2::ggsave(plot = figure2$plotB, filename = file.path(dirImages, "Figure2B_Prototype9.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotC, filename = file.path(dirImages, "Figure2C_Prototype9.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotD, filename = file.path(dirImages, "Figure2D_Prototype9.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotE, filename = file.path(dirImages, "Figure2E_Prototype9.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
