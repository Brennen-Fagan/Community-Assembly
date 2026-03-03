# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use, when
# species preferences are 100% 0. Also functinally an overview plot of network
# structure. This version combines aspects of figure 3.

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
  abundlog = TRUE,
  pref = "100% 0"#"Uniform(0, 1)"
)

# if (figure2$abundlog) {
#   figure2$abundLimits <- c(1e-1, 4e4)
# } else {
#   figure2$abundLimits <- c(0, 3.7e4)
# }

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

figure2$dataBC <- tidytable::bind_rows(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure2$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0",
    !is.na(Subset)
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure2$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance",
    !is.na(Subset)
  )
) |> tidytable::mutate(
  Subset = factor(Subset, levels = c("Consumer_0", "Basal_0"),
                  labels = c("Consumer", "Basal"), ordered = TRUE),
  Metric = factor(Metric, levels = c("Alpha Hill:0", "Alpha Abundance"),
                  labels = c("Richness", "Abundance"), ordered = TRUE)
) |> tidytable::filter(
  Time > Start, Time < Stop
) |> tidytable::group_by(
  PoolPatchSeed, Intervention, SpeciesAffinity, Metric, Subset
) |> tidytable::summarise(
  Value = mean(Value), .groups = "drop"
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
  xlim = c(0, 31000), ylim = c(0, richnessYMax), expand = FALSE
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
figure2$plotNetworks <- figure2$graph$networks$Plot + ggplot2::facet_grid(
  # Reverse order
  factor(Intervention, levels = c("(1)", "(0.5)", "(0)"), ordered = T) ~ .
) + ggplot2::theme(
  axis.title.x = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_blank(),
  panel.border = ggplot2::element_rect(color = "black", fill = NA)
)

##### c: ######################################################################
# Richness varies with land-use type for our fixed land-use preference (0).
figure2$plotB <- ggplot2::ggplot(
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
  y = "Richness",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank()
) + ggplot2::annotate(
  "text", x = c(1.5, 4.5), y = c(30, 15),
  label = c("Well\nAdapted", "Poorly\nAdapted"),
  size = 3
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
)

##### d: ######################################################################
# Abundance has a complex relationship with land-use type for fixed preference.
figure2$plotC <- ggplot2::ggplot(
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
  y = "Abundance",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
  # ) + ggplot2::annotate(
  #   "text", x = c(1.5, 4.5), y = 3e3, label = c("Well\nAdapted", "Poorly\nAdapted"),
  #   size = 3
) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank()
# ) + ggplot2::coord_cartesian(
#   ylim = figure2$abundLimits
)

if (figure2$abundlog) {
  figure2$plotC <- figure2$plotC + ggplot2::scale_y_log10(
    label = scales::label_log(digits = 2)
  )
}

##### f: ######################################################################
figure2$plotD <- ggplot2::ggplot(
  figure2$dataBC |> tidytable::filter(
    Metric == "Richness"
  ) |> tidytable::mutate( # Left-Right ordering, not Top-Bottom
    Subset = factor(Subset, ordered = TRUE, levels = c("Basal", "Consumer"))
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention,
    group = interaction(Intervention, Subset)
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 1,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::geom_text(
  data = function(x) {
    x |> tidytable::mutate(
      Offset = (max(Value) - min(Value))* 0.08
    ) |> tidytable::group_by(
      Intervention, Subset
    ) |> tidytable::summarise(
      Value = max(Value) + Offset,
      Label = substr(Subset, 0, 1),
      .groups = "drop"
    )
  },
  mapping = ggplot2::aes(label = Label),
  position = ggplot2::position_dodge(0.9)
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank()
) + ggplot2::labs(
  x = "Habitat Type", y = "Richness"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
)

##### g: ######################################################################
figure2$plotE <- ggplot2::ggplot(
  figure2$dataBC |> tidytable::filter(
    Metric == "Abundance"
  ) |> tidytable::mutate( # Left-Right ordering, not Top-Bottom
    Subset = factor(Subset, ordered = TRUE, levels = c("Basal", "Consumer"))
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention,
    group = interaction(Intervention, Subset)
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 1,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank()
) + ggplot2::labs(
  x = "Habitat Type", y = "Abundance"
# ) + ggplot2::coord_cartesian(
#   ylim = figure2$abundLimits
)

if (figure2$abundlog) {
  figure2$plotE <- figure2$plotE + ggplot2::scale_y_log10(
    label = scales::label_log(digits = 2)
  ) + ggplot2::geom_text(
    data = function(x) {
      x |> tidytable::mutate(
        Offset = (max(log10(Value)) - min(log10(Value))) * 0.08
      ) |> dplyr::group_by(
        Intervention, Subset
      ) |> dplyr::summarise(
        Value = ifelse(min(log10(Value)) >= 0,
                       10^min(log10(Value) - Offset),
                       10^max(log10(Value) + Offset)),
        Label = substr(Subset[1], 0, 1),
        .groups = "drop"
      )
    },
    mapping = ggplot2::aes(label = Label),
    position = ggplot2::position_dodge(0.9)
  )
} else {
  figure2$plotE <- figure2$plotE + ggplot2::geom_text(
    data = function(x) {
      x |> tidytable::mutate(
        Offset = (max(Value) - min(Value))* 0.08
      ) |> tidytable::group_by(
        Intervention, Subset
      ) |> tidytable::summarise(
        Value = max(Value) + Offset,
        Label = substr(Subset, 0, 1),
        .groups = "drop"
      )
    },
    mapping = ggplot2::aes(label = Label),
    position = ggplot2::position_dodge(0.9)
  )
}



##### Combine: ################################################################
figure2$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure2$plotA,
    figure2$plotNetworks + ggplot2::theme(legend.position = "none"),
    cowplot::plot_grid(plotlist = list(
      figure2$plotB, figure2$plotD,
      figure2$plotC, figure2$plotE
    ), ncol = 2)
  ), nrow = 1, widths = c(2/5, 1/5, 2/5), common.legend = TRUE
)
figure2$plotNoB <- ggpubr::ggarrange(
  plotlist = list(
    figure2$plotA,
    cowplot::plot_grid(plotlist = list(
      figure2$plotB, figure2$plotD,
      figure2$plotC, figure2$plotE
    ), ncol = 2)
  ), nrow = 1, widths = c(2/5, 3/5), common.legend = TRUE
)

ggplot2::ggsave(plot = figure2$plot,
                filename = file.path(dirImages, "Figure2_Prototype1.pdf"),
                units = "cm", width = 6.5*5, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plot,
                filename = file.path(dirImages, "Figure2_Prototype1.png"),
                units = "cm", width = 6.5*5, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotNoB,
                filename = file.path(dirImages, "Figure2noB_Prototype1.pdf"),
                units = "cm", width = 6.5*5, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotA,
                filename = file.path(dirImages, "Figure2A_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*3)
ggplot2::ggsave(plot = figure2$plotNetworks,
                filename = file.path(dirImages, "Figure2Networks_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotB,
                filename = file.path(dirImages, "Figure2B_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotC,
                filename = file.path(dirImages, "Figure2C_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotD,
                filename = file.path(dirImages, "Figure2D_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotE,
                filename = file.path(dirImages, "Figure2E_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)

