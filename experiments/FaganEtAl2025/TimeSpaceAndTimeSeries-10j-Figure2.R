# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use,
# when species preferences are 100% 0.
# Also functinally an overview plot of network structure.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")
source(file.path("R", "flattenDiversity.R")) # Req'd by below
source(file.path("R", "generateNetworks.R")) # To create inset graphs.

# This is better as an environment, but that's more opaque.
figure2 <- list(
  graph = list(
    seed = "2", # "11", "17", "2"!,
    time = 25000
  )
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
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  PoolPatchSeed %in% figure2$graph$seed,
  Time == figure2$graph$time
) |> tidytable::distinct(
)

figure2$graph$networks <- generateNetworks(figure2$graph$specification)

# Main Plots: #################################################################
### Plot 2:####################################################################
# a=>b&c
figure2$dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset)
)

figure2$indices <- figure2$graph$networks$Index |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds
) |> tidytable::arrange(
  Intervention
)

figure2$dataC <- Pers |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds
) |> tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) |> tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) |> tidytable::group_by(
  Species:Affinity, AffinityBins,
  PoolPatch:InterventionNicheDistance,
  Intervention, SpeciesPreferences, Start, Stop
) |> tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
)

##### a: ######################################################################
####### Core Plot: ############################################################
figure2$plotA <- plotMeanAndInner(
  rbind(
    figure2$dataA |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    figure2$dataA |> tidytable::filter(
      PoolPatchSeed == figure2$graph$seed,
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - figure2$graph$time) == min(abs(Time - figure2$graph$time))
    ) |> tidytable::mutate(
      Value = -100 # coord_cartesian will eliminate these points.
    )
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = function(x) {x |> tidytable::filter(
    PoolPatchSeed == figure2$graph$seed,
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    abs(Time - figure2$graph$time) == min(abs(Time - figure2$graph$time))
  )}
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), ylim = c(0, richnessYMax), expand = FALSE
  ####### Insets: #############################################################
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    figure2$graph$networks$Envs[[
      figure2$indices$ID[1]
      ]]$singletonGraphs[[1]] +
      ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      ) + ggplot2::coord_cartesian(
        xlim = c(NA, NA), ylim = c(-2, 0.5)
      ) # Easiest to probably just not worry about comparing between.
  ),
  xmin = 30500, xmax = 40000, ymin = 7, ymax = 17
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    figure2$graph$networks$Envs[[
      figure2$indices$ID[2]
      ]]$singletonGraphs[[1]] +
      ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      ) + ggplot2::coord_cartesian(
        xlim = c(NA, NA), ylim = c(-2, 0.5)
      ) # Easiest to probably just not worry about comparing between.
  ),
  xmin = 30500, xmax = 40000, ymin = 18, ymax = 28
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    figure2$graph$networks$Envs[[
      figure2$indices$ID[3]
      ]]$singletonGraphs[[1]] +
      ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      ) + ggplot2::coord_cartesian(
        xlim = c(NA, NA), ylim = c(-2, 0.5)
      ) # Easiest to probably just not worry about comparing between.
  ),
  xmin = 30500, xmax = 40000, ymin = 29, ymax = 39
  ####### Annotations: ########################################################
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(figure2$dataA$Start),
  xmax = max(figure2$dataA$Stop),
  ymin = 0, ymax = richnessYMax,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_segment( # Arrows to the boxes.
  data = data.frame(
    x = figure2$graph$time+250,
    y = c(12, 26, 35),
    xend = 30500,
    yend = c(11, 22, 36),
    Intervention = c("(0)", "(0.5)", "(1)")
  ),
  mapping = ggplot2::aes(
    x = x, y = y, xend = xend, yend = yend, color = Intervention
  ),
  inherit.aes = FALSE,
  arrow = arrow(length = unit(0.03, "npc"))
) + ggplot2::annotate(
  "text", x = 36500, y = c(16, 38), size = 3, lineheight = 0.7,
  label = c("Fully\nAdapted", "Poorly\nAdapted")
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
figure2$plotB <- ggplot2::ggplot(
  figure2$dataA |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
  ###### Background Annotation: ###############################################
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = 0, ymax = richnessYMax,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
  ####### Core Plot: ##########################################################
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
) + ggplot2::geom_line(
  data = ~ summarise(group_by(.x, Intervention), Value = mean(Value)),
  color = "black", group = 1
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "b)"
  ####### Annotations: ########################################################
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 0.95)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::annotate(
  "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted")
)

##### c: ######################################################################
figure2$plotC <- ggplot2::ggplot(
  figure2$dataC,
  ggplot2::aes(
    y = Persistence,
    x = Intervention,
    color = Intervention,
    group = interaction(Intervention, SpeciesType),
    fill = SpeciesType
  )
  ###### Background Annotation: ###############################################
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = -Inf, ymax = Inf,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
  ####### Core Plot: ##########################################################
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9), show.legend = FALSE,
  scale = "count", draw_quantiles = 0.5, linewidth = 1.3
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("darkgreen", "burlywood4")
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  tag = "c)", x = "Habitat"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 0.95)
) + ggplot2::facet_grid(
  factor(SpeciesType, levels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::scale_x_discrete(
  breaks = c("(0)", "(0.5)", "(1)")
) + ggplot2::labs(
  x = "Habitat's Land-use"
)

figure2$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure2$plotA,
    figure2$plotB,
    figure2$plotC
  ), nrow = 1, widths = c(0.5, 0.27, 0.23)
)

figure2$plot <- figure2$plot + ggplot2::annotate(
  "curve",
  x = 0.35, y = 0.97, xend = c(0.57, 0.87), yend = 0.97,
  curvature = -0.075,
  arrow = arrow(length = unit(0.03, "npc"))
)

ggplot2::ggsave(plot = figure2$plot, filename = "Figure2_Prototype7.png",
                units = "cm", width = 6.5*3, height = 6.5*2)
