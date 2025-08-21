# Setup: ######################################################################
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source(file.path("R", "generateNetworks.R")) # To create inset graphs.

# This is better as an environment, but that's more opaque.
figure2 <- list(
  graph = list(
    seed = "1",
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
  # "SpeciesAffinity",
  "SpeciesAffinitySeed", "PatchAffinity", "PatchAffinitySeed",
  # Which File (Intervention):
  "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
  "InterventionTimeSeed", "InterventionDispersal", "InterventionNicheDistance",
  # Ease of Use
  "SpeciesAffinity", "Intervention"
)) |> tidytable::filter(
  SpeciesAffinity == "100% 0" &
    NicheDistance == defaultNicheDistance &
    Intervention %in% c("(0)", "(0.5)", "(1)") &
    PoolPatchSeed %in% as.character(figure2$graph$seed) &
    Time == figure2$graph$time
) |> tidytable::distinct(
)

figure2$graph$exampleNetworks <- generateNetworks(figure2$graph$specification)

# Main Plots: #################################################################
### Plot 2:####################################################################
# a=>b&c
newplot2_dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot2_dataAS <- diversitiesRichness |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  !is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot2_indices <- exampleNetworks$Index |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)

newplot2_dataC <- Pers |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
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
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, AffinityBins,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, DispersalParam, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) |> tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
)

##### a: ######################################################################
newplot2_a <- plotMeanAndInner(
  rbind(
    newplot2_dataA |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    newplot2_dataA |> tidytable::filter(
      PoolPatchSeed == newplot2_a_seed,
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time))
    ) |> tidytable::mutate(
      Value = -100
    )
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = newplot2_dataA |> tidytable::filter(
    PoolPatchSeed == newplot2_a_seed,
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time))
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), ylim = c(0, 42), expand = FALSE
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    exampleNetworks$Envs[[newplot2_indices$ID[1]]]$singletonGraphs[[1]] +
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
    exampleNetworks$Envs[[newplot2_indices$ID[2]]]$singletonGraphs[[1]] +
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
    exampleNetworks$Envs[[newplot2_indices$ID[3]]]$singletonGraphs[[1]] +
      ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      ) + ggplot2::coord_cartesian(
        xlim = c(NA, NA), ylim = c(-2, 0.5)
      ) # Easiest to probably just not worry about comparing between.
  ),
  xmin = 30500, xmax = 40000, ymin = 29, ymax = 39
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(newplot2_dataA$Start),
  xmax = max(newplot2_dataA$Stop),
  ymin = 0, ymax = max(newplot2_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_segment(
  data = data.frame(
    x = newplot2_a_time+250,
    y = c(10, 22, 39),
    xend = 30500,
    yend = c(11, 22, 36),
    Intervention = c("(0)", "(0.5)", "(1)")
  ),
  mapping = ggplot2::aes(
    x = x, y = y, xend = xend, yend = yend, color = Intervention
  ),
  inherit.aes = FALSE,
  arrow = arrow(length = unit(0.03, "npc"))
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
) + ggplot2::annotate(
  "text", x = 36500, y = c(16, 38), size = 3, lineheight = 0.7,
  label = c("Fully\nAdapted", "Poorly\nAdapted")
)

##### b: ######################################################################
newplot2_b <- ggplot2::ggplot(
  newplot2_dataA |> tidytable::filter(
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
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = 0, ymax = max(newplot2_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
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
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 0.95)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::annotate(
  "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted")
)

##### c: ######################################################################
newplot2_c <- ggplot2::ggplot(
  newplot2_dataC,
  ggplot2::aes(
    y = Persistence,
    x = Intervention,
    color = Intervention,
    group = interaction(Intervention, SpeciesType),
    fill = SpeciesType
  )
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

newplot2 <- ggpubr::ggarrange(
  plotlist = list(
    newplot2_a,
    newplot2_b,
    newplot2_c
  ), nrow = 1, widths = c(0.5, 0.27, 0.23)
)

newplot2 <- newplot2 + ggplot2::annotate(
  "curve",
  x = 0.35, y = 0.97, xend = c(0.57, 0.87), yend = 0.97,
  curvature = -0.075,
  arrow = arrow(length = unit(0.03, "npc"))
)

ggplot2::ggsave(plot = newplot2, filename = "Figure2_Prototype7.png",
                units = "cm", width = 6.5*3, height = 6.5*2)
