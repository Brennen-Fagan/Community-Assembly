
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

figure3 <- list()

### Plot 3:####################################################################
figure3$dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences != "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(1:44)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

figure3$dataB <- Pers |> tidytable::filter(
  SpeciesPreferences != "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(1:44))
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
newplot3_a <- ggplot2::ggplot(
  newplot3_dataA |> tidytable::filter(
    Time > Start, Time < Stop, SpeciesAffinity == "50% 0, 50% 1"
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
  position = ggplot2::position_dodge(0.9), scale = "count"
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesAffinity ~ .
)
newplot3_b <- ggplot2::ggplot(
  newplot3_dataA |> tidytable::filter(
    Time > Start, Time < Stop, SpeciesAffinity == "Uniform(0, 1)"
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
  position = ggplot2::position_dodge(0.9), scale = "count"
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesAffinity ~ .
)

# Iteration 10 will have the actual affinities, rather than the affinitybins
# available, which will allow us to look at P/CDFs rather than bar charts.
# It still lets us think about how to quantify the distribution of affinities.
# I'm thinking Persistence as a weight, then by species aggregation. That way
# we get something like if I pick a random simulation, a random time, and then
# a random species, the plot shows the probability we would get a certain
# land-use preference out. (Weight by abundance as well for individuals, but
# that skews even more heavily towards basal species.)

newplot3_inset1 <- ggplot2::ggplot(
  newplot3_dataB |> tidytable::filter(SpeciesAffinity == "50% 0, 50% 1"),
  ggplot2::aes(
    x = AffinityBins,
    weight = Persistence,
    fill = Intervention
  )
) + ggplot2::geom_bar(
  show.legend = FALSE
) + ggplot2::facet_grid(
  . ~ Intervention
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_void(
) + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
) + ggplot2::coord_cartesian(
  expand = FALSE
)
newplot3_inset2 <- ggplot2::ggplot(
  newplot3_dataB |> tidytable::filter(SpeciesAffinity == "Uniform(0, 1)"),
  ggplot2::aes(
    x = AffinityBins,
    weight = Persistence,
    fill = Intervention
  )
) + ggplot2::geom_bar(
  show.legend = FALSE
) + ggplot2::facet_grid(
  . ~ Intervention
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_void(
) + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
) + ggplot2::coord_cartesian(
  expand = FALSE
)

newplot3 <- ggpubr::ggarrange(
  plotlist = list(
    newplot3_a + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(newplot3_inset1),
      xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
    ),
    newplot3_b + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(newplot3_inset2),
      xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
    )
  ), nrow = 1
)

ggplot2::ggsave(plot = newplot3, filename = "Figure3_Prototype6.png",
                units = "cm", width = 6.5*3, height = 6.5*2)
