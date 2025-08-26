# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use,
# when species preferences are 50% 50% or Uniform.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")

figure3 <- list()

# Main Plots: #################################################################
### Plot 3:####################################################################
figure3$dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences != "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset)
)

figure3$dataB <- Pers |> tidytable::filter(
  SpeciesPreferences != "100% 0",
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
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed,
  Affinity, AffinityBins,
  PoolPatch:InterventionNicheDistance,
  Intervention, SpeciesPreferences, Start, Stop
) |> tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
)

##### a: ######################################################################
figure3$plotA <- ggplot2::ggplot(
  figure3$dataA |> tidytable::filter(
    Time > Start, Time < Stop, SpeciesPreferences == "50% 0, 50% 1"
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
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
  ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesPreferences ~ .
)

figure3$plotB <- ggplot2::ggplot(
  figure3$dataA |> tidytable::filter(
    Time > Start, Time < Stop, SpeciesPreferences == "Uniform(0, 1)"
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
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
  ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesPreferences ~ .
)

# Iteration 10 will have the actual affinities, rather than the affinitybins
# available, which will allow us to look at P/CDFs rather than bar charts.
# It still lets us think about how to quantify the distribution of affinities.
# I'm thinking Persistence as a weight, then by species aggregation. That way
# we get something like if I pick a random simulation, a random time, and then
# a random species, the plot shows the probability we would get a certain
# land-use preference out. (Weight by abundance as well for individuals, but
# that skews even more heavily towards basal species.)

figure3$insetA <- ggplot2::ggplot(
  figure3$dataB |> tidytable::filter(SpeciesPreferences == "50% 0, 50% 1"),
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
figure3$insetB <- ggplot2::ggplot(
  figure3$dataB |> tidytable::filter(SpeciesPreferences == "Uniform(0, 1)"),
  ggplot2::aes(
    x = Affinity,
    weight = Persistence,
    fill = Intervention,
    group = Intervention
  )
  #   # ECDF's not very clear here, and we'd need to recompute to get it proper
  #   # (weights not included by default).
  # ) + ggplot2::stat_ecdf(
  #   geom = "step",
  # Density Plot: I think better?
) + ggplot2::geom_density(
  adjust = 1/2,
  # # Histogram: Okay
  # ) + ggplot2::geom_histogram(
  #   bins = 100
  # # Bar Chart: Okay
  # ) + ggplot2::geom_bar(
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

figure3$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure3$plotA + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(figure3$insetA),
      xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
    ),
    figure3$plotB + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(figure3$insetB),
      xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
    )
  ), nrow = 1
)

ggplot2::ggsave(plot = figure3$plot, filename = "Figure3_Prototype7.png",
                units = "cm", width = 6.5*3, height = 6.5*2)
