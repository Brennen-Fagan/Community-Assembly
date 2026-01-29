# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use,
# when species preferences are 100% 0.
# Also functinally an overview plot of network structure.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")
library(ggbreak)

# This is better as an environment, but that's more opaque.
figure6 <- list(
  pref = c("Uniform(0, 1)", "50% 0, 50% 1")
)

# Main Plots: #################################################################
### Plot 6: ###################################################################
##### Data: ###################################################################
# Richness data: should be straightforward.
figure6$dataRich <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences %in% figure6$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset)
)

# Persistence data: why? Because we're using persistence as a weight, followed
# with by species aggregation. That way approximately we are picking a random
# simulation, a random time, and then a random species, the plot shows the
# probability we would get a certain land-use preference out.
figure6$dataPers <- Pers |> tidytable::filter(
  SpeciesPreferences %in% figure6$pref,
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
# Richness through time across simulations, showing stability and separation.
figure6$plotA <- plotMeanAndInner(
  rbind(
    figure6$dataRich |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    figure6$dataRich |> tidytable::filter(
      PoolPatchSeed == figure6$dataRich$PoolPatchSeed[1],
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - 10000) == min(abs(Time - 10000))
    ) |> tidytable::mutate(
      Value = 10 # coord_cartesian will eliminate these points.
    )
  ), CIs = 0.75, facets = as.formula(SpeciesPreferences ~ .)
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::theme(
  legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
) + ggbreak::scale_x_break(
  c(5000, 20000), expand = FALSE
) + ggplot2::scale_y_continuous(
  ylim = c(0, richnessYMax)
)

##### b: Violins ##############################################################
figure6$plotB <- ggplot2::ggplot(
  figure6$dataRich |> tidytable::filter(
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
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9), scale = "count"
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesPreferences ~ .
)

figure6$plotC <- ggplot2::ggplot(
  figure6$dataRich |> tidytable::filter(
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
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9), scale = "count"
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesPreferences ~ .
)

##### b: Insets ###############################################################
figure6$insetB <- ggplot2::ggplot(
  figure6$dataPers |> tidytable::filter(SpeciesPreferences == "50% 0, 50% 1"),
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
  name = "Habitat Type"
) + ggplot2::theme_void(
) + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
) + ggplot2::coord_cartesian(
  expand = FALSE
)
figure6$insetC <- ggplot2::ggplot(
  figure6$dataPers |> tidytable::filter(SpeciesPreferences == "Uniform(0, 1)"),
  ggplot2::aes(
    x = Affinity,
    weight = Persistence,
    fill = Intervention,
    group = Intervention
  )
) + ggplot2::geom_density(
  adjust = 1/2,
  show.legend = FALSE
) + ggplot2::facet_grid(
  . ~ Intervention
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_void(
) + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
) + ggplot2::coord_cartesian(
  expand = FALSE
)

##### Combine: ################################################################
figure6$plot <- ggpubr::ggarrange(
  plotlist = list(
    print(figure6$plotA), # github.com/YuLab-SMU/ggbreak/issues/36 to fix.
    ggpubr::ggarrange(
      plotlist = list(
        figure6$plotB + ggplot2::annotation_custom(
          ggplot2::ggplotGrob(figure6$insetB),
          xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
        ),
        figure6$plotC + ggplot2::annotation_custom(
          ggplot2::ggplotGrob(figure6$insetC),
          xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
        )
      ), ncol = 1
    )
  ),
  nrow = 1, common.legend = TRUE
)

ggplot2::ggsave(plot = figure6$plot, filename = file.path(dirImages, "Figure6_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure6$plot, filename = file.path(dirImages, "Figure6_Prototype1.png"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure6$plotA, filename = file.path(dirImages, "Figure6A_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure6$plotB, filename = file.path(dirImages, "Figure6B_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure6$plotC, filename = file.path(dirImages, "Figure6C_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure6$insetB, filename = file.path(dirImages, "Figure6BInset_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure6$insetC, filename = file.path(dirImages, "Figure6CInset_Prototype1.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
