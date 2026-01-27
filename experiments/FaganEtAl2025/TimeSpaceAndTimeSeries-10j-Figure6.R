# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use,
# when species preferences are 100% 0.
# Also functinally an overview plot of network structure.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
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

)

##### a: ######################################################################
# Richness through time across simulations, showing stability and separation.
figure6$plotA <- plotMeanAndInner(
  rbind(
    figure6$dataRich |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    figure6$data |> tidytable::filter(
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
  plot.tag.position = c(0.025, 0.95),
  axis.text.x = ggplot2::element_text(hjust = 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
) + ggbreak::scale_x_break(
  c(5000, 20000), expand = FALSE
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

