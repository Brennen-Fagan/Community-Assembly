# Setup: ######################################################################
# Figure 2 Equivalents, but for other species preferences and all
# intervention cases, which are handled via facetting to the ending type.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

# This is better as an environment, but that's more opaque.
supplement14 <- list()
supplement14$preferences <- 3
supplement14$yMax <- diversitiesRichness$Value |> max()

# Main Plots: #################################################################
### Plot 2:####################################################################
supplement14$data <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == switch(supplement14$preferences,
                               "100% 0",
                               "50% 0, 50% 1",
                               "Uniform(0, 1)"
  ),
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset)
)

##### a: ######################################################################
####### Core Plot: ############################################################
supplement14$plot <- plotMeanAndInner(
  supplement14$data, CIs = 0.75, facets = as.formula(. ~ InterventionFinal)
) + ggplot2::geom_line(
  data = function(x) x |> tidytable::filter(
    InterventionInitial == InterventionFinal
    ) |> tidytable::mutate(
      Time = round(Time, digits = -2)
    ) |> tidytable::group_by(
      Time, Subset,
      Intervention, InterventionInitial, InterventionFinal, SpeciesPreferences
    ) |> tidytable::summarise(
      Value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    ) |> tidytable::select(
      -InterventionFinal
    ), # Force to be plotted on each
  inherit.aes = FALSE,
  mapping = ggplot2::aes(
    x = Time, y = Value, color = Intervention
  ),
  alpha = 0.5
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 30000), ylim = c(0, supplement14$yMax), expand = TRUE
  ####### Annotations: ########################################################
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(supplement14$data$Start),
  xmax = max(supplement14$data$Stop),
  ymin = 0, ymax = supplement14$yMax,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::theme(
  legend.position = "right"
)

ggplot2::ggsave(plot = supplement14$plot,
                filename = paste0("Figure_supplement14_v1_",
                                  supplement14$preferences,
                                  ".pdf"),
                units = "cm", width = 6.5*9, height = 6.5*2)
