# Setup: ######################################################################
# Figure 2 Equivalents, but for other species preferences and all
# no-intervention cases.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

# This is better as an environment, but that's more opaque.
supplement13 <- list()
supplement13$preferences <- 1
supplement13$yMax <- diversitiesRichness$Value |> max()

# Main Plots: #################################################################
### Plot 2:####################################################################
supplement13$data <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == switch(supplement13$preferences,
    "100% 0",
    "50% 0, 50% 1",
    "Uniform(0, 1)"
  ),
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset)
)

##### a: ######################################################################
####### Core Plot: ############################################################
supplement13$plot <- plotMeanAndInner(
    supplement13$data, CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 30000), ylim = c(0, supplement13$yMax), expand = TRUE
  ####### Annotations: ########################################################
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(supplement13$data$Start),
  xmax = max(supplement13$data$Stop),
  ymin = 0, ymax = supplement13$yMax,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::theme(
  legend.position = "bottom"
)

ggplot2::ggsave(plot = supplement13$plot,
                filename = paste0("Figure_supplement13_v1_",
                                  supplement13$preferences,
                                  ".pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
