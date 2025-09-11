# Setup: ######################################################################
# Plot of long-term and short-term responses of species richness to land-use
# change, taking 0.5 as a base case for sake of argument.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

figure4 <- list()

# Main Plots: #################################################################
### Plot 4:####################################################################
figure4$dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset)
)

figure4$dataB <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  InterventionInitial == "(0.5)",
  InterventionFinal %in% c("(0)", "(0.5)", "(1)"),
  SpeciesPreferences == "100% 0",
  !is.na(Subset)
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal#,
  # Time %in% c(Time[1:140]) # 30 is 1:10, 12:20:2, 23:50:3, 100 gets to 600
) |> tidytable::summarise(
  Time = Time - Time[1],
  Value = Value - Value[1],
  Method = "Temporal",
  .groups = "drop"
  # ) |> tidytable::filter(
  #   Time <= 1100
) |> tidytable::mutate(
  Weight = ifelse(Time < 1e-6, 1e9, 1), # loess in geom_smooth to anchor to 0.
  Alpha = ifelse(Time <= 10, 0.1, 0)
)

##### a: ######################################################################
figure4$plotA <- plotMeanAndInner(
  rbind(
    figure4$dataA |> tidytable::filter(
      Intervention %in% c(
        #"(0)",
        "(0.5)->(0)", "(0.5)", "(0.5)->(1)"
        #"(1)"
      ),
      ifelse(Intervention == "(0.5)", Time <= 16300, TRUE)
    # ),
    # # We want to appear in the legend but not on the plot!
    # figure4$dataA |> tidytable::filter(
    #   PoolPatchSeed == "1",
    #   Intervention %in% c("(0.5)->(0.25)", "(0.5)->(0.75)"),
    #   abs(Time - 20000) == min(abs(Time - 20000))
    # ) |> tidytable::mutate(
    #   Value = -100
    )
  ),
  CIs = 0.75, facets = as.formula(. ~ .)
  ###### Annotation: ##########################################################
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 7),
  fill = ggplot2::guide_legend(ncol = 7)
) + ggplot2::coord_cartesian(
  xlim = c(0, 30500), ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 16300,
  xmax = 21300,
  ymin = 0, ymax = richnessYMax,
  fill = "grey",
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
figure4$plotB <- ggplot2::ggplot(
  figure4$dataA |> tidytable::filter(
    Time >= 16250, Time <= 21250,
    Intervention %in% c("(0.5)->(0)", "(0.5)", "(0.5)->(1)")
  ),
  ggplot2::aes(
    x = Time, y = Value,
    color = Intervention,
    group = interaction(PoolPatchSeed, Intervention),
    alpha = (PoolPatchSeed == "1")
  )
) + ggplot2::geom_line(
  show.legend = FALSE
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme(
  plot.tag.position = c(0.025, 1)
) + ggplot2::coord_cartesian(
  xlim = c(16500, 21000)
)

##### c: ######################################################################
figure4$plotC <- ggplot2::ggplot(
  figure4$dataB,
  aes(x = Time, y = Value,
      group = interaction(SpeciesPreferences, Intervention),
      color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE, alpha = 0.02
  # mapping = ggplot2::aes(alpha = Alpha)
) + ggplot2::geom_smooth(
  show.legend = FALSE,
  # ggplot2::aes(weight = Weight),
  n = 1000
  # method = "loess",
  # formula = "y~x"
  ###### Annotation: ##########################################################
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  tag = "c)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
  # ) + ggplot2::scale_alpha(
  #   range = c(0, 0.1)
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::theme(
  plot.background = ggplot2::element_rect(linetype = "solid")
) + ggplot2::coord_cartesian(
  xlim = c(0, 5000)
)

# # colorPalette was not being picked up correctly when providing all values.
# figure4$colorPalette <- colorPalette[
#   grepl(x = names(colorPalette), pattern = "^[(]0.5[])]")
#   ]
#
#
# ###### Labels: ################################################################
# # Witchcraft from stackoverflow.com/a/6675163
# # works by pre-building the plot and then extracting coordinates.
# figure4$smoothsC <- ggplot2::ggplot_build(
#   figure4$plotC
# )$data[[2]] |> dplyr::group_by(
#   group
# ) |> dplyr::filter(
#   x == x[min(abs(x - 5000)) == abs(x - 5000)]
# ) |> dplyr::ungroup(
# ) |> dplyr::mutate(
#   Subset = rev(levels(factor(figure4$dataB$Subset)))[PANEL],
#   Intervention = names(figure4$colorPalette)[
#     apply(outer(colour, figure4$colorPalette, `==`), 1, which)
#     ],
#   yshift = y + c(-1, +2.5, -2, -4, -1, -3, -2, +2.5) # DONE MANUALLY!!!
# )
#
# figure4$plotC <-
#   figure4$plotC + ggplot2::geom_segment(
#     data = figure4$smoothsC,
#     mapping = ggplot2::aes(x = x+1, y = yshift, xend = x, yend = y,
#                            color = Intervention),
#     show.legend = FALSE,
#     inherit.aes = FALSE
#   ) + ggplot2::geom_label(
#     data = figure4$smoothsC,
#     mapping = ggplot2::aes(x = x+2.5, y = yshift,
#                            label = Intervention, color = Intervention),
#     show.legend = FALSE,
#     inherit.aes = FALSE
#   )

figure4$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure4$plotA,
    figure4$plotB,
    figure4$plotC
  ),
  ncol = 1, common.legend = TRUE, heights = c(0.4, 0.2, 0.5)
  # nrow = 1, common.legend = TRUE #, widths = c(0.5, 0.27, 0.23)
# ) + ggplot2::annotate(
#   "curve", x = 0.29, y = 0.1, xend = 0.5, yend = 0.15,
#   arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
)

ggplot2::ggsave(plot = figure4$plot, filename = "Figure4_Prototype3.pdf",
                units = "cm", width = 6.5*3, height = 6.5*2)
