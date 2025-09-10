# Setup: ######################################################################
# Figure 2 but for basal and consumer richness instead.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

supplement2 <- list()

### 2 Supplement: #############################################################
supplement2$dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  !is.na(Subset)
)

supplement2$dataB <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  !is.na(Subset)
)

supplement2$dataC <- with(supplement2, {
  # Average Abundance of species in the patch
  rbind(dataA, dataB) |> tidytable::filter(
    Start <= Time, Time <= Stop
  ) |> tidytable::pivot_wider(
    names_from = Metric, values_from = Value
  ) |> tidytable::mutate(
    `Average Abundance` =
      ifelse(`Alpha Hill:0` > 0, `Alpha Abundance` / `Alpha Hill:0`, 0)
  ) |> tidytable::group_by(
    # Average over relevant times
    Environment1:DispersalParam
  ) |> tidytable::summarise(
    `Average Abundance` = mean(`Average Abundance`)
  )
}) # 44 replicates, 5 treatments, 2 subsets.

##### a: ######################################################################
supplement2$plotA <- plotMeanAndInner(
  rbind(
    supplement2$dataA |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    supplement2$dataA |> tidytable::filter(
      PoolPatchSeed == "1",
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - 0) == min(abs(Time - 0)),
      !is.na(Subset)
    ) |> tidytable::mutate(
      Value = -100
    )
  ), CIs = 0.75, facets = as.formula(
    factor(Subset, levels = c("Consumer_0", "Basal_0"),
           labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 31000), ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(supplement2$dataA$Start),
  xmax = max(supplement2$dataA$Stop),
  ymin = 0, ymax = max(supplement2$dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.5, 0.9),
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
supplement2$plotB <- ggplot2::ggplot(
  supplement2$dataA |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences, Subset
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
  ymin = 0, ymax = richnessYMax,
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
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
  # ) + ggplot2::annotate(
  #   "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted")
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
)

##### c: ######################################################################
supplement2$plotC <- ggplot2::ggplot(
  # supplement2$dataB |> tidytable::filter(
  #   Time > Start, Time < Stop
  # ) |> tidytable::group_by(
  #   PoolPatchSeed, Intervention, SpeciesPreferences, Subset
  # ) |> tidytable::summarise(
  #   Value = mean(Value)
  # ),
  supplement2$dataC,
  ggplot2::aes(
    x = Intervention,
    # y = Value,
    y = `Average Abundance`,
    color = Intervention
  )
# ) + ggplot2::geom_rect(
#   data = data.frame(
#     1 # 1 rectangle per row, so dummy df to prevent overplotting
#   ),
#   xmin = 0,
#   xmax = 6,
#   # ymin = 0,
#   # ymax = max(supplement2$dataB$Value),
#   ymin = 0.01,
#   ymax = 4000,
#   fill = "grey",
#   alpha = 0.2,
#   inherit.aes = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28, alpha = 0.5
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "c)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1),
  panel.background = ggplot2::element_rect(
    fill = "grey90", color = NA
  )
) + ggplot2::labs(
  # y = "Avg. Total Abundance",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
  # ) + ggplot2::annotate(
  #   "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted")
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .,
  scales = "free_y"
) + ggplot2::scale_y_log10(
)

supplement2$plot <- ggpubr::ggarrange(
  plotlist = list(
    supplement2$plotA,
    supplement2$plotB,
    supplement2$plotC
  ), nrow = 1, widths = c(0.5, 0.4, 0.4)
)

ggplot2::ggsave(plot = supplement2$plot, filename = "Figure_supplement2_v3.pdf",
                units = "cm", width = 6.5*4, height = 6.5*2)
