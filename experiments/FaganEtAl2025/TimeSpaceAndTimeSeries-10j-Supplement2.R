### 2 Supplement: #############################################################
##### a: ######################################################################
newplot2_as <- plotMeanAndInner(
  rbind(
    newplot2_dataAS |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    newplot2_dataAS |> tidytable::filter(
      PoolPatchSeed == newplot2_a_seed,
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time)),
      !is.na(Subset)
    ) |> tidytable::mutate(
      Value = -100
    )
  ), CIs = 0.75, facets = as.formula(
    factor(Subset, levels = c("Consumer_0", "Basal_0"), ordered = TRUE) ~ .
  )
) + ggplot2::geom_point(
  data = newplot2_dataAS |> tidytable::filter(
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
  xlim = c(0, 31000), ylim = c(0, 42), expand = FALSE
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
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.5, 0.9),
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
newplot2_bs <- ggplot2::ggplot(
  newplot2_dataAS |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity, Subset
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
  ylim = c(0, 42), expand = FALSE
  # ) + ggplot2::annotate(
  #   "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted")
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"), ordered = TRUE) ~ .
)

newplot2s <- ggpubr::ggarrange(
  plotlist = list(
    newplot2_as,
    newplot2_bs
  ), nrow = 1, widths = c(0.5, 0.4)
)

ggplot2::ggsave(plot = newplot2s, filename = "Figure2s1_Prototype2.png",
                units = "cm", width = 6.5*3, height = 6.5*2)