old4 <- list()

old4$newplot4_dataA <- diversitiesRichness |> tidytable::filter(
  # SpeciesPreferences == "100% 0",
  # SpeciesPreferences == "50% 0, 50% 1",
  SpeciesPreferences == "Uniform(0, 1)",
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(1:44)),
  Metric == "Alpha Hill:0",
  is.na(Subset),
  (InterventionInitial == "(0.5)") #| (Time > 20000 & Time < 30000))
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

# old4$newplot4_dataB <- diversitiesRichness |> tidytable::filter(
#   NicheDistance == defaultNicheDistance,
#   (PoolPatchSeed %in% as.character(1:44)),
#   Metric == "Alpha Hill:0",
#   InterventionInitial == "(0.5)",
#   SpeciesPreferences == "100% 0",
#   !is.na(Subset)
# ) |> tidytable::left_join(
#   endTimes |> dplyr::select(-Times)
# ) |> tidytable::group_by(
#   SpeciesPreferences, Intervention, PoolPatchSeed,
#   InterventionInitial, InterventionFinal, Subset
# ) |> tidytable::arrange(
#   Time
# ) |> tidytable::filter(
#   InterventionInitial != InterventionFinal,
#   Time == Time[1] | Time == Time[2]
# ) |> tidytable::summarise(
#   Time = Time[2] - Time[1],
#   Value = Value[2] - Value[1],
#   Method = "Temporal",
#   Weight = 1, # for loess in geom_smooth
#   .groups = "drop"
# ) |> tidytable::right_join(
#   tidytable::expand(
#     diversitiesRichness,
#     tidytable::nesting(
#       SpeciesPreferences, Intervention, # SpeciesAffinity not working???
#       InterventionInitial, InterventionFinal,
#       Subset
#     )
#   ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
# ) |> tidytable::mutate(
#   Time = ifelse(is.na(Time), 0, Time),
#   Value = ifelse(is.na(Value), 0, Value),
#   Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
#   SpeciesPreferences = ifelse(is.na(SpeciesPreferences), "100% 0", SpeciesPreferences)
# ) |> tidytable::filter(
#   InterventionInitial == "(0.5)"
# )

old4$newplot4_dataB <- diversitiesRichness |> tidytable::filter( # taken from current 4
  # SpeciesPreferences == "100% 0",
  # SpeciesPreferences == "50% 0, 50% 1",
  SpeciesPreferences == "Uniform(0, 1)",
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  !is.na(Subset), # Basals and Consumers
  InterventionInitial == "(0.5)",
  InterventionFinal %in% c("(0)", "(0.5)", "(1)")
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal#,
  # Time %in% c(Time[1:140]) # 30 is 1:10, 12:20:2, 23:50:3, 100 gets to 600
) |> tidytable::summarise(
  Time = Time - Time[2],
  Value = Value - Value[2],
  Method = "Temporal",
  .groups = "drop"
  ) |> tidytable::filter(
    Time <= 15
) |> tidytable::mutate(
  Weight = ifelse(Time < 1e-6, 1e9, 1), # loess in geom_smooth to anchor to 0.
  Alpha = ifelse(Time <= 10, 0.1, 0)
)

##### a: ######################################################################
old4$newplot4_a <- plotMeanAndInner(
  rbind(
    old4$newplot4_dataA |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)->(0)", "(0.5)", "(0.5)->(1)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    old4$newplot4_dataA |> tidytable::filter(
      PoolPatchSeed == 1,
      Intervention %in% c("(0.5)->(0.25)", "(0.5)->(0.75)"),
      abs(Time - 20000) == min(abs(Time - 20000))
    ) |> tidytable::mutate(
      Value = -100
    )
  ),
  CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 7),
  fill = ggplot2::guide_legend(ncol = 7)
) + ggplot2::coord_cartesian(
  xlim = c(0, 30500), ylim = c(0, 42), expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 16300,
  xmax = 16500,
  ymin = 0, ymax = 42,
  fill = "grey",
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  # legend.position = c(0.35, 0.09),
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
  # ) + ggforce::facet_zoom(
  #   xlim = c(16000, 17000),
  #   shrink = FALSE
)

##### b: ######################################################################
# Could probably add in the counterfactual as well, but might be too messy?
old4$newplot4_b <- ggplot2::ggplot(
  old4$newplot4_dataB,
  aes(x = Time, y = Value,
      group = interaction(SpeciesPreferences, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::theme(
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

old4$colorPalette_0.5 <- colorPalette[
  grepl(x = names(colorPalette), pattern = "^[(]0.5[])]")
  ]
old4$newplot4_b_smooths <- ggplot2::ggplot_build(
  old4$newplot4_b
)$data[[2]] |> dplyr::group_by(
  group
) |> dplyr::filter(
  x == max(x)
) |> dplyr::ungroup(
) |> dplyr::mutate(
  Subset = rev(levels(factor(old4$newplot4_dataB$Subset)))[PANEL],
  Intervention = names(old4$colorPalette_0.5)[
    apply(outer(colour, old4$colorPalette_0.5, `==`), 1, which)
    ],
  yshift = y + c(-1, +2.5)#, -2, -4, -1, -3, -2, +2.5)
)

old4$newplot4_b <-
  # Witchcraft from stackoverflow.com/a/6675163
  # works by pre-building the plot and then extracting coordinates.
  old4$newplot4_b + ggplot2::coord_cartesian(
    xlim = c(0, 14), clip = "off"
  ) + ggplot2::geom_segment(
    data = old4$newplot4_b_smooths,
    mapping = ggplot2::aes(x = x+1, y = yshift, xend = x, yend = y,
                           color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + ggplot2::geom_label(
    data = old4$newplot4_b_smooths,
    mapping = ggplot2::aes(x = x+2.5, y = yshift,
                           label = Intervention, color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  )

old4$newplot4 <- ggpubr::ggarrange(
  plotlist = list(
    old4$newplot4_a,
    old4$newplot4_b
  ), nrow = 1, common.legend = TRUE #, widths = c(0.5, 0.27, 0.23)
) + ggplot2::annotate(
  "curve", x = 0.29, y = 0.1, xend = 0.5, yend = 0.15,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
)

ggplot2::ggsave(plot = old4$newplot4, 
                # filename = "Figure4_Prototype2-1000.pdf",
                filename = "Figure4_Prototype2-unif.pdf",
                units = "cm", width = 6.5*3, height = 6.5*2)
