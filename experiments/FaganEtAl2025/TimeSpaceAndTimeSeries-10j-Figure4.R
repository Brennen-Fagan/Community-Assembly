### Plot 4:####################################################################
# Need to contrast with 2a (Richness). Long and short time scales.
# So starting from the top again, we want to construct the image to lead
# naturally to the comparison of the two different time spans for the same
# statistic, which probably looks more like differences of various sorts.
# Leaning on some of what we had originally set off to do, we might be able
# to characterise it as
#     time before-after short term,
#     counterfactual before-after short term,
#     time before-after long term (slope 0),
#     counterfactual before-after long term.
# Regardless, we probably need to re-run things so we have consistent
# comparisons that we are making, specifically for the time before-after short.
# So:
#  a => b + c, with d "contained within" a
#  d => e + f
# long time a: regular time; b: temporal comparison; c: counterfactual compare
# short time d: rescaled time; e: temporal compare; f: counterfactual compare
#
# Then it might be a good idea to summarise the b vs c and e vs f in the text
# via the differences.

newplot4_dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  # SpeciesAffinity == "50% 0, 50% 1",
  # SpeciesAffinity == "Uniform(0, 1)",
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot4_dataB <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  InterventionInitial == "(0.5)",
  SpeciesAffinity == "100% 0",
  !is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  Weight = 1, # for loess in geom_smooth
  .groups = "drop"
) |> tidytable::right_join(
  tidytable::expand(
    diversitiesRichness,
    tidytable::nesting(
      SpeciesAffinity, Intervention, # SpeciesAffinity not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
  ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesAffinity = ifelse(is.na(SpeciesAffinity), "100% 0", SpeciesAffinity)
) |> tidytable::filter(
  InterventionInitial == "(0.5)"
)

##### a: ######################################################################
newplot4_a <- plotMeanAndInner(
  rbind(
    newplot4_dataA |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)->(0)", "(0.5)", "(0.5)->(1)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    newplot4_dataA |> tidytable::filter(
      PoolPatchSeed == newplot2_a_seed,
      Intervention %in% c("(0.5)->(0.25)", "(0.5)->(0.75)"),
      abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time))
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

# ##### Temporal Vs. Counterfactual Statistics: #################################
# temporalVCounterfactualStats <- rbind(
#   # Temporal Substitution
#   diversitiesRichness |> tidytable::filter(
#     NicheDistance == defaultNicheDistance,
#     (PoolPatchSeed %in% as.character(343:386)),
#     Metric == "Alpha Hill:0",
#     is.na(Subset)
#   ) |> tidytable::left_join(endTimes |> dplyr::select(-Times)) |> tidytable::group_by(
#     SpeciesAffinity, Intervention, PoolPatchSeed
#   ) |> tidytable::filter(
#     Time == max(16300, min(Time)) | Time == 30000
#   ) |> tidytable::group_by(
#     SpeciesAffinity, PoolPatchSeed,
#     Intervention, InterventionInitial, InterventionFinal
#   ) |> tidytable::arrange(
#     Time
#   ) |> tidytable::summarise(
#     Value = Value[2] - Value[1],
#     Method = "Temporal"
#   ),
#   # Counterfactual comparison TODO, ONLY 100% VALID AFTER RE-RUNS WITH FIXED
#   #                                 POOL PREFERENCE ASSIGNMENTS.
#   diversitiesRichness |> tidytable::filter(
#     NicheDistance == defaultNicheDistance,
#     (PoolPatchSeed %in% as.character(343:386)),
#     Metric == "Alpha Hill:0",
#     is.na(Subset)
#   ) |> tidytable::left_join(endTimes |> dplyr::select(-Times)) |> tidytable::group_by(
#     SpeciesAffinity, Intervention, PoolPatchSeed
#   ) |> tidytable::filter(
#     Time == 30000
#   ) |> tidytable::group_by(
#     SpeciesAffinity, PoolPatchSeed,
#     # InterventionFinal # What if you had always been in your final state?
#     InterventionInitial # What if you had stayed in your initial state?
#   ) |> tidytable::mutate(
#     Value = Value - Value[InterventionInitial == InterventionFinal],
#     Method = "Counterfactual"
#   ) |> tidytable::select(
#     SpeciesAffinity, PoolPatchSeed,
#     Intervention, InterventionInitial, InterventionFinal,
#     Value, Method
#   )
#   # ) |> tidytable::filter(
#   #   InterventionFinal == "(0)"
# ) |> tidytable::group_by(
#   SpeciesAffinity:InterventionFinal
# ) |> tidytable::mutate(
#   Deviation = Value - Value[Method=="Counterfactual"]
# ) |> tidytable::group_by(
#   SpeciesAffinity, Intervention:InterventionFinal, Method
# ) |> tidytable::summarise(
#   Mean = mean(Value),
#   StDev = sd(Value),
#   Bias = mean(Deviation),
#   MeanAbsDev = mean(abs(Deviation)),
#   PADGT1 = sum(abs(Deviation) > 1)/tidytable::n(),
#   PADGT3 = sum(abs(Deviation) > 3)/tidytable::n(),
#   PADGT5 = sum(abs(Deviation) > 5)/tidytable::n()
# )
#
# temporalVCounterfactualStats |> pivot_wider(
#   names_from = Method, values_from = StDev,
#   id_cols = c(SpeciesAffinity,
#               Intervention, InterventionInitial, InterventionFinal)
# ) |> filter(
#   Counterfactual != 0
# ) |> mutate(
#   out = (Temporal - Counterfactual)
# ) |> pull(out) |> quantile(probs = (seq(from = 0.05, by = 0.05, to = 0.95)))
#
# temporalVCounterfactualStats |> filter(
#   Method != "Counterfactual"
# ) |> pull(Bias) |> summary()
# temporalVCounterfactualStats |> filter(
#   Method != "Counterfactual"
# ) |> pull(Bias) |> quantile(probs = (seq(from = 0.05, by = 0.05, to = 0.95)))

##### b: ######################################################################
# Could probably add in the counterfactual as well, but might be too messy?
newplot4_b <- ggplot2::ggplot(
  newplot4_dataB,
  aes(x = Time, y = Value,
      group = interaction(SpeciesAffinity, Intervention),
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

colorPalette_0.5 <- colorPalette[
  grepl(x = names(colorPalette), pattern = "^[(]0.5[])]")
  ]
newplot4_b_smooths <- ggplot2::ggplot_build(
  newplot4_b
)$data[[2]] |> dplyr::group_by(
  group
) |> dplyr::filter(
  x == max(x)
) |> dplyr::ungroup(
) |> dplyr::mutate(
  Subset = rev(levels(factor(newplot4_dataB$Subset)))[PANEL],
  Intervention = names(colorPalette_0.5)[
    apply(outer(colour, colorPalette_0.5, `==`), 1, which)
    ],
  yshift = y + c(-1, +2.5, -2, -4, -1, -3, -2, +2.5)
)

newplot4_b <-
  # Witchcraft from stackoverflow.com/a/6675163
  # works by pre-building the plot and then extracting coordinates.
  newplot4_b + ggplot2::coord_cartesian(
    xlim = c(0, 14), clip = "off"
  ) + ggplot2::geom_segment(
    data = newplot4_b_smooths,
    mapping = ggplot2::aes(x = x+1, y = yshift, xend = x, yend = y,
                           color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + ggplot2::geom_label(
    data = newplot4_b_smooths,
    mapping = ggplot2::aes(x = x+2.5, y = yshift,
                           label = Intervention, color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  )

newplot4 <- ggpubr::ggarrange(
  plotlist = list(
    newplot4_a,
    newplot4_b
  ), nrow = 1, common.legend = TRUE #, widths = c(0.5, 0.27, 0.23)
) + ggplot2::annotate(
  "curve", x = 0.29, y = 0.1, xend = 0.5, yend = 0.15,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
)

ggplot2::ggsave(plot = newplot4, filename = "Figure4_Prototype2.png",
                units = "cm", width = 6.5*3, height = 6.5*2)