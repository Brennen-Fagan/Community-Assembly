##### bs2: ####################################################################
newplot4_bs2 <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  # InterventionInitial %in% c("(0)", "(1)"),
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
  )
  # ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesAffinity =
    ifelse(is.na(SpeciesAffinity), "100% 0", SpeciesAffinity)
) |> tidytable::filter(
  # InterventionInitial %in% c("(0)", "(1)")
  Subset %in% c("Basal_0", "Consumer_0"),
  !is.na(InterventionInitial), !is.na(InterventionFinal),
  InterventionInitial != InterventionFinal
) |> ggplot2::ggplot(
  aes(x = Time, y = Value,
      group = interaction(SpeciesAffinity, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  # show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x",
  show.legend = FALSE
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  subtitle =
    "Columns = Final Land-use, Rows = Initial Land-use and Species Type",
  # tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  InterventionInitial +
    factor(Subset, levels = c("Consumer_0", "Basal_0"),
           labels = c("C0", "B0"),
           ordered = TRUE) ~ InterventionFinal
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  # strip.text.x = ggplot2::element_blank()
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

ggplot2::ggsave(plot = newplot4_bs2, filename = "Figure4s2_Prototype1.png",
                units = "cm", width = 6.5*5, height = 6.5*4)