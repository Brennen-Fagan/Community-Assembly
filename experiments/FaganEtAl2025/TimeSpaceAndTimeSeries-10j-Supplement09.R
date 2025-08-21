##### Examine the difference in distributions of sizes visually. ##############
newplot6_Data <- Pers |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  # SpeciesAffinity == "100% 0",
  # SpeciesAffinity == "50% 0, 50% 1",
  # SpeciesAffinity == "Uniform(0, 1)",
  # InterventionInitial == InterventionFinal,
  Stop > In, Start < Out
) |> tidytable::group_by(
  Intervention, InterventionInitial, InterventionFinal, Size, SpeciesAffinity
) |> tidytable::summarise(
  Weight = sum(Out - In),
  .groups = "drop"
) |> dplyr::group_by(
  Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
) |> dplyr::mutate(
  Weight = Weight / sum(Weight)
) |> dplyr::ungroup(
)
newplot_6a <- ggplot2::ggplot(
  newplot6_Data,
  ggplot2::aes(x = Size, y = Weight, color = Intervention)
) + ggplot2::geom_col(
  show.legend = FALSE
  # ) + ggplot2::geom_line(
  #   data = ~ .x |> dplyr::arrange(Size) |> dplyr::group_by(
  #     SpeciesAffinity, Intervention, InterventionInitial, InterventionFinal
  #   ) |> dplyr::mutate(Weight = cumsum(Weight)),
  #   show.legend = FALSE
) + ggplot2::facet_grid(
  SpeciesAffinity + InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_x_log10(
) + ggplot2::geom_vline(
  xintercept = 0.1
)
newplot_6b <- ggplot2::ggplot(
  newplot6_Data,
  ggplot2::aes(x = Size, y = Weight, color = Intervention)
  # ) + ggplot2::geom_col(
  #   show.legend = FALSE
) + ggplot2::geom_line(
  data = ~ .x |> dplyr::arrange(Size) |> dplyr::group_by(
    SpeciesAffinity, Intervention, InterventionInitial, InterventionFinal
  ) |> dplyr::mutate(Weight = cumsum(Weight)),
  show.legend = FALSE
) + ggplot2::facet_grid(
  SpeciesAffinity ~.#+ InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_x_log10(
) + ggplot2::geom_vline(
  xintercept = 0.1
)

newplot6 <- ggarrange(newplot_6a, newplot_6b)
ggplot2::ggsave(
  newplot6,
  filename = "Figure6s_Prototype1.png", # Uniform(0, 1)
  units = "cm", width = 20*3, height = 20*2
)