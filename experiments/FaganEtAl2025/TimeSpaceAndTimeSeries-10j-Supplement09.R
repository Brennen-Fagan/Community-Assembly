# Setup: ######################################################################


source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")

supplement9 <- list()

### 9 Supplement: #############################################################
##### Examine the difference in distributions of sizes visually. ##############
supplement9$data <- Pers |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Stop > In, Start < Out
) |> tidytable::group_by(
  Intervention, InterventionInitial, InterventionFinal, Size, SpeciesPreferences
) |> tidytable::summarise(
  Weight = sum(Out - In),
  .groups = "drop"
) |> dplyr::group_by(
  Intervention, InterventionInitial, InterventionFinal, SpeciesPreferences
) |> dplyr::mutate(
  Weight = Weight / sum(Weight)
) |> dplyr::ungroup(
)
supplement9$plotA <- ggplot2::ggplot(
  supplement9$data,
  ggplot2::aes(x = Size, y = Weight, color = Intervention)
) + ggplot2::geom_col(
  show.legend = FALSE
) + ggplot2::facet_grid(
  SpeciesPreferences + InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_x_log10(
) + ggplot2::geom_vline(
  xintercept = 0.1
)
supplement9$plotB <- ggplot2::ggplot(
  supplement9$data,
  ggplot2::aes(x = Size, y = Weight, color = Intervention)
) + ggplot2::geom_line(
  data = ~ .x |> dplyr::arrange(Size) |> dplyr::group_by(
    SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal
  ) |> dplyr::mutate(Weight = cumsum(Weight)),
  show.legend = FALSE
) + ggplot2::facet_grid(
  SpeciesPreferences ~.#+ InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_x_log10(
) + ggplot2::geom_vline(
  xintercept = 0.1
)

supplement9$plot <- ggarrange(supplement9$plotA, supplement9$plotB)
ggplot2::ggsave(
  supplement9$plot,
  filename = "Figure_supplement9_v1.png",
  units = "cm", width = 20*3, height = 20*2
)
