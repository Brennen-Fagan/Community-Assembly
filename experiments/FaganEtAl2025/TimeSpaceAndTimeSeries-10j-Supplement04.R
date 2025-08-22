# Setup: ######################################################################

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")


supplement4 <- list()

### 4 Supplement: ############################################################
##### bs1: ####################################################################
newplot4_bs <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  InterventionInitial %in% c("(0)", "(1)"),
  SpeciesPreferences == "100% 0",
  !is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, PoolPatchSeed,
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
      SpeciesPreferences, Intervention, # SpeciesAffinity not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
    # )
  ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesPreferences = ifelse(is.na(SpeciesPreferences), "100% 0", SpeciesPreferences)
) |> tidytable::filter(
  InterventionInitial %in% c("(0)", "(1)")
) |> ggplot2::ggplot(
  aes(x = Time, y = Value,
      group = interaction(SpeciesPreferences, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  # show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  # tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  InterventionInitial +
    factor(Subset, levels = c("Consumer_0", "Basal_0"),
           labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::theme(
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

colorPalette_01 <- colorPalette[
  grepl(x = names(colorPalette), pattern = "^[(](0|1)[])]")
  ]
newplot4_bs_smooths <- ggplot2::ggplot_build(
  newplot4_bs
)$data[[2]] |> dplyr::group_by(
  group
) |> dplyr::filter(
  x == max(x)
) |> dplyr::ungroup(
) |> dplyr::mutate(
  Subset = rev(levels(factor(newplot4_dataB$Subset)))[
    ((as.numeric(PANEL) - 1) %% 2) + 1
    ],
  Intervention = names(colorPalette_01)[
    apply(outer(colour, colorPalette_01, `==`), 1, which)
    ],
  InterventionInitial = substr(Intervention, 1, 3),
  yshift = y + c(+2.5, 0, -5, -10,
                 -15, -10, -5, +2,
                 +2, -5, -10, -15,
                 -6, -2, 0, 0)
)

newplot4_bs <-
  # Witchcraft from stackoverflow.com/a/6675163
  # works by pre-building the plot and then extracting coordinates.
  newplot4_bs + ggplot2::coord_cartesian(
    xlim = c(0, 13), clip = "off"
  ) + ggplot2::geom_segment(
    data = newplot4_bs_smooths,
    mapping = ggplot2::aes(x = x+1, y = yshift, xend = x, yend = y,
                           color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + ggplot2::geom_label(
    data = newplot4_bs_smooths,
    mapping = ggplot2::aes(x = x+2, y = yshift,
                           label = Intervention, color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  )

ggplot2::ggsave(plot = newplot4_bs, filename = "Figure_supplement4_v1.png",
                units = "cm", width = 6.5*3, height = 6.5*2)
