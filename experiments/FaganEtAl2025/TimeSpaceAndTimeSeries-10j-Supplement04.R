# Setup: ######################################################################
# Fig 4b but for from (0) and from (1).

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

supplement4 <- list()

### 4 Supplement: ############################################################
##### bs1: ####################################################################
supplement4$dataA <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  InterventionInitial %in% c("(0)", "(1)"),
  SpeciesPreferences == "100% 0",
  !is.na(Subset)
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time %in% c(Time[1:30])
) |> tidytable::summarise(
  Time = Time - Time[1],
  Value = Value - Value[1],
  Method = "Temporal",
  .groups = "drop"
) |> tidytable::mutate(
  Weight = ifelse(Time < 1e-6, 1e9, 1), # loess in geom_smooth to anchor to 0.
  Alpha = ifelse(Time <= 10, 0.1, 0)
)

supplement4$plotA <- ggplot2::ggplot(
  supplement4$dataA,
  aes(x = Time, y = Value,
      group = interaction(SpeciesPreferences, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE,
  mapping = ggplot2::aes(alpha = Alpha)
) + ggplot2::geom_smooth(
  show.legend = FALSE,
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

supplement4$colorPalette <- colorPalette[
  grepl(x = names(colorPalette), pattern = "^[(](0|1)[])]")
  ]
supplement4$smoothsA <- ggplot2::ggplot_build(
  supplement4$plotA
)$data[[2]] |> dplyr::group_by(
  group
) |> dplyr::filter(
  x == x[min(abs(x - 10)) == abs(x - 10)]
) |> dplyr::ungroup(
) |> dplyr::mutate(
  Subset = rev(levels(factor(supplement4$dataA$Subset)))[
    ((as.numeric(PANEL) - 1) %% 2) + 1
    ],
  Intervention = names(supplement4$colorPalette)[
    apply(outer(colour, supplement4$colorPalette, `==`), 1, which)
    ],
  InterventionInitial = substr(Intervention, 1, 3),
  yshift = y + c(+3, -2, -8, -15, # Consumer from 0
                 -24, -15, -7, +2, # Basal from 0
                 +3, -6, -16, -25, # Consumer from 1
                 -12, -5, 0, +3) # Basal from 1
)

supplement4$plotA <-
  # Witchcraft from stackoverflow.com/a/6675163
  # works by pre-building the plot and then extracting coordinates.
  supplement4$plotA + ggplot2::coord_cartesian(
    xlim = c(0, 14)#, clip = "off"
  ) + ggplot2::geom_segment(
    data = supplement4$smoothsA,
    mapping = ggplot2::aes(x = x+1, y = yshift, xend = x, yend = y,
                           color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + ggplot2::geom_label(
    data = supplement4$smoothsA,
    mapping = ggplot2::aes(x = x+2, y = yshift,
                           label = Intervention, color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  )

ggplot2::ggsave(plot = supplement4$plotA,
                filename = "Figure_supplement4_v1.png",
                units = "cm", width = 6.5*3, height = 6.5*2)
