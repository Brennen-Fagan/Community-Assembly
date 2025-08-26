# Setup: ######################################################################
# Overview of all interventions for species preferences 100% 0.
# Looking specifically at short term responses.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

supplement5 <- list()

### 5 Supplement: #############################################################
##### bs2: ####################################################################
supplement5$dataA <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  SpeciesPreferences == "100% 0",
  !is.na(Subset)
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time %in% c(Time[1:50])
) |> tidytable::summarise(
  Time = Time - Time[1],
  Value = Value - Value[1],
  Method = "Temporal",
  .groups = "drop"
) |> tidytable::mutate(
  Weight = ifelse(Time < 1e-6, 1e2, 1), # loess in geom_smooth to anchor to 0.
  Alpha = ifelse(Time <= 10, 0.1, 0)
)

supplement5$plotA <- ggplot2::ggplot(
  supplement5$dataA,
  aes(x = Time, y = Value,
      group = interaction(SpeciesPreferences, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_hline(
  yintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_point(
  show.legend = FALSE,
  mapping = ggplot2::aes(alpha = Alpha)
) + ggplot2::geom_smooth(
  # show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  # method = "loess",
  # formula = "y~x",
  show.legend = FALSE,
  color = "black"
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
) + ggplot2::scale_alpha(
  range = c(0, 0.1)
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
) + ggplot2::coord_cartesian(
  xlim = c(0, 10)
)

ggplot2::ggsave(plot = supplement5$plotA,
                filename = "Figure_supplement5_v1.png",
                units = "cm", width = 6.5*5, height = 6.5*4)
