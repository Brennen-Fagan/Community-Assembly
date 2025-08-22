# Setup: ######################################################################


source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

### 6 Supplement: #############################################################
##### bs4: ####################################################################
newplot4_dataBS4 <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  # InterventionInitial %in% c("(0)", "(1)"),
  SpeciesPreferences == "Uniform(0, 1)",
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
  AffinityBins = gsub(pattern = "(Basal|Consumer)[_]", replacement = "",
                      x = Subset, perl = TRUE),
  SpeciesGuild = gsub(pattern = "(?=_).+", replacement = "",
                      x = Subset, perl = TRUE),
  .groups = "drop"
  #####
) |> unifyAffinityBins(
) |> tidytable::mutate(
  Subset = paste0(SpeciesGuild, "_", AffinityBins)
)

newplot4_bs4 <- newplot4_dataBS4 |> tidytable::right_join(
  tidytable::expand(
    newplot4_dataBS4,
    tidytable::nesting(
      SpeciesPreferences, Intervention, # SpeciesPreferences not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
  )
  # ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesPreferences =
    ifelse(is.na(SpeciesPreferences), "Uniform(0, 1)", SpeciesPreferences)
) |> tidytable::filter(
  # InterventionInitial %in% c("(0)", "(1)")
  Subset %in% c("Basal_(0, 0.2]", "Consumer_(0, 0.2]",
                "Basal_(0.2, 0.4]", "Consumer_(0.2, 0.4]",
                "Basal_(0.4, 0.6]", "Consumer_(0.4, 0.6]",
                "Basal_(0.6, 0.8]", "Consumer_(0.6, 0.8]",
                "Basal_(0.8, 1]", "Consumer_(0.8, 1]"),
  !is.na(InterventionInitial), !is.na(InterventionFinal),
  InterventionInitial != InterventionFinal,
  InterventionInitial == "(0.5)"
) |> ggplot2::ggplot(
  aes(x = Time, y = Value,
      group = interaction(SpeciesPreferences, Intervention),
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
    factor(Subset#,
           # levels = c("Consumer_0", "Basal_0",
           #                    "Consumer_1", "Basal_1"),
           # labels = c("C0", "B0", "C1", "B1"),
           # ordered = TRUE
    ) ~ InterventionFinal
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  # strip.text.x = ggplot2::element_blank()
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

ggplot2::ggsave(plot = newplot4_bs4, filename = "Figure4s4_Prototype2.png",
                units = "cm", width = 6.5*5, height = 6.5*4)
