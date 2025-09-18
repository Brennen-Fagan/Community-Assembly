# Setup: ######################################################################
# Plot of long-term and short-term responses of species richness to land-use
# change, taking 0.5 as a base case for sake of argument.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

supplement18 <- list()
supplement18$interventions <- c(
  "(0)", "(0)->(0.5)", "(0)->(1)",
  "(0.5)->(0)", "(0.5)", "(0.5)->(1)",
  "(1)->(0)", "(1)->(0.5)", "(1)"
)

# Main Plots: #################################################################
### Data:######################################################################
# We need richness, abundance post-intervention, aggregated to guild level, for
# each of our mainland scenarios as they pass through the interventions.

supplement18$interventionTimes <- diversitiesAbund |> tidytable::filter(
  # SpeciesPreferences == supplement17$basePreference,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  !is.na(Subset), # Basals and Consumers
  InterventionInitial != InterventionFinal
) |> tidytable::select(
  PoolPatch, PoolPatchSeed, Time
) |> tidytable::group_by(
  PoolPatch, PoolPatchSeed
) |> tidytable::summarise(
  Time = min(
    Time[round(Time, 6)!=round(min(Time), 6)]
  ), # Not min(Time) but the next time.
  .groups = "drop"
)

supplement18$dataA <- diversitiesAbund |> tidytable::filter(
  # SpeciesPreferences == supplement17$basePreference,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  !is.na(Subset), # Aggregates
  Intervention %in% supplement18$interventions
) |> tidytable::left_join(
  supplement18$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ),
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::mutate(
  Time = Time - InterventionTime
) |> tidytable::filter(
  Time > -50
) |> tidytable::separate(
  Subset, into = c("Guild", "AffinityBins"), sep = "_"
) |> unifyAffinityBins(
)

supplement18$dataB <- diversitiesRichness |> tidytable::filter(
  # SpeciesPreferences == supplement16$basePreference,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  !is.na(Subset), # Aggregates
  Intervention %in% supplement18$interventions
) |> tidytable::left_join(
  supplement18$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ),
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::mutate(
  Time = Time - InterventionTime
) |> tidytable::filter(
  Time > -50
) |> tidytable::separate(
  Subset, into = c("Guild", "AffinityBins"), sep = "_"
) |> unifyAffinityBins(
)

# Plot: #######################################################################
supplement18$plotA <- ggplot2::ggplot(
  supplement18$dataA |> tidytable::filter(
    abs(Time - round(Time)) < 1e-6 | Time >= 55
  ) |> tidytable::mutate(
    Time = tidytable::case_when(
      Time <= 50 ~ round(Time, 0),
      Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
      Time < 16350 ~ round(Time, -2),
      TRUE ~ Time
    )
  ) |> tidytable::group_by(
    # Average Over the now grouped times.
    PoolPatchSeed, SpeciesPreferences, 
    Intervention, InterventionInitial, InterventionFinal,
    Guild, AffinityBins, Time
  ) |> tidytable::summarise(
    Value = mean(Value)
  ) |> tidytable::group_by(
    # Average Over the Affinity Bins.
    PoolPatchSeed, SpeciesPreferences, 
    Intervention, InterventionInitial, InterventionFinal,
    Guild, Time
  ) |> tidytable::summarise(
    Value = sum(Value)
  ),
  aes(x = Time, y = Value,
      linetype = SpeciesPreferences,
      group = interaction(SpeciesPreferences, Intervention),
      color = Intervention
  )
# ) + ggplot2::geom_point(
#   show.legend = FALSE, alpha = 0.02
) + ggplot2::geom_line(
  data = function(x) x |> tidytable::filter(
    abs(Time - round(Time)) < 1e-6 | Time >= 55
  ) |> tidytable::mutate(
    Time = tidytable::case_when(
      Time <= 50 ~ round(Time, 0),
      Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
      Time < 16350 ~ round(Time, -2),
      TRUE ~ Time
    )
  ) |> tidytable::group_by(
    SpeciesPreferences, # Summarise over simulations
    Intervention, InterventionInitial, InterventionFinal,
    Guild, Time
  ) |> tidytable::summarise(
    Value = mean(Value)
  )
) + ggplot2::geom_line(
  data = function(x) x |> tidytable::filter(
    abs(Time - round(Time)) < 1e-6 | Time >= 55
  ) |> tidytable::mutate(
    Time = tidytable::case_when(
      Time <= 50 ~ round(Time, 0),
      Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
      Time < 16350 ~ round(Time, -2),
      TRUE ~ Time
    )
  ) |> tidytable::group_by(
    SpeciesPreferences, # Summarise over simulations
    Intervention, InterventionInitial, InterventionFinal,
    Guild, Time
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  color = "black", linetype = "dotted",
  show.legend = FALSE
  ###### Annotation: ##########################################################
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  caption = "Log(1+Time) Axis to show differing regimes.",
  x = "Time Since Intervention", #"Time since Land-use Change",
  y = "Total Abundance"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  plot.background = ggplot2::element_rect(linetype = "solid")
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  Guild + InterventionInitial ~ SpeciesPreferences,
  scales = "free_y"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
) + ggplot2::scale_y_log10(
) + ggplot2::coord_cartesian(
  xlim = c(0, 16000),
  expand = FALSE # Instability after 16000 because sims start ending.
)

supplement18$plotB <- ggplot2::ggplot(
  supplement18$dataB |> tidytable::filter(
    abs(Time - round(Time)) < 1e-6 | Time >= 55
  ) |> tidytable::mutate(
    Time = tidytable::case_when(
      Time <= 50 ~ round(Time, 0),
      Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
      Time < 16350 ~ round(Time, -2),
      TRUE ~ Time
    )
  ) |> tidytable::group_by(
    # Average Over the now grouped times.
    PoolPatchSeed, SpeciesPreferences, 
    Intervention, InterventionInitial, InterventionFinal,
    Guild, AffinityBins, Time
  ) |> tidytable::summarise(
    Value = mean(Value)
  ) |> tidytable::group_by(
    # Average Over the Affinity Bins.
    PoolPatchSeed, SpeciesPreferences, 
    Intervention, InterventionInitial, InterventionFinal,
    Guild, Time
  ) |> tidytable::summarise(
    Value = sum(Value)
  ),
  aes(x = Time, y = Value,
      linetype = SpeciesPreferences,
      group = interaction(SpeciesPreferences, Intervention),
      color = Intervention
  )
# ) + ggplot2::geom_point(
#   show.legend = FALSE, alpha = 0.02
) + ggplot2::geom_line(
  data = function(x) x |> tidytable::filter(
    abs(Time - round(Time)) < 1e-6 | Time >= 55
  ) |> tidytable::mutate(
    Time = tidytable::case_when(
      Time <= 50 ~ round(Time, 0),
      Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
      Time < 16350 ~ round(Time, -2),
      TRUE ~ Time
    )
  ) |> tidytable::group_by(
    SpeciesPreferences, # Summarise over simulations
    Intervention, InterventionInitial, InterventionFinal,
    Guild, Time
  ) |> tidytable::summarise(
    Value = mean(Value)
  )
) + ggplot2::geom_line(
  data = function(x) x |> tidytable::filter(
    abs(Time - round(Time)) < 1e-6 | Time >= 55
  ) |> tidytable::mutate(
    Time = tidytable::case_when(
      Time <= 50 ~ round(Time, 0),
      Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
      Time < 16350 ~ round(Time, -2),
      TRUE ~ Time
    )
  ) |> tidytable::group_by(
    SpeciesPreferences, # Summarise over simulations
    Intervention, InterventionInitial, InterventionFinal,
    Guild, Time
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  color = "black", linetype = "dotted",
  show.legend = FALSE
  ###### Annotation: ##########################################################
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  caption = "Log(1+Time) Axis to show differing regimes.",
  x = "Time Since Intervention", #"Time since Land-use Change",
  y = "Total Richness"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  plot.background = ggplot2::element_rect(linetype = "solid")
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  Guild + InterventionInitial ~ SpeciesPreferences,
  scales = "free_y"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
# ) + ggplot2::scale_y_log10(
) + ggplot2::coord_cartesian(
  xlim = c(0, 16000),
  expand = FALSE # Instability after 16000 because sims start ending.
)


ggplot2::ggsave(plot = supplement18$plotA,
                filename = paste0(
                  "Figure_supplement18_v2_a.pdf"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = supplement18$plotA,
                filename = paste0(
                  "Figure_supplement18_v2_a.png"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)

ggplot2::ggsave(plot = supplement18$plotB,
                filename = paste0(
                  "Figure_supplement18_v2_b.pdf"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = supplement18$plotB,
                filename = paste0(
                  "Figure_supplement18_v2_b.png"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)
