# Setup: ######################################################################
# Plot of long-term and short-term responses of species richness to land-use
# change, taking 0.5 as a base case for sake of argument.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

figure4 <- list()
figure4$baseCaseVersion <- 1

figure4$baseCase <- switch(
  figure4$baseCaseVersion,
  "(0)", "(0.5)", "(1)"
)

# Main Plots: #################################################################
### Plot 4:####################################################################
figure4$dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset) # Aggregates
)

figure4$rangeXMin <- figure4$dataA |> tidytable::group_by(
  PoolPatchSeed, Intervention
) |> tidytable::filter(
  Time == min(Time) & Time != 0
) |> tidytable::pull(
  Time
) |> range(
) # Range of times over which interventions occur.

figure4$interventionTimes <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  !is.na(Subset), # Basals and Consumers
  InterventionInitial != InterventionFinal
) |> tidytable::select(
  PoolPatch, PoolPatchSeed, Time
) |> tidytable::group_by(
  PoolPatch, PoolPatchSeed
) |> tidytable::summarise(
  Time = min( # Smallest time that's not the smallest time.
    Time[round(Time, 6)!=round(min(Time), 6)]
  ), # Not min(Time) but the next time.
  .groups = "drop"
)

figure4$dataB <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset),
  Intervention == figure4$baseCase
) |> tidytable::left_join(
  figure4$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ), 
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::mutate(
  Time = Time - InterventionTime
) |> tidytable::filter(
  Time > -50
)

figure4$dataC <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == "100% 0",
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  !is.na(Subset), # Basals and Consumers
  InterventionInitial == figure4$baseCase,
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
  # ) |> tidytable::filter(
  #   Time <= 1100
)

figure4$rangeXMax <- figure4$dataB |> tidytable::group_by(
  PoolPatchSeed, Intervention
) |> tidytable::filter(
  Time == max(Time) & Time != 0
) |> tidytable::pull(
  Time
) |> range(
)

##### c: ######################################################################
figure4$plotC <- plotMeanAndInner(
  rbind(
    figure4$dataA |> tidytable::filter(
      Intervention %in% paste0(
        figure4$baseCase, "->", c("(0)", "(0.5)", "(1)")
      )
    ) |> tidytable::group_by(
      SpeciesPreferences, Intervention, PoolPatchSeed,
      InterventionInitial, InterventionFinal, Subset
    ) |> tidytable::arrange(
      Time
    ) |> tidytable::filter(
      InterventionInitial != InterventionFinal#,
      # Time %in% c(Time[1:140]) # 30 is 1:10, 12:20:2, 23:50:3, 100 gets to 600
    ) |> tidytable::mutate(
      Time = round(Time - Time[2], 7),
      Method = "Temporal"
    ),
    figure4$dataB |> tidytable::mutate(
      Time = round(Time, 7),
      Method = "Temporal"
    ) |> tidytable::select(
      -InterventionTime
    )
    # Group so that largest rounding without collisions
  ) |> tidytable::filter(
    # Special rounding early where there are ~2 contemporaneous time series.
    Time > 50 | (abs(Time - round(Time)) < 1e-6)
  ) |> tidytable::group_by(
    SpeciesPreferences, Intervention, PoolPatchSeed,
    InterventionInitial, InterventionFinal, Subset
  ) |> tidytable::mutate( 
    Time = tidytable::case_when(
      Time <= 50 ~ round(Time, 0), 
      Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
      Time < 16350 ~ round(Time, -2),
      TRUE ~ Time
    )
  ) |> tidytable::filter(
    (abs(Time - round(Time)) < 1e-6)
  ) |> tidytable::group_by(
    SpeciesPreferences, Intervention, PoolPatchSeed,
    InterventionInitial, InterventionFinal, Subset,
    Time
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  CIs = 0.75, facets = as.formula(. ~ .), digits = 0
  ###### Annotation: ##########################################################
) + ggplot2::labs(
  y = "Richness", x = "Time Since Intervention"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 7),
  fill = ggplot2::guide_legend(ncol = 7)
) + ggplot2::coord_cartesian(
  xlim = c(0, figure4$rangeXMax[1]-1),
  ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 5000,
  ymin = 0, ymax = richnessYMax,
  fill = "grey",
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "c)"
) + ggplot2::theme(
  plot.tag.position = c(0.025, 1),
  axis.text.y = ggplot2::element_blank(),
  axis.title.y = ggplot2::element_blank()
)