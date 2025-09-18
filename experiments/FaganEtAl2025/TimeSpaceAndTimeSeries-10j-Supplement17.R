# Setup: ######################################################################
# Plot of long-term and short-term responses of species richness to land-use
# change, taking 0.5 as a base case for sake of argument.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

supplement17 <- list()
supplement17$baseCaseVersion <- 5

supplement17$baseCase <-
  rep(c("(0)", "(0.5)", "(1)"), 3)[
    supplement17$baseCaseVersion
    ]

supplement17$basePreference <-
  rep(c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"), each = 3)[
    supplement17$baseCaseVersion
    ]

# Main Plots: #################################################################
### Plot 4:####################################################################
supplement17$dataA <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences == supplement17$basePreference,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  is.na(Subset) # Aggregates
)

supplement17$rangeXMin <- supplement17$dataA |> tidytable::group_by(
  PoolPatchSeed, Intervention
) |> tidytable::filter(
  Time == min(Time) & Time != 0
) |> tidytable::pull(
  Time
) |> range(
)
supplement17$rangeY <- c(
  1e-5, 
  supplement17$dataA |> tidytable::filter(
    InterventionInitial == supplement17$baseCase
  ) |> tidytable::pull(Value) |> max()
)

supplement17$interventionTimes <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences == supplement17$basePreference,
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

supplement17$dataB <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences == supplement17$basePreference,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  is.na(Subset),
  Intervention == supplement17$baseCase
) |> tidytable::left_join(
  supplement17$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ),
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::mutate(
  Time = Time - InterventionTime
) |> tidytable::filter(
  Time > -50
)

supplement17$dataC <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences == supplement17$basePreference,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  !is.na(Subset), # Basals and Consumers
  InterventionInitial == supplement17$baseCase,
  InterventionFinal %in% c("(0)", "(0.5)", "(1)"),
  InterventionInitial != InterventionFinal#,
  # Time %in% c(Time[1:140]) # 30 is 1:10, 12:20:2, 23:50:3, 100 gets to 600
) |> tidytable::separate(
  Subset, into = c("Guild", "AffinityBins"), sep = "_"
  # ) |> unifyAffinityBins(# Strictly unnecessary here.
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Time, Guild
) |> tidytable::summarise( # Across AffinityBins
  Value = sum(Value)
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Guild
) |> tidytable::arrange(
  Time
) |> tidytable::summarise(
  Time = Time - Time[2],
  Value = Value - Value[2],
  Method = "Temporal",
  .groups = "drop"
  # ) |> tidytable::filter(
  #   Time <= 1100
) |> tidytable::mutate(
  Weight = ifelse(Time < 1e-6 & Time > 1e-6, 1e9, 1)
)

supplement17$rangeXMax <- supplement17$dataB |> tidytable::group_by(
  PoolPatchSeed, Intervention
) |> tidytable::filter(
  Time == max(Time) & Time != 0
) |> tidytable::pull(
  Time
) |> range(
)

supplement17$dataD <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences == supplement17$basePreference,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  !is.na(Subset), # Basals and Consumers
  InterventionInitial == supplement17$baseCase,
  InterventionFinal %in% c("(0)", "(0.5)", "(1)"),
  InterventionInitial != InterventionFinal
) |> tidytable::mutate(
  Time = round(Time, 6) # Avoid numerical problems
) |> tidytable::separate(
  Subset, into = c("Guild", "AffinityBins"), sep = "_"
  # ) |> unifyAffinityBins(# Strictly unnecessary here.
) |> tidytable::group_by(
  PoolPatch, SpeciesPreferences, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Time, Guild
) |> tidytable::summarise( # Across AffinityBins
  ValueIntervention = sum(Value)
) |> tidytable::left_join(
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == supplement17$basePreference,
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance",
    !is.na(Subset), # Basals and Consumers
    Intervention == supplement17$baseCase
  )  |> tidytable::mutate(
    Time = round(Time, 6) # Avoid numerical problems
  ) |> tidytable::separate(
    Subset, into = c("Guild", "AffinityBins"), sep = "_"
    # ) |> unifyAffinityBins(# Strictly unnecessary here.
  ) |> tidytable::group_by(
    PoolPatch, SpeciesPreferences, PoolPatchSeed,
    InterventionInitial, Time, Guild
  ) |> tidytable::summarise( # Across AffinityBins
    ValueBase = sum(Value)
  ),
  by = c("PoolPatch", "SpeciesPreferences", "PoolPatchSeed",
         "InterventionInitial", "Time", "Guild")
) |> tidytable::left_join(
  supplement17$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ),
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::mutate(
  Time = round(Time - InterventionTime, 6),
  Value = ValueIntervention - ValueBase, # Avoid numerical problems
  Method = "Counterfactual"
)

##### a: ######################################################################
supplement17$plotA <- plotMeanAndInner(
  rbind(
    supplement17$dataA |> tidytable::filter(
      Intervention %in% c(
        supplement17$baseCase
      )
    ),
    # We want to appear in the legend but not on the plot!
    supplement17$dataA |> tidytable::filter(
      PoolPatchSeed == "1",
      Intervention %in% paste0(
        supplement17$baseCase, "->", c("(0)", "(0.5)", "(1)")
      ),
      abs(Time - 20000) == min(abs(Time - 20000))
    ) |> tidytable::mutate(
      Value = -100
    )
  ),
  CIs = 0.75, facets = as.formula(. ~ .)
  ###### Annotation: ##########################################################
) + ggplot2::labs(
  y = "Total Abundance", x = "Time Since Simulation Start"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 7),
  fill = ggplot2::guide_legend(ncol = 7)
) + ggplot2::coord_cartesian(
  xlim = c(0, supplement17$rangeXMin[1]-1), 
  ylim = supplement17$rangeY,
  expand = FALSE
  # ) + ggplot2::geom_rect(
  #   data = data.frame(
  #     1 # 1 rectangle per row, so dummy df to prevent overplotting
  #   ),
  #   xmin = 16300,
  #   xmax = 21300,
  #   ymin = 0, ymax = richnessYMax,
  #   fill = "grey",
  #   alpha = 0.4,
  #   inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:4)*4000
)

##### b: ######################################################################
supplement17$plotB <- ggplot2::ggplot(
  supplement17$dataA |> tidytable::filter(
    Time >= supplement17$rangeXMin[1]-10, Time <= supplement17$rangeXMin[2]+10,
    Intervention %in% c(
      supplement17$baseCase,
      paste0(
        supplement17$baseCase, "->", c("(0)", "(0.5)", "(1)")
      )
    )
  ),
  ggplot2::aes(
    x = Time, y = Value,
    color = Intervention,
    group = interaction(PoolPatchSeed, Intervention),
    alpha = (PoolPatchSeed == "1")
  )
) + ggplot2::geom_line(
  show.legend = FALSE
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::labs(
  y = "Total Abundance",
  x = "Time Since Simulation Start"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.025, 1),
  axis.text.y = ggplot2::element_blank(),
  axis.title.y = ggplot2::element_blank()
) + ggplot2::coord_cartesian(
  xlim = supplement17$rangeXMin,
  ylim = supplement17$rangeY,
  expand = FALSE
)

##### c: ######################################################################
supplement17$plotC <- plotMeanAndInner(
  rbind(
    supplement17$dataA |> tidytable::filter(
      Intervention %in% paste0(
        supplement17$baseCase, "->", c("(0)", "(0.5)", "(1)")
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
    supplement17$dataB |> tidytable::mutate(
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
  y = "Total Abundance", x = "Time Since Intervention"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 7),
  fill = ggplot2::guide_legend(ncol = 7)
) + ggplot2::coord_cartesian(
  xlim = c(0, supplement17$rangeXMax[1]-1),
  ylim = supplement17$rangeY, 
  expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 5000,
  ymin = 0, ymax = Inf,
  # ymax = richnessYMax,
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

##### Top: #####################################################################

supplement17$plotTOP <- ggpubr::ggarrange(
  plotlist = list(
    supplement17$plotA + ggplot2::scale_y_log10(),
    supplement17$plotB + ggplot2::scale_y_log10(),
    supplement17$plotC + ggplot2::scale_y_log10()
  ),
  nrow = 1, common.legend = TRUE, heights = c(0.3, 0.2, 0.3)
)

##### c: ######################################################################
supplement17$plotD <- ggplot2::ggplot(
  supplement17$dataC,
  aes(x = Time, y = Value,
      group = interaction(SpeciesPreferences, Intervention),
      color = Intervention
  )
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6000,
  ymin = -Inf, ymax = Inf,
  # ymin = -richnessYMax, ymax = richnessYMax,
  fill = "grey",
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::geom_point(
  show.legend = FALSE, alpha = 0.02
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
    SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal,
    Guild, Time, Method
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  show.legend = FALSE
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
    SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal,
    Guild, Time, Method
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  color = "black", linetype = "dotted",
  show.legend = FALSE
  ###### Annotation: ##########################################################
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales After Interventions",
  caption = "Log(1+Time) Axis to show differing regimes.",
  tag = "d)",
  x = "Time Since Intervention", #"Time since Land-use Change",
  y = "Impact - Control (Total Abundance)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  factor(Guild, levels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::theme(
  plot.background = ggplot2::element_rect(linetype = "solid")
) + ggplot2::coord_cartesian(
  xlim = c(0, 5000), expand = FALSE
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 2000, 5000),
  transform = "log1p"
) + ggplot2::scale_y_continuous(
  transform = scales::transform_modulus(0),
  breaks = round(c(
    -rev(10^seq(
      from = 1, to = round(log10(supplement17$rangeY[2]), 1), length.out = 3)),
    0,
    10^seq(
      from = 1, to = round(log10(supplement17$rangeY[2]), 1), length.out = 3)
  ))
)

##### Composite: ##############################################################

supplement17$plot <- ggpubr::ggarrange(
  plotlist = list(
    supplement17$plotTOP,
    supplement17$plotD
  ),
  ncol = 1, common.legend = TRUE, heights = c(0.4, 0.6)
)

ggplot2::ggsave(plot = supplement17$plot,
                filename = paste0(
                  "Figure_supplement17_v2_", supplement17$baseCaseVersion, ".pdf"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = supplement17$plot,
                filename = paste0(
                  "Figure_supplement17_v2_", supplement17$baseCaseVersion, ".png"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)
