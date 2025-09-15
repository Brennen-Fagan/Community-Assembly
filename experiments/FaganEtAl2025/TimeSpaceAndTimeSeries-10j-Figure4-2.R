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
)

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
  Time = min(
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
) |> tidytable::mutate(
  Weight = ifelse(Time < 1e-6, 1e9, 1), # loess in geom_smooth to anchor to 0.
  Alpha = ifelse(Time <= 10, 0.1, 0)
)

figure4$rangeXMax <- figure4$dataB |> tidytable::group_by(
  PoolPatchSeed, Intervention
) |> tidytable::filter(
  Time == max(Time) & Time != 0
) |> tidytable::pull(
  Time
) |> range(
)

##### a: ######################################################################
figure4$plotA <- plotMeanAndInner(
  rbind(
    figure4$dataA |> tidytable::filter(
      Intervention %in% c(
        figure4$baseCase
      )
    ),
    # We want to appear in the legend but not on the plot!
    figure4$dataA |> tidytable::filter(
      PoolPatchSeed == "1",
      Intervention %in% paste0(
        figure4$baseCase, "->", c("(0)", "(0.5)", "(1)")
      ),
      abs(Time - 20000) == min(abs(Time - 20000))
    ) |> tidytable::mutate(
      Value = -100
    )
  ),
  CIs = 0.75, facets = as.formula(. ~ .)
  ###### Annotation: ##########################################################
) + ggplot2::labs(
  y = "Richness", x = "Time Since Simulation Start"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 7),
  fill = ggplot2::guide_legend(ncol = 7)
) + ggplot2::coord_cartesian(
  xlim = c(0, figure4$rangeXMin[1]-1), ylim = c(0, richnessYMax), expand = FALSE
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
figure4$plotB <- ggplot2::ggplot(
  figure4$dataA |> tidytable::filter(
    Time >= figure4$rangeXMin[1]-10, Time <= figure4$rangeXMin[2]+10,
    Intervention %in% c(
      figure4$baseCase,
      paste0(
        figure4$baseCase, "->", c("(0)", "(0.5)", "(1)")
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
  y = "Richness", 
  x = "Time Since Simulation Start"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.025, 1),
  axis.text.y = ggplot2::element_blank(),
  axis.title.y = ggplot2::element_blank()
) + ggplot2::coord_cartesian(
  xlim = figure4$rangeXMin,
  ylim = c(0, richnessYMax), expand = FALSE
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

##### Top: #####################################################################

figure4$plotTOP <- ggpubr::ggarrange(
  plotlist = list(
    figure4$plotA,
    figure4$plotB,
    figure4$plotC
  ),
  nrow = 1, common.legend = TRUE, heights = c(0.3, 0.2, 0.3)
)

##### c: ######################################################################
figure4$plotD <- ggplot2::ggplot(
  figure4$dataC,
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
  ymin = -richnessYMax, ymax = richnessYMax,
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
    Subset, Time, Method, Alpha
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  show.legend = FALSE
  ###### Annotation: ##########################################################
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales After Interventions", 
  caption = "Log(1+Time) Axis to show differing regimes.",
  tag = "d)",
  x = "Time Since Intervention", #"Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::theme(
  plot.background = ggplot2::element_rect(linetype = "solid")
) + ggplot2::coord_cartesian(
  xlim = c(0, 5000), expand = FALSE
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 2000, 5000),
  transform = "log1p"
)

##### Composite: ##############################################################

figure4$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure4$plotTOP,
    figure4$plotD
  ),
  ncol = 1, common.legend = TRUE, heights = c(0.4, 0.6)
)

ggplot2::ggsave(plot = figure4$plot, 
                filename = paste0(
                  "Figure4_Prototype4_", figure4$baseCaseVersion, ".pdf"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure4$plot, 
                filename = paste0(
                  "Figure4_Prototype4_", figure4$baseCaseVersion, ".png"
                ),
                units = "cm", width = 6.5*3, height = 6.5*2)
