# Setup: ######################################################################
# Figure 2 but for basal and consumer richness instead.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

supplement2 <- list()
supplement2$preferences <- 2
supplement2$initial <- 1

### 2 Supplement: #############################################################
supplement2$dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == switch(supplement2$preferences,
                               "100% 0",
                               "50% 0, 50% 1",
                               "Uniform(0, 1)"
  ),
  NicheDistance == defaultNicheDistance,
  InterventionInitial %in% switch(
    supplement2$initial,
    c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"), # But only initial == final
    c("(0)"),
    c("(0.25)"),
    c("(0.5)"),
    c("(0.75)"),
    c("(1)")
  ),
  if(supplement2$initial == 1) {
    InterventionInitial == InterventionFinal
  } else {
    InterventionFinal %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)")
  },
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  !is.na(Subset)
) |> tidytable::separate_wider_delim(
  Subset, delim = "_", names = c("Guild", "AffinityBins"), cols_remove = FALSE
) |> unifyAffinityBins(
) |> tidytable::mutate(
  Subset = Guild # this interacts with plotMeanAndInner.
)

supplement2$dataB <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences == switch(supplement2$preferences,
                               "100% 0",
                               "50% 0, 50% 1",
                               "Uniform(0, 1)"
  ),
  NicheDistance == defaultNicheDistance,
  InterventionInitial %in% switch(
    supplement2$initial,
    c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"), # But only initial == final
    c("(0)"),
    c("(0.25)"),
    c("(0.5)"),
    c("(0.75)"),
    c("(1)")
  ),
  if(supplement2$initial == 1) {
    InterventionInitial == InterventionFinal
  } else {
    InterventionFinal %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)")
  },
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance",
  !is.na(Subset)
) |> tidytable::separate_wider_delim(
  Subset, delim = "_", names = c("Guild", "AffinityBins"), cols_remove = FALSE
) |> unifyAffinityBins(
) |> tidytable::mutate(
  Subset = Guild # this interacts with plotMeanAndInner.
)

supplement2$dataC <- with(supplement2, {
  # Average Abundance of species in the patch
  rbind(dataA, dataB) |> tidytable::filter(
    Start <= Time, Time <= Stop
  ) |> tidytable::pivot_wider(
    names_from = Metric, values_from = Value
    # ) |> tidytable::mutate(
    #   `Average Abundance` =
    #     ifelse(`Alpha Hill:0` > 0, `Alpha Abundance` / `Alpha Hill:0`, 0)
    # ) |> tidytable::group_by(
    #   # Average over relevant times
    #   Environment1:DispersalParam, Guild
    # ) |> tidytable::summarise(
    #   `Average Abundance` = mean(`Average Abundance`)
  )
}) # 44 replicates, 5 treatments, 2 guilds, 1 2 or 5 affinity bins.

##### a: ######################################################################
# Richness through time for each guild.
supplement2$plotA <- plotMeanAndInner(
  # rbind(
  supplement2$dataA |> tidytable::group_by(
    PoolPatchSeed, Intervention, InterventionInitial, InterventionFinal,
    SpeciesPreferences, Subset, Guild, Time
  ) |> tidytable::summarise(
    Value = sum(Value) # Overall richness across guilds at a time
  ) ,
  # |> tidytable::filter(
  # Intervention %in% c("(0)", "(0.5)", "(1)")
  # ),
  # We want to appear in the legend but not on the plot!
  # supplement2$dataA |> tidytable::filter(
  #   PoolPatchSeed == "1",
  #   Intervention %in% c("(0.25)", "(0.75)"),
  #   abs(Time - 0) == min(abs(Time - 0))
  # ) |> tidytable::mutate(
  #   Value = -100
  # )
  # ),
  CIs = 0.75, facets = as.formula(
    factor(Subset, levels = c("Consumer", "Basal"),
           labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none", color = "none", fill = "none"
  # color = ggplot2::guide_legend(ncol = 5),
  # fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 31000),# ylim = c(0, richnessYMax),
  expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(supplement2$dataA$Start),
  xmax = max(supplement2$dataA$Stop),
  ymin = 0, ymax = richnessYMax,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  # legend.position = c(0.5, 0.9),
  legend.position = "none",
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
# Total richness averaged through time for each guild.
supplement2$plotB <- ggplot2::ggplot(
  supplement2$dataA |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, InterventionFinal,
    SpeciesPreferences, Subset, Guild, Time
  ) |> tidytable::summarise(
    Value = sum(Value) # Overall richness across guilds at a time
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, InterventionFinal,
    SpeciesPreferences, Subset, Guild
  ) |> tidytable::summarise(
    Value = mean(Value) # Average richness over time.
  ),
  ggplot2::aes(
    x = InterventionFinal,
    y = Value,
    color = Intervention
  )
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = 0, ymax = richnessYMax,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::geom_line(
  data =
    ~ summarise(group_by(.x, Intervention, InterventionFinal, Subset, Guild),
                Value = mean(Value),
                .groups = "drop"),
  color = "black", group = 1
) + ggplot2::geom_line(
  # Summarise as above, but for each individual species preference
  # (so we don't need to total).
  data = supplement2$dataA |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    # Average across times
    PoolPatchSeed, Intervention, InterventionFinal,
    SpeciesPreferences, Subset, Guild, AffinityBins
  ) |> tidytable::summarise(
    Value = mean(Value), .groups = "drop"
  ) |> tidytable::group_by(
    Intervention, InterventionFinal, Subset, Guild, AffinityBins
  ) |> tidytable::summarise(
    # Average across simulations.
    Value = mean(Value), .groups = "drop"
  ),
  color = "grey50",
  mapping = ggplot2::aes(
    x = InterventionFinal,
    y = Value,
    linetype = AffinityBins,
    group = interaction(AffinityBins, Subset, Guild)
  )
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
) + ggplot2::labs(
  y = "(Avg. over Time) Richness",
  x = "Habitat's Land-use"
# ) + ggplot2::guides(
#   color = "none",
#   fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, NA),
  expand = FALSE
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer", "Basal"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
)

##### c: ######################################################################
# Total abundance averaged through time for each guild.
supplement2$plotC <- ggplot2::ggplot(
  supplement2$dataB |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, InterventionFinal,
    SpeciesPreferences, Subset, Guild, Time
  ) |> tidytable::summarise(
    Value = sum(Value) # Overall value across guilds at a time
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, InterventionFinal,
    SpeciesPreferences, Subset, Guild
  ) |> tidytable::summarise(
    Value = mean(Value) # Average values over time.
    # Present across simulations
  ),
  ggplot2::aes(
    x = InterventionFinal,
    y = Value,
    color = Intervention
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28, alpha = 0.5
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
) + ggplot2::geom_line(
  data =
    ~ summarise(group_by(.x, Intervention, InterventionFinal, Subset, Guild),
                Value = mean(Value),
                .groups = "drop"),
  color = "black", group = 1
) + ggplot2::geom_line(
  # Summarise as above, but for each individual species preference
  # (so we don't need to total).
  data = supplement2$dataB |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    # Average across times
    PoolPatchSeed, Intervention, InterventionFinal,
    SpeciesPreferences, Subset, Guild, AffinityBins
  ) |> tidytable::summarise(
    Value = mean(Value), .groups = "drop"
  ) |> tidytable::group_by(
    Intervention, InterventionFinal, Subset, Guild, AffinityBins
  ) |> tidytable::summarise(
    # Average across simulations.
    Value = mean(Value), .groups = "drop"
  ),
  color = "grey50",
  mapping = ggplot2::aes(
    x = InterventionFinal,
    y = Value,
    linetype = AffinityBins,
    group = interaction(AffinityBins, Subset, Guild)
  )
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "c)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1),
  panel.background = ggplot2::element_rect(
    fill = "grey95", color = NA
  )
) + ggplot2::labs(
  y = "(Avg. over Time) Total Abundance (Density)",
  x = "Habitat's Land-use"
# ) + ggplot2::guides(
#   color = "none",
#   fill = "none"
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer", "Basal"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .,
  scales = "free_y"
) + ggplot2::scale_y_log10(
)

##### d: ######################################################################
# Average abundance averaged through time for each guild.
supplement2$plotD <- ggplot2::ggplot(
  supplement2$dataC |> tidytable::group_by(
    # For each time, collect total abundance and richness across affinity bins.
    Environment1:DispersalParam, Guild, Time
  ) |> tidytable::summarise(
    `Alpha Abundance` = sum(`Alpha Abundance`),
    `Alpha Hill:0` = sum(`Alpha Hill:0`)
  ) |> tidytable::mutate(
    # Then calculate average abundance at that time
    `Average Abundance` =
      ifelse(`Alpha Hill:0` > 0, `Alpha Abundance` / `Alpha Hill:0`, 0)
  ) |> tidytable::group_by(
    # Average over relevant times
    Environment1:DispersalParam, Guild
  ) |> tidytable::summarise(
    `Average Abundance` = mean(`Average Abundance`)
  ),
  ggplot2::aes(
    x = InterventionFinal,
    y = `Average Abundance`,
    color = Intervention
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28, alpha = 0.5
) + ggplot2::geom_line(
  data =
    ~ summarise(group_by(.x, Intervention, InterventionFinal, Subset, Guild),
                `Average Abundance` = mean(`Average Abundance`),
                .groups = "drop"),
  color = "black", group = 1
) + ggplot2::geom_line(
  # Summarise as above, but for each individual species preference
  # (so we don't need to total).
  data = supplement2$dataC |> tidytable::group_by(
    # For each time, collect total abundance and richness across affinity bins.
    Environment1:DispersalParam, Guild, AffinityBins, Time
  ) |> tidytable::summarise(
    `Alpha Abundance` = sum(`Alpha Abundance`),
    `Alpha Hill:0` = sum(`Alpha Hill:0`)
  ) |> tidytable::mutate(
    # Then calculate average abundance at that time
    `Average Abundance` =
      ifelse(`Alpha Hill:0` > 0, `Alpha Abundance` / `Alpha Hill:0`, 0)
  ) |> tidytable::group_by(
    # Average over relevant times
    Environment1:DispersalParam, Guild, AffinityBins
  ) |> tidytable::summarise(
    `Average Abundance` = mean(`Average Abundance`)
  ) |> tidytable::group_by(
    # Average across simulations.
    Intervention, InterventionFinal, Subset, Guild, AffinityBins
  ) |> tidytable::summarise(
    `Average Abundance` = mean(`Average Abundance`), .groups = "drop"
  ),
  color = "grey50",
  mapping = ggplot2::aes(
    x = InterventionFinal,
    y = `Average Abundance`,
    linetype = AffinityBins,
    group = interaction(AffinityBins, Subset, Guild)
  )
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "d)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1),
  panel.background = ggplot2::element_rect(
    fill = "grey95", color = NA
  )
) + ggplot2::labs(
  y = "(Avg. over Time and Species) Abundance (Density)",
  x = "Habitat's Land-use"
# ) + ggplot2::guides(
#   color = "none",
#   fill = "none"
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer", "Basal"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .,
  scales = "free_y"
) + ggplot2::scale_y_log10(
)

supplement2$plot <- ggpubr::ggarrange(
  plotlist = list(
    supplement2$plotA,
    supplement2$plotB,
    supplement2$plotC,
    supplement2$plotD
  ), nrow = 1, widths = c(0.3, 0.4, 0.4, 0.4), common.legend = TRUE
) |> ggpubr::annotate_figure(
  top = ggpubr::text_grob(
    paste0(unique(supplement2$dataA$SpeciesPreferences),
           " Scenarios: Richness and Abundance")
  )
)

ggplot2::ggsave(plot = supplement2$plot,
                filename = paste0(
                  "Figure_supplement2_v4_",
                  supplement2$preferences, "-", supplement2$initial,
                  ".pdf"),
                units = "cm", width = 6.5*5, height = 6.5*2)
ggplot2::ggsave(plot = supplement2$plot,
                filename = paste0(
                  "Figure_supplement2_v4_",
                  supplement2$preferences, "-", supplement2$initial,
                  ".png"),
                units = "cm", width = 6.5*5, height = 6.5*2)
