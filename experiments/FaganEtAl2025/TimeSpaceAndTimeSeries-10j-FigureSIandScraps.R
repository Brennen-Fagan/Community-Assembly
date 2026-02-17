### Figure 2: #################################################################
##### e: ######################################################################
# Richness and abundance co-vary for our scenarios.
figure2$plotE <- ggplot2::ggplot(
  figure2$data |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by( # Reduce to per run (x44 sims for param combns)
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) |> tidytable::summarise(
    `Alpha Hill:0` = mean(`Alpha Hill:0`),
    `Alpha Abundance` = mean(`Alpha Abundance`)
  ),
  ggplot2::aes(
    x = `Alpha Abundance`,
    y = `Alpha Hill:0`,
    fill = Intervention
  )
) + ggplot2::geom_point(
  shape = 21, color = "white" # circles (21), squares (22), triangles (24)
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  x = "Abundance",
  y = "Richness"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax),
  xlim = figure2$abundLimits
) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank()
)

if (figure2$abundlog) {
  figure2$plotE <- figure2$plotE + ggplot2::scale_x_log10()
}

##### h: ######################################################################
# Basal and Consumer Richness And Abundance co-vary for our scenarios.
# Note separate from Abundance to have two separate grobs to inset.
figure2$plotH <- ggplot2::ggplot(
  figure2$dataBC |> tidytable::pivot_wider(
    names_from = Metric, values_from = Value
  ),
  ggplot2::aes(
    x = Abundance,
    y = Richness,
    fill = Intervention,
    shape = Subset
  )
) + ggplot2::geom_point(
  color = "white"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  color = "none",
  fill = "none",
  shape = "none"
) + ggplot2::scale_shape_manual(
  values = c(22, 24) # circles (21), squares (22), triangles (24)
) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank()
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax),
  xlim = figure2$abundLimits
)

if (figure2$abundlog) {
  figure2$plotH <- figure2$plotH + ggplot2::scale_x_log10()
}


##### i: ######################################################################
# Basal and Consumer Richness And Abundance co-vary for our scenarios.
# Note separate from Abundance to have two separate grobs to inset.
figure2$plotI <- ggplot2::ggplot(
  figure2$dataBC |> tidytable::filter(
    Subset == "Consumer"
  ) |> tidytable::pivot_wider(
    names_from = Metric, values_from = Value
  ),
  ggplot2::aes(
    x = Abundance,
    y = Richness,
    fill = Intervention,
    shape = Subset
  )
) + ggplot2::geom_point(
  color = "white"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  color = "none",
  fill = "none",
  shape = "none"
) + ggplot2::scale_x_log10(
) + ggplot2::scale_shape_manual(
  values = c(22, 24) # circles (21), squares (22), triangles (24)
) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank()
)

### Figure 4: ################################################################
##### Data: ##################################################################


# Ratios need to be handled slightly differently due to consumer/basal
# resulting in row changes.
figure4$dataBCSummary <- figure4$dataBase |> tidytable::filter(
  Time > -1000,
  Metric %in% c("Richness", "Abundance"),
  !is.na(Subset) # Not overall values
) |> tidytable::separate_wider_delim(
  delim = "_", cols = Subset, names = c("Guild", "AffinityBins")
) |> unifyAffinityBins( # if many preference types.
) |> tidytable::group_by(
  # Aggregate Over the AffinityBins.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Guild, Time
) |> tidytable::summarise(
  Value = sum(Value)
) |> tidytable::pivot_wider(
  names_from = Guild, values_from = Value
) |> tidytable::mutate(
  Time = tidytable::case_when( # Create groupings for times.
    Time < -50 ~ round(Time, -2),
    Time < 0 ~ -25, # In the last bin before regime change.
    Time <= 50 ~ round(Time, 0),
    Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
    Time < 16350 ~ round(Time, -2),
    TRUE ~ Time
  ),
  Subset = NA,
  Value = Consumer/Basal
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, Time
) |> tidytable::summarise(
  Lower = quantile(Value, p = (1 - figure4$CI) - (1 - figure4$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure4$CI + (1 - figure4$CI)/2, na.rm = TRUE)
)

# As in dataBCSummary, but broken up by AffinityBins
figure4$dataBCSupplement2 <- figure4$dataBase |> tidytable::filter(
  Time > -1000,
  Metric %in% c("Richness", "Abundance"),
  !is.na(Subset) # Not overall values
) |> tidytable::separate_wider_delim(
  delim = "_", cols = Subset, names = c("Guild", "AffinityBins")
) |> unifyAffinityBins( # if many preference types.
) |> tidytable::pivot_wider(
  names_from = Guild, values_from = Value
) |> tidytable::mutate(
  Time = tidytable::case_when( # Create groupings for times.
    Time < -50 ~ round(Time, -2),
    Time < 0 ~ -25, # In the last bin before regime change.
    Time <= 50 ~ round(Time, 0),
    Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
    Time < 16350 ~ round(Time, -2),
    TRUE ~ Time
  ),
  Subset = NA,
  Value = Consumer/Basal
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, AffinityBins, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, AffinityBins, Time
) |> tidytable::summarise(
  Lower = quantile(Value, p = (1 - figure4$CI) - (1 - figure4$CI)/2, na.rm = TRUE),
  Average = mean(Value),
  Upper = quantile(Value, p = figure4$CI + (1 - figure4$CI)/2, na.rm = TRUE)
)

##### a: ######################################################################
# Summarised time-series plot for overall richness.
figure4$plotA <- ggplot2::ggplot(
  figure4$dataOverallSummary |> tidytable::filter(
    Metric == "Richness"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax),
  expand = FALSE
)



##### b: ######################################################################
# Summarised time-series plot for short time scale richness ratio.
# Note this should be Consumer / Basal (because Consumers often -> 0).

figure4$plotB <- ggplot2::ggplot(
  figure4$dataBCSummary |> tidytable::filter(
    Metric == "Richness",
    Time > -200, Time < 100
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::coord_cartesian(
  xlim = c(0, 50), expand = FALSE
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness Ratio"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
)

##### c: ######################################################################
# Summarised time-series plot for short time scale abundance ratio.
# Note this should be Consumer / Basal (because Consumers often -> 0).

figure4$plotC <- ggplot2::ggplot(
  figure4$dataBCSummary |> tidytable::filter(
    Metric == "Abundance",
    Time > -200, Time < 100
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::coord_cartesian(
  xlim = c(0, 50), expand = FALSE
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Abundance Ratio"
) + ggplot2::scale_y_log10(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
)

##### SUPPLEMENT: #############################################################
# Short-long term transition for richness and abundance on log(1+Time) scale
# expecting to hit the ratios as well as the base values.

figure4$SupplementRichness <- ggplot2::ggplot(
  figure4$dataBCSupplement |> tidytable::filter(
    Metric == "Richness"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  factor(Guild, levels = c("Consumer", "Basal"), ordered = TRUE) ~
    AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
)

figure4$SupplementAbundance <- ggplot2::ggplot(
  figure4$dataBCSupplement |> tidytable::filter(
    Metric == "Abundance"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Abundance"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  factor(Guild, levels = c("Consumer", "Basal"), ordered = TRUE) ~
    AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
) + ggplot2::scale_y_log10(
)

figure4$SupplementRichnessRatio <- ggplot2::ggplot(
  figure4$dataBCSupplement2 |> tidytable::filter(
    Metric == "Richness"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Richness Ratio"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  . ~
    AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
)

figure4$SupplementAbundanceRatio <- ggplot2::ggplot(
  figure4$dataBCSupplement2 |> tidytable::filter(
    Metric == "Abundance"
  ),
  aes(x = Time, y = Average,
      color = Intervention,
      fill = Intervention
  )
) + ggplot2::geom_vline(
  xintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_hline(
  yintercept = 0, color = "black", linetype = "dashed"
) + ggplot2::geom_line(
) + ggplot2::geom_ribbon(
  ggplot2::aes(ymin = Lower, ymax = Upper),
  alpha = 0.25, linewidth = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Abundance Ratio"
) + ggplot2::coord_cartesian(
  expand = FALSE
) + ggplot2::facet_grid(
  . ~ AffinityBins, scales = "free"
) + ggplot2::scale_x_continuous(
  breaks = c(0, 1, 10, 100, 1000, 10000),
  transform = "log1p"
) + ggplot2::scale_y_log10(
)
