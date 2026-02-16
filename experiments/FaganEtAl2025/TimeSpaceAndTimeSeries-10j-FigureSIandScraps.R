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