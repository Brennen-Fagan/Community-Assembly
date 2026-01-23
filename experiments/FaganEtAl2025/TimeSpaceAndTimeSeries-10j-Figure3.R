# Setup: ######################################################################
# Plot of Richness and Abundance split up by basal and consumer species.
# Also functinally an overview plot of network structure.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")
source(file.path("R", "flattenDiversity.R")) # Req'd by below
source(file.path("R", "generateNetworks.R")) # To create inset graphs.

# This is better as an environment, but that's more opaque.
figure3 <- list(
  pref = "100% 0"#"Uniform(0, 1)"
)

# Main Plots: #################################################################
### Plot 3: ###################################################################
##### Data: ###################################################################
figure3$data <- tidytable::bind_rows(
  diversitiesRichness |> tidytable::filter(
    SpeciesPreferences == figure3$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0",
    !is.na(Subset)
  ),
  diversitiesAbund |> tidytable::filter(
    SpeciesPreferences == figure3$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance",
    !is.na(Subset)
  )
) |> tidytable::mutate(
  Subset = factor(Subset, levels = c("Consumer_0", "Basal_0"), 
                  labels = c("Consumer", "Basal"), ordered = TRUE),
  Metric = factor(Metric, levels = c("Alpha Hill:0", "Alpha Abundance"), 
                  labels = c("Richness", "Abundance"), ordered = TRUE)
) |> tidytable::filter(
  Time > Start, Time < Stop
) |> tidytable::group_by(
  PoolPatchSeed, Intervention, SpeciesAffinity, Metric, Subset
) |> tidytable::summarise(
  Value = mean(Value), .groups = "drop"
)

##### Core: ###################################################################
# 2x2 facet plot of violins over LU, Richness - Abundance vs Consumer - Basal.
figure3$plotA <- ggplot2::ggplot(
  figure3$data,
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 1,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
) + ggplot2::geom_line(
  data = ~ summarise(group_by(.x, Intervention, Metric, Subset), 
                     Value = mean(Value), # Avg. over sims.
                     .groups = "drop"), 
  color = "black", group = 1
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::facet_wrap(
  Subset ~ Metric,
  scales = "free"
)

##### B: ######################################################################
# Basal and Consumer Richness co-vary for our scenarios.
# Note separate from Abundance to have two separate grobs to inset.
figure3$plotB <- ggplot2::ggplot(
  figure3$data |> tidytable::filter(
    Metric == "Richness"
    ) |> tidytable::pivot_wider(
    names_from = Subset, values_from = Value
  ),
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = Intervention
  )
) + ggplot2::geom_point(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  x = "Basal Richness",
  y = "Consumer Richness"
) + ggplot2::guides(
  color = "none",
  fill = "none"
  # ) + ggplot2::scale_x_log10(
)

##### C: ######################################################################
# Basal and Consumer Abundance co-vary for our scenarios.
# Note separate from Richness to have two separate grobs to inset.
figure3$plotC <- ggplot2::ggplot(
  figure3$data |> tidytable::filter(
    Metric == "Abundance"
  ) |> tidytable::pivot_wider(
    names_from = Subset, values_from = Value
  ),
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = Intervention
  )
) + ggplot2::geom_point(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  x = "Basal Abundance",
  y = "Consumer Abundance"
) + ggplot2::guides(
  color = "none",
  fill = "none"
  # ) + ggplot2::scale_x_log10(
)

##### Combine: ################################################################
figure3$plot <- figure3$plotA + ggplot2::annotation_custom(
  grob = ggplot2::ggplotGrob(figure3$plotB),
  xmin = , xmax = , ymin = , ymax = 
) + ggplot2::annotation_custom(
  grob = ggplot2::ggplotGrob(figure3$plotC),
  xmin = , xmax = , ymin = , ymax = 
)

ggplot2::ggsave(plot = figure2$plot, filename = file.path(dirImages, "Figure3_Prototype8.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plot, filename = file.path(dirImages, "Figure3_Prototype8.png"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotA, filename = file.path(dirImages, "Figure3A_Prototype8.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotB, filename = file.path(dirImages, "Figure3B_Prototype8.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotC, filename = file.path(dirImages, "Figure3C_Prototype8.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)