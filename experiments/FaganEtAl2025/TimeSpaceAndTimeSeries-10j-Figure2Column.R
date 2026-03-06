# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use,
# when species preferences are 100% 0.
# Also functinally an overview plot of network structure.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

library(cowplot) # plot arrangement with pagination (ggarrange).
library(tibble) # data.frames that can incorporate grobs.
library(ggpp) # for the geom_grob function -> factor placement of grobs.

# This is better as an environment, but that's more opaque.
figure2 <- list(
  graph = list(
    seed = "2", # "11", "17", "2"!,
    time = 25000
  ),
  abundlog = TRUE,
  pref = "100% 0"#c("Uniform(0, 1)", "50% 0, 50% 1")
)

figure2$graph$specification <- diversitiesRichness |> tidytable::select(c(
  # Which network:
  "Time", "Environment1",
  # Which File (Base):
  "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
  "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
  "Dispersal", "NicheDistance",
  # Oops, there was a collision causing human readable to replace machine.
  # Will be replaced SpeciesAffinity#2 will -> SpeciesPreferences.
  "SpeciesAffinity",
  "SpeciesAffinitySeed", "PatchAffinity", "PatchAffinitySeed",
  # Which File (Intervention):
  "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
  "InterventionTimeSeed", "InterventionDispersal", "InterventionNicheDistance",
  # Ease of Use
  "SpeciesPreferences", "Intervention"
)) |> tidytable::filter(
  SpeciesPreferences == figure2$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  PoolPatchSeed %in% figure2$graph$seed,
  Time == figure2$graph$time
) |> tidytable::distinct(
)

figure2$graph$networks <- generateNetworks(figure2$graph$specification,
                                           Date = "2025-07-30", split = TRUE)
# Main Plots: #################################################################
### Plot 2: ###################################################################
##### Data: ###################################################################
# Richness data: should be straightforward.
figure2$dataRich <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences %in% figure2$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0"
)
figure2$dataAbund <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences %in% figure2$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance"
)

figure2$indices <- figure2$graph$networks$Index |> tidytable::filter(
  SpeciesPreferences == figure2$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds
) |> tidytable::arrange(
  Intervention
)

##### a: ######################################################################
# Richness through time across simulations, showing stability and separation.
figure2$plotA <- plotMeanAndInner(
  rbind(
    figure2$dataRich |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)"),
      is.na(Subset)
    ),
    # We want to appear in the legend but not on the plot!
    figure2$dataRich |> tidytable::filter(
      PoolPatchSeed == figure2$dataRich$PoolPatchSeed[1],
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - 10000) == min(abs(Time - 10000)),
      is.na(Subset)
    ) |> tidytable::mutate(
      Value = 10 # coord_cartesian will eliminate these points.
    )
  ) |> tidytable::mutate(
    SpeciesPreferences = tidytable::case_when(
      SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
      SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
      SpeciesPreferences == "Uniform(0, 1)" ~ "Multiple Adaptation Types",
      TRUE ~ SpeciesPreferences
    )
  ), CIs = 0.75
) + ggplot2::geom_point(
  data = function(x) {x |> tidytable::filter(
    PoolPatchSeed == figure2$graph$seed,
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    abs(Time - figure2$graph$time) == min(abs(Time - figure2$graph$time))
  )},
  mapping = ggplot2::aes(fill = Intervention),
  shape = 21,
  color = "black"
)  + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5, title = "Habitat Type"),
  fill = ggplot2::guide_legend(ncol = 5, title = "Habitat Type")
) + ggplot2::theme(
  legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95),
  panel.spacing = ggplot2::unit(1, "lines"),
  strip.text = ggplot2::element_text(size = 12)
) + ggplot2::coord_cartesian(
  xlim = c(0, 31000), ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
) + ggplot2::facet_grid(
  # switch = "y",
  cols = ggplot2::vars(SpeciesPreferences)
)

##### Networks: ###############################################################
# Example networks from different scenarios of the same simulation, showing
# effects of the current habitat type through time on network shape.
# Previously, these were independent panels, but I'm switching to a facets.
# figure2$plotNetworks <- figure2$graph$networks$Plot + ggplot2::facet_grid(
#   # Reverse order
#   factor(Intervention, levels = c("(1)", "(0.5)", "(0)"), ordered = T) ~ .
# ) + ggplot2::theme(
#   axis.title.x = ggplot2::element_blank(),
#   panel.border = ggplot2::element_rect(color = "black", fill = NA),
#   legend.position = "none",
#   axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = -30),
#                                       size = 12),
#   axis.text.x = ggplot2::element_blank()
# ) + ggplot2::coord_cartesian(xlim = c(-0.25, 1))
figure2$plotNetworks <- dplyr::bind_rows(lapply(
  figure2$graph$networks$Envs, function(env) {
    value <- if (env$Row$Intervention == "(0)") {
      10
    } else if (env$Row$Intervention == "(0.5)") {
      22.5
    } else if (env$Row$Intervention == "(1)") {
      35
    } else {
      {warning("env$Row$Intervention not handled."); NA}
    }
    tibble::tibble(
      grob = list(
        (env$singletonGraphs[[1]] + ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          legend.position = "none",
          axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = -23),
                                              size = 9, vjust = 0.2),
          axis.text.x = ggplot2::element_blank()
        ) + ggplot2::coord_cartesian(
          xlim = c(-0.25, 1), ylim = c(0.01, 3)
        ) + ggplot2::ylab(
          ""
        ) + ggplot2::scale_size(
          range = c(0.1, 2)
        )) |> ggplot2::ggplotGrob() 
      ),
      Time = 20000, 
      Value = value
    )
  } 
))

##### b: Violins ##############################################################
figure2$plotB <- ggplot2::ggplot(
  figure2$dataRich |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences, Subset
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention,
    group = paste(Intervention, Subset)
  )
  # OVERALL Violins
) + ggplot2::geom_violin(
  data = function(x) x |> tidytable::filter(is.na(Subset)),
  position = ggplot2::position_dodge(0.9), scale = "count"
) + ggplot2::geom_boxplot(
  data = function(x) x |> tidytable::filter(is.na(Subset)),
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
  # BASAL Violin
) + geom_flat_violin(
  data = function(x) x |> tidytable::filter(
    !is.na(Subset)
  ) |> tidytable::separate_wider_delim(
    cols = Subset, delim = "_", names = c("Subset", "Bins")
  ) |> tidytable::filter(
    Subset == "Basal"
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences, Subset
  ) |> tidytable::summarise(
    Value = sum(Value)
  ),
  flip = 1, color = "black",
  ggplot2::aes(fill = Intervention)
  # CONSUMER Violin
) + geom_flat_violin(
  data = function(x) x |> tidytable::filter(
    !is.na(Subset)
  ) |> tidytable::separate_wider_delim(
    cols = Subset, delim = "_", names = c("Subset", "Bins")
  ) |> tidytable::filter(
    Subset == "Consumer"
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences, Subset
  ) |> tidytable::summarise(
    Value = sum(Value)
  ),
  flip = 2, color = "black",
  ggplot2::aes(fill = Intervention)
  # LABELS
) + ggplot2::geom_text(
  data = function(x) {
    x |> tidytable::filter(
      !is.na(Subset)
    ) |> tidytable::separate_wider_delim(
      cols = Subset, delim = "_", names = c("Subset", "Bins")
    ) |> tidytable::group_by(
      PoolPatchSeed, Intervention, SpeciesPreferences, Subset
    ) |> tidytable::summarise(
      Value = sum(Value),
      .groups = "drop"
    ) |> tidytable::mutate(
      Offset = (max(Value) - min(Value)) * 0.08
    ) |> dplyr::group_by(
      Intervention, SpeciesPreferences, Subset
    ) |> dplyr::summarise(
      Value = min(Value - Offset),
      Value = ifelse(Value > 0, Value, 5),
      Label = substr(Subset[1], 0, 1),
      .groups = "drop"
    )
  },
  mapping = ggplot2::aes(label = Label),
  position = ggplot2::position_dodge(0.9)
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::coord_cartesian(
  ylim = c(0, richnessYMax), expand = FALSE
) + ggplot2::facet_grid(
  . ~ SpeciesPreferences
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  panel.grid.minor = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  strip.text = ggplot2::element_blank(),
  panel.spacing = ggplot2::unit(1, "lines")
) + ggplot2::labs(
  y = "Richness",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
)

figure2$plotC <- ggplot2::ggplot(
  figure2$dataAbund |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences, Subset
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention,
    group = paste(Intervention, Subset)
  )
  # OVERALL Violins
) + ggplot2::geom_violin(
  data = function(x) x |> tidytable::filter(is.na(Subset)),
  position = ggplot2::position_dodge(0.9), scale = "count"
) + ggplot2::geom_boxplot(
  data = function(x) x |> tidytable::filter(is.na(Subset)),
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
  # BASAL Violins
) + geom_flat_violin(
  data = function(x) x |> tidytable::filter(
    !is.na(Subset)
  ) |> tidytable::separate_wider_delim(
    cols = Subset, delim = "_", names = c("Subset", "Bins")
  ) |> tidytable::filter(
    Subset == "Basal"
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences, Subset
  ) |> tidytable::summarise(
    Value = sum(Value)
  ),
  flip = 1, color = "black",
  ggplot2::aes(fill = Intervention)
  # CONSUMER Violins
) + geom_flat_violin(
  data = function(x) x |> tidytable::filter(
    !is.na(Subset)
  ) |> tidytable::separate_wider_delim(
    cols = Subset, delim = "_", names = c("Subset", "Bins")
  ) |> tidytable::filter(
    Subset == "Consumer"
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesPreferences, Subset
  ) |> tidytable::summarise(
    Value = sum(Value)
  ),
  flip = 2, color = "black",
  ggplot2::aes(fill = Intervention)
  # LABELS
) + ggplot2::geom_text(
  data = function(x) {
    x |> tidytable::filter(
      !is.na(Subset)
    ) |> tidytable::separate_wider_delim(
      cols = Subset, delim = "_", names = c("Subset", "Bins")
    ) |> tidytable::group_by(
      PoolPatchSeed, Intervention, SpeciesPreferences, Subset
    ) |> tidytable::summarise(
      Value = sum(Value),
      .groups = "drop"
    ) |> tidytable::mutate(
      Offset = (max(log10(Value)) - min(log10(Value))) * 0.04
    ) |> dplyr::group_by(
      Intervention, SpeciesPreferences, Subset
    ) |> dplyr::summarise(
      Value = 10^min(log10(Value) - Offset),
      Label = substr(Subset[1], 0, 1),
      .groups = "drop"
    )
  },
  mapping = ggplot2::aes(label = Label),
  position = ggplot2::position_dodge(0.9)
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
# ) + ggplot2::facet_grid(
#   !(is.na(Subset) | Subset == "Basal") ~ SpeciesPreferences,
#   scales = "free"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  panel.grid.minor = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  strip.text = ggplot2::element_blank(),
  panel.spacing = ggplot2::unit(0.25, "lines")
) + ggplot2::labs(
  y = "Abundance",
  x = "Habitat Type"
) + ggplot2::guides(
  color = "none",
  fill = "none"
# ) + ggplot2::scale_y_log10(
) + ggplot2::scale_y_continuous(
  transform = "pseudo_log", breaks = 10^(0:4),
  label = scales::label_log(digits = 2)
)

##### Combine: ################################################################
figure2$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure2$plotA + ggpp::geom_grob(
      data = figure2$plotNetworks,
      mapping = ggplot2::aes(x = Time, y = Value, label = grob),
      vp.height = 0.35, vp.width = 0.3, default.alpha = 0.25
    ),
    cowplot::plot_grid(
      plotlist = list(
        figure2$plotB,
        figure2$plotC
      ), nrow = 2
    )
  ),
  ncol = 1, heights = c(2/5, 3/5), common.legend = TRUE
)

ggplot2::ggsave(plot = figure2$plot, filename = file.path(dirImages, "Figure2_Prototype2.pdf"),
                units = "cm", width = 6.5*2, height = 6.5*3)
ggplot2::ggsave(plot = figure2$plot, filename = file.path(dirImages, "Figure2_Prototype2.png"),
                units = "cm", width = 6.5*2, height = 6.5*3)
ggplot2::ggsave(plot = figure2$plotA, filename = file.path(dirImages, "Figure2A_Prototype2.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
# ggplot2::ggsave(plot = figure2$plotNetworks,
#                 filename = file.path(dirImages, "Figure2Networks_Prototype1.pdf"),
#                 units = "cm", width = 6.5*1, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotB, filename = file.path(dirImages, "Figure2B_Prototype2.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)
ggplot2::ggsave(plot = figure2$plotC, filename = file.path(dirImages, "Figure2C_Prototype2.pdf"),
                units = "cm", width = 6.5*3, height = 6.5*2)

