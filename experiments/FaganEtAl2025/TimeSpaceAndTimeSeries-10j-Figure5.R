# Setup: ######################################################################
# Progression of network change as we undergo intervention. For orientation,
# we use the richness changes of two experiments that are trading places.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source(file.path("R", "flattenDiversity.R")) # Req'd by below
source(file.path("R", "generateNetworks.R")) # To create inset graphs.
load("TSTS_Interventions_10a1.RData")

# This is better as an environment, but that's more opaque.
figure5 <- list(
  graph = list(
    seed = "27" # "7", "9", "15", "27", "35", "38", "41" are clearest.
  )
)
# If you want weird, try "25", "30", "43", "44"
figure5$graph$timeIntervention <- with(
  InterventionTimes |> tidytable::filter(
    PoolPatchSeed == figure5$graph$seed,
    NicheDistance == defaultNicheDistance,
    SpeciesAffinity == "1"
  ),
  unique(TimeIntervention)
)
stopifnot(length(figure5$graph$timeIntervention) == 1)

# Guaranteed Options: Intervention + c(0:10, 2*(6:10), 20+3*(1:10))
#                     round(Intervention, digits = -1) + multiples of 10
figure5$graph$time <- c(
  floor(figure5$graph$timeIntervention/10)*10, # Just before!
  ceiling(figure5$graph$timeIntervention/10)*10 + c(20, 200, 500, 900) # After
)

# # Try something like the below to identify a good option:
# diversitiesRichness |> tidytable::filter(
#   NicheDistance == defaultNicheDistance,
#   # (PoolPatchSeed %in% as.character(19)),#:386)),
#   Metric == "Alpha Hill:0",
#   SpeciesPreferences == "100% 0",
#   Intervention %in% c("(0)->(0.5)", "(0.5)->(0)", "(0)", "(0.5)"),
#   Time > 15700, Time < 18300,
#   is.na(Subset)
# ) |> ggplot(
#   aes(x = Time, y = Value, color = Intervention)
# ) + ggplot2::geom_line(
# ) + ggplot2::facet_wrap(
#   PoolPatchSeed ~ .
# ) + coord_cartesian(
#   xlim = c(16000, 18000)
# )

# Main Plots: #################################################################
### Plot 5: ###################################################################

figure5$graph$specification <- rbind(diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% figure5$graph$seed,
  Time == figure5$graph$time[1],
  Metric == "Alpha Hill:0",
  SpeciesPreferences == "100% 0",
  Intervention %in% c("(0)", "(0.5)"),
  is.na(Subset)
), diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% figure5$graph$seed,
  round(Time, 1) %in% round(figure5$graph$time[-1], 1),
  Metric == "Alpha Hill:0",
  SpeciesPreferences == "100% 0",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)"),
  is.na(Subset)
))

figure5$graph$networks <- generateNetworks(figure5$graph$specification)

##### a: ######################################################################
figure5$plotA <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed == figure5$graph$seed,
  Metric == "Alpha Hill:0",
  SpeciesPreferences == "100% 0",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)", "(0)", "(0.5)"),
  Time >= floor(min(figure5$graph$time)/100)*100 - 100,
  Time <= ceiling(max(figure5$graph$time)/100)*100 + 100,
  is.na(Subset)
) |> ggplot(
  aes(x = Time, y = Value, color = Intervention)
) + ggplot2::geom_line(
) + coord_cartesian(
  xlim = c(
    floor(min(figure5$graph$time)/100)*100 - 10,
    ceiling(max(figure5$graph$time)/100)*100 + 10
  )#,
  # ylim = c(6, 28)
) + ggplot2::geom_rect(
  data = figure5$graph$specification |> tidytable::mutate(
    Time = round(Time)
  ) |> tidytable::group_by(
    Time
  ) |> tidytable::summarise(
    xmin = Time - 2, xmax = Time + 2,
    ymin = 5, ymax = 27
  ),
  mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  color = "grey", inherit.aes = FALSE, alpha = 0.2
) + ggplot2::geom_point(
  show.legend = FALSE,
  data = figure5$graph$specification
) + ggplot2::geom_vline(
  xintercept = figure5$graph$timeIntervention, linetype = "dashed"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = ""#"Habitat's Land-use"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 4),
  fill = ggplot2::guide_legend(ncol = 4)
) + ggplot2::labs(
  tag = "a)"
  #   y = "Richness "
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  legend.position = c(0.5, 1.2),
  plot.tag.position = c(0, 1),
  axis.title.y = ggplot2::element_blank()
)

##### b: ######################################################################
# For some reason, this is returning a list of two plots rather than a single
# plot when used with ncol or nrow.
figure5$plotB <- ggarrange(plotlist = ggarrange(
  plotlist = lapply(
    seq_along(figure5$graph$networks$Envs),
    function(i, e) {
      e <- e[[i]]
      g <- e$singletonGraphs[[1]] + ggplot2::theme_void(
      ) + ggplot2::annotate(
        "text", x = -0.05, y = 0.5, label = e$Row$Time, hjust = 0
      ) + ggplot2::theme(
        panel.border = ggplot2::element_rect(
          # Ordinal in different order than color palette, and indexing
          # defaults to as.numeric for ordinal rather than as.character.
          color = colorPalette[as.character(e$Row$Intervention)],
          fill = NA
        ),
        plot.title = ggplot2::element_text(size = 13, hjust = 1)
      ) + ggplot2::ggtitle(
        if (i %in% c(1, 2, 6, 7)) {e$Row$Intervention} else {""}
      )

      return(g)
    },
    e = figure5$graph$networks$Envs[c(1, 3:6, 2, 7:10)]
  ),
  ncol = 5
), nrow = 2, labels = list("b)", "c)"),
font.label = list(face = "plain"),
vjust = 1.4, hjust = 0)

figure5$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure5$plotA,
    figure5$plotB
  ), nrow = 2, heights = c(0.2/0.9, 0.7/0.9)
)

# Decorate with additional indicators.
figure5$plot <-
  # Arrows to the corresponding panel. # MANUAL CORRECTIONS.
  figure5$plot + ggplot2::annotate(
    "curve", x = 0.1675, y = 0.885, xend = 0.1, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = 0.1
  ) + ggplot2::annotate(
    "curve", x = 0.19, y = 0.885, xend = 0.3, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = -0.3
  ) + ggplot2::annotate(
    "curve", x = 0.34, y = 0.885, xend = 0.5, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = -0.3
  ) + ggplot2::annotate(
    "curve", x = 0.6, y = 0.885, xend = 0.7, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = -0.3
  ) + ggplot2::annotate(
    "curve", x = 0.94, y = 0.885, xend = 0.9, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = 0.3
    # Intervention Line
  ) + ggplot2::annotate(
    "segment", linetype = "dashed",
    x = 0.17, y = 0.885, xend = 0.2, yend = 0.745
  ) + ggplot2::annotate(
    "segment", linetype = "dashed",
    x = 0.2, y = 0.745, xend = 0.2, yend = 0
  )

ggplot2::ggsave(plot = figure5$plot, filename = "Figure5_Prototype3.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

