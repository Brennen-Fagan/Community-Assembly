### Plot 5: ###################################################################
# Progression of network change as we undergo intervention. As a base plot
# we use the richness changes of two experiments that are trading places.
# Both of these should probably be 100% 0. Maybe (0) -> (0.5) and (0.5) -> (0).
# Because of the size of the plots, we can't actually show them along the
# richness plots directly. We could potentially put labeled points instead,
# and then time staggered facets. I think we'll need higher resolution evals
# in order to capture the level of detail we're describing in the main text.
# We also need to convert the existing code for creating the networks into
# more general code, since we'll need to make a few here as well...

# Try something like the below to identify a good option:
# diversitiesRichness |> tidytable::filter(
#   NicheDistance == defaultNicheDistance,
#   (PoolPatchSeed %in% as.character(383)),#:386)),
#   Metric == "Alpha Hill:0",
#   SpeciesAffinity == "100% 0",
#   Intervention %in% c("(0)->(0.5)", "(0.5)->(0)", "(0)", "(0.5)"),
#   Time > 16000, Time < 18000,
#   is.na(Subset)
# ) |> ggplot(
#   aes(x = Time, y = Value, color = Intervention)
# ) + ggplot2::geom_line(
# ) + ggplot2::facet_wrap(
#   PoolPatchSeed ~ .
# ) + coord_cartesian(
#   xlim = c(16600, 17300)
# )

newplot5_a_Specification <- rbind(diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(383)),#:386)),
  Metric == "Alpha Hill:0",
  SpeciesAffinity == "100% 0",
  Intervention %in% c("(0)", "(0.5)"),
  Time %in% c(16700),
  is.na(Subset)
), diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(383)),#:386)),
  Metric == "Alpha Hill:0",
  SpeciesAffinity == "100% 0",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)"),
  Time %in% c(16720, 16800, 16900, 17100),
  is.na(Subset)
))

newplot5_a_Networks <- generateNetworks(newplot5_a_Specification)

newplot5_a <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(383)),#:386)),
  Metric == "Alpha Hill:0",
  SpeciesAffinity == "100% 0",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)", "(0)", "(0.5)"),
  Time > 16000, Time < 18000,
  is.na(Subset)
) |> ggplot(
  aes(x = Time, y = Value, color = Intervention)
) + ggplot2::geom_line(
  # show.legend = FALSE
) + coord_cartesian(
  xlim = c(16600, 17300),
  ylim = c(6, 28)
) + ggplot2::geom_rect(
  data = newplot5_a_Specification |> tidytable::group_by(
    Time
  ) |> tidytable::summarise(
    xmin = Time - 2, xmax = Time + 2,
    ymin = 5, ymax = 27
  ),
  mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  color = "grey", inherit.aes = FALSE, alpha = 0.2
) + ggplot2::geom_point(
  show.legend = FALSE,
  data = newplot5_a_Specification
) + ggplot2::geom_vline(
  xintercept = (
    newplot5_a_Networks$Envs[[3]]$result$Ellipsis$Affinity$TimeIntervention
    / newplot5_a_Networks$Envs[[3]]$result$ReactionTime
  ), linetype = "dashed"
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

# For some reason, this is returning a list of two plots rather than a single
# plot when used with ncol or nrow.
newplot5_b <- ggarrange(plotlist = ggarrange(
  plotlist = lapply(
    seq_along(newplot5_a_Networks$Envs),
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
    e = newplot5_a_Networks$Envs[c(1, 3:6, 2, 7:10)]
  ),
  ncol = 5
), nrow = 2, labels = list("b)", "c)"),
font.label = list(face = "plain"),
vjust = 1.4, hjust = 0)


newplot5 <- ggpubr::ggarrange(
  plotlist = list(
    newplot5_a,
    newplot5_b
  ), nrow = 2, heights = c(0.2/0.9, 0.7/0.9)
)

# Decorate with additional indicators.
newplot5 <-
  newplot5 + ggplot2::annotate(
    "curve", x = 0.2, y = 0.865, xend = 0.1, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
  ) + ggplot2::annotate(
    "curve", x = 0.23, y = 0.865, xend = 0.3, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = -0.3
  ) + ggplot2::annotate(
    "curve", x = 0.33, y = 0.865, xend = 0.5, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = -0.3
  ) + ggplot2::annotate(
    "curve", x = 0.45, y = 0.865, xend = 0.7, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = -0.3
  ) + ggplot2::annotate(
    "curve", x = 0.7, y = 0.865, xend = 0.9, yend = 0.73,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
    curvature = -0.2
  ) + ggplot2::annotate(
    "segment", linetype = "dashed", x = 0.215, y = 0.865, xend = 0.2, yend = 0.745
  ) + ggplot2::annotate(
    "segment", linetype = "dashed", x = 0.2, y = 0.745, xend = 0.2, yend = 0
  )

ggplot2::ggsave(plot = newplot5, filename = "Figure5_Prototype2.png",
                units = "cm", width = 6.5*3, height = 6.5*2)