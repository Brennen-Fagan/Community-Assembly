# Setup: ######################################################################

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")

### 3 Supplement: #############################################################
newplot3_bs <- ggplot2::ggplot(
  newplot3_dataB,
  ggplot2::aes(
    y = Persistence,
    x = Intervention,
    color = Intervention,
    group = interaction(Intervention, SpeciesType, AffinityBins),
    fill = SpeciesType
  )
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = -Inf, ymax = Inf,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9), show.legend = FALSE,
  scale = "count", draw_quantiles = 0.5
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("darkgreen", "burlywood4")
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::facet_grid(
  factor(SpeciesType, levels = c("Consumer", "Basal"), ordered = TRUE) ~
    SpeciesPreferences
) + ggplot2::labs(
  x = "Habitat's Land-use,\nsubdivided by Species Land-use Preference"
) + ggplot2::geom_text(
  data = rbind(
    data.frame(
      x = c(1:5 - 0.22, 1:5 + 0.22),
      y = 12000,
      lab = c(rep("0", 5), rep("1", 5)),
      SpeciesPreferences = "50% 0, 50% 1"
    ),
    data.frame(
      # Approximately the "right" spacing when blown up to large scales...
      x = c(1:5 - 0.36, 1:5 - 0.18, 1:5 - 0, 1:5 + 0.18, 1:5 + 0.36),
      y = 12000,
      lab = c(rep("0.1", 5), rep("0.3", 5), rep("0.5", 5),
              rep("0.7", 5), rep("0.9", 5)),
      SpeciesPreferences = "Uniform(0, 1)"
    )
  ),
  inherit.aes = FALSE,
  mapping = ggplot2::aes(
    x = x, y = y, label = lab
  )
)

# newplot3 <- ggpubr::ggarrange(
#   plotlist = list(
#     newplot3_a,
#     newplot3_b
#   ), nrow = 2, widths = c(0.5, 0.5)
# )

ggplot2::ggsave(plot = newplot3_bs, filename = "Figure3s1_Prototype2.png",
                units = "cm", width = 6.5*6, height = 6.5*2)
