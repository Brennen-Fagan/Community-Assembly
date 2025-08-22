# Setup: ######################################################################

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsTimeBC.R")

### 11 Supplement: ############################################################
##### Turnover related statistics: ############################################
# Problem here: data is bimodal, so central tendency isn't quite catching the
# right information.

load("diversitiesFlattened9a9_subset4TimeBC.RData")

# newplot5_dataB <- diversitiesAll %>% newplot5_filtration(
# ) %>% tidytable::filter(
#   Metric == "TimeBrayCurtis",
#   is.na(Subset)
# ) %>% tidytable::left_join(endTimes %>% dplyr::select(-Times))
#
# newplot5_ba <- plotMeanAndInner(
#   newplot5_dataB, CIs = 0.75, facets = as.formula(. ~ .)
#   # ) + ggplot2::geom_rect(
#   #   data = data.frame(
#   #     1 # 1 rectangle per row, so dummy df to prevent overplotting
#   #   ),
#   #   xmin = min(newplot5_dataB$Start),
#   #   xmax = max(newplot5_dataB$Stop),
#   #   ymin = 0, ymax = max(newplot5_dataB$Value, na.rm = TRUE),
#   #   fill = "grey",
#   #   alpha = 0.2,
#   #   inherit.aes = FALSE
# ) + ggplot2::labs(
#   y = "Time-BC",
#   tag = "b)"
# ) + ggplot2::theme(
#   legend.position = "none",
#   plot.tag.position = c(0.05, 1)
# ) + ggplot2::coord_cartesian(
#   xlim = c(0, 30000), ylim = c(0, 0.5), expand = FALSE
# )
# newplot5_bb <- ggplot2::ggplot(
#   newplot5_dataB %>% tidytable::filter(
#     Time > Start, Time < Stop
#   ),
#   ggplot2::aes(
#     y = Value, color = Intervention
#   )
# ) + ggplot2::geom_density(
#   adjust = 1/2, show.legend = FALSE
# ) + ggplot2::coord_cartesian(
#   ylim = c(0, 0.5), expand = FALSE
# ) + ggplot2::theme_minimal(
# ) + ggplot2::scale_color_manual(
#   values = colorPalette, aesthetics = c("color", "fill"),
#   name = "Habitat Land-use"
# ) + ggplot2::scale_x_continuous(
#   name = "", labels = function(x) rep("", length(x))
# ) + ggplot2::theme(
#   axis.title.y = ggplot2::element_blank(),
#   axis.text.y = ggplot2::element_blank(),
#   axis.ticks.y = ggplot2::element_blank()#,
#   # axis.text.x = ggplot2::element_blank()
# )
# newplot5_b <- ggpubr::ggarrange(
#   newplot5_ba, newplot5_bb,
#   nrow = 1, widths = c(0.9, 0.1), align = "h"
# )
