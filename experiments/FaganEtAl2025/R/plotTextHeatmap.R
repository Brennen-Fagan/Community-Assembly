##### FUNCTION: ###############################################################
# Different Scales mean we have to separate out the data, so we define a
# function to perform the plotting repeatedly/consistently.
plotTextHeatmap <- function(
  data, legendName,
  legendtrans = "identity",
  scalestrans = scales::label_percent(accuracy = 0.1)
) {
  cutbreaks <- with(
    data[is.finite(data$Average),], {
      if (!any(Average < 0) && mean(Average > 0 & Average < 2) > 0.25) {
        # Relative/Multiplicative Scale likely, centre around 1
        c(0,
          round(10^c(-0.5, -0.1, 0.1, 0.5), digits = 1),
          Inf)
      } else if (any(Average < 0) && mean(Average > -1 & Average < 1) > 0.25) {
        # Relative/Multiplicative scale centred around 0 likely (% difference)
        c(-1,
          round(10^c(-0.3, -0.05, 0.05, 0.3), digits = 1)-1,
          Inf)
      } else { # Raw values
        c(floor(min(Average)),
          quantile(Average, p = c(0.2, 0.4, 0.6, 0.8)),
          ceiling(max(Average)))
      }
    }
  )
  cutbreaksLabels <- scalestrans(cutbreaks)
  cutbreaksLabels <- paste(cutbreaksLabels[1:(length(cutbreaksLabels) - 1)],
                           cutbreaksLabels[2:length(cutbreaksLabels)],
                           sep = " - ")

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = InterventionInitial,
      y = InterventionFinal,
      fill = factor(findInterval(Average, vec = cutbreaks, all.inside = TRUE),
                    levels = 1:length(cutbreaksLabels),
                    labels = cutbreaksLabels)
    )
  ) + ggplot2::geom_tile(
    width = 1, height = 1, color = NA
  ) + ggplot2::geom_tile(
    data = function(x) x |> tidytable::filter(Emphasis),
    fill = NA, color = "black", linewidth = 1
  ) + ggplot2::geom_text(
    ggplot2::aes(label = scalestrans(Average))
    # ) + ggplot2::facet_grid(#
  ) + ggh4x::facet_nested_wrap(
    . ~ Metric + Time
  ) + ggplot2::theme_minimal(
  ) + ggplot2::theme(
    strip.text = ggplot2::element_text(size = 12)
  ) + ggplot2::labs(
    fill = legendName,
    x = "Initial Habitat Type",
    y = "Final Habitat Type"
  ) + ggplot2::coord_fixed(
  )
}
