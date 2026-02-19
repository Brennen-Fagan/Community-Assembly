
# https://stackoverflow.com/a/52035174
# which cites https://gist.github.com/dgrtwo
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             flip = 1,
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      flip = flip,
      ...
    )
  )
}


GeomFlatViolinSD <- function(data, params) {
  data$width <- data$width %||%
    params$width %||% (resolution(data$x, FALSE) * 0.9)

  # print(params)
  # print(str(data))

  # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
  data %>%
    group_by(group) %>%
    mutate(ymin = min(y),
           ymax = max(y),
           xmin = x - (params$flip == 1) * width / 2,
           xmax = x + (params$flip == 2) * width / 2)
}

GeomFlatViolinDG <- function(data, panel_scales, coord) {
  # print(str(data))
  # print(panel_scales)
  # print(coord)

  # Find the points for the line to go all the way around
  data <- transform(data,
                    xmaxv = x + violinwidth * (flip == 2) * (xmax - x),
                    xminv = x + violinwidth * (flip == 1) * (xmin - x))

  # Make sure it's sorted properly to draw the outline
  newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                   plyr::arrange(transform(data, x = xmaxv), -y))

  # Close the polygon: set first and last point the same
  # Needed for coord_polar and such
  newdata <- rbind(newdata, newdata[1,])

  ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = GeomFlatViolinSD,
            # function(data, params) {
            #   data$width <- data$width %||%
            #     params$width %||% (resolution(data$x, FALSE) * 0.9)
            #
            #   # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            #   data %>%
            #     group_by(group) %>%
            #     mutate(ymin = min(y),
            #            ymax = max(y),
            #            xmin = x - width / 2,
            #            xmax = x)
            # },

          draw_group = GeomFlatViolinDG,
            # function(data, panel_scales, coord) {
            #   # Find the points for the line to go all the way around
            #   data <- transform(data,
            #                     xmaxv = x,
            #                     xminv = x + violinwidth * (xmin - x))
            #
            #   # Make sure it's sorted properly to draw the outline
            #   newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
            #                    plyr::arrange(transform(data, x = xmaxv), -y))
            #
            #   # Close the polygon: set first and last point the same
            #   # Needed for coord_polar and such
            #   newdata <- rbind(newdata, newdata[1,])
            #
            #   ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
            # },

          draw_key = draw_key_polygon,

          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid", flip = 1),

          required_aes = c("x", "y")
  )
