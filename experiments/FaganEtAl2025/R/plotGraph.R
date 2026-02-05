plotGraph <- function(graph, mainLayout, legends = FALSE) {
  obj <- ggraph::ggraph(
    graph = graph %N>% mutate(
      nodesize = (log10(N)+5)/10+1
    ),
    layout = graph %N>% data.frame(
    ) |> select(node) |> left_join(
      mainLayout |> select(x, y, node)
    )
  ) + ggraph::geom_edge_diagonal(
    mapping = aes(
      color = Type,
      linetype = Type,
      alpha = log10(effectNormalised),
      end_cap = circle(node2.nodesize+2, 'pt')
    ),
    arrow = arrow(length = unit(2, 'mm')),
    show.legend = legends
  ) + ggraph::geom_node_point(
    mapping = aes(
      color = Type,
      size = nodesize
    ),
    show.legend = legends
    # ) + ggplot2::geom_hline(
    #   yintercept = -1, linetype = "dashed", color = "black"
  ) + theme_minimal(
  ) + ylab(
    "Size" # "Log10(Size)"
  ) + xlab(
    "Land-use Preference"
  ) + scale_color_manual(
    values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
  ) + ggraph::scale_edge_color_manual(
    values =
      c("Exploit+" = "darkgreen", "Exploit-" = "burlywood4")
    # c("Basal" = "darkgreen", "Consumer" = "burlywood4")
  ) + ggplot2::scale_y_log10(
  ) + scale_size(
    range = c(0.5, 4) # size of the resulting points, not the inputs.
    # limits = c(10^-5, 10^5)#, trans = "log10"
  ) + coord_cartesian(
    x = c(0, 1), y = 10^c(-2, 0.5),
    clip = "off"
  )

  # if (!legends) {
  #   obj <- obj + ggplot2::theme(
  #     legend.position = "none"
  #   )
  # }

  return(obj)
}
