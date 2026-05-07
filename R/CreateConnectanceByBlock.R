CreateConnectanceByBlock <- function(
  BlockSizes = c(34, 66),
  Connectances = matrix(c(1, 0.2, NA, 1),
                        byrow = TRUE, nrow = 2),
  Directed = FALSE
) {
  stopifnot(length(BlockSizes) == ncol(Connectances),
            length(BlockSizes) == nrow(Connectances))
  # Won't use lower-triangle if !Directed, hence the defaults.

  adjmat <- Matrix::Matrix(nrow = sum(BlockSizes), ncol = sum(BlockSizes))

  indices <- cumsum(c(0, BlockSizes))
  for (i in seq_along(BlockSizes)) {
    is <- (1 + indices[i]):indices[i + 1]
    for (j in seq_along(BlockSizes)) {
      if (!Directed && j < i) next()

      js <- (1 + indices[j]):indices[j + 1]

      if (i != j) {
        graph <- igraph::as_adjacency_matrix(igraph::bipartite.random.game(
          BlockSizes[i], BlockSizes[j], p = Connectances[i, j]
        ))[1:BlockSizes[i], BlockSizes[i] + 1:BlockSizes[j]]
      } else {
        graph <- igraph::as_adjacency_matrix(igraph::random.graph.game(
          BlockSizes[i], p = Connectances[i, j]
        ))
      }

      adjmat[is, js] <- graph
      if(!Directed) adjmat[js, is] <- Matrix::t(graph)
    }
  }

  return(adjmat)
}
