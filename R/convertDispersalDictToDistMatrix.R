convertDispersalDictToDistMatrix <- function(dispersalDictionary, nEnv) {
  with(c(
    dispersalDictionary,
    Environments = nEnv
  ), {
    if (Configuration == "None") {
      DistanceMatrix <- Matrix::sparseMatrix(
        i = Environments, j = Environments, x = 0)
    }
    if (Configuration == "Ring" || Configuration == "Line")
      DistanceMatrix <- Matrix::bandSparse(
        Environments, k = c(-1, 1),
        diagonals = list(rep(Resistance, Environments - 1),
                         rep(Resistance, Environments - 1))
      )
    if (Configuration == "Ring") {
      DistanceMatrix[Environments, 1] <- Resistance
      DistanceMatrix[1, Environments] <- Resistance
    }
    if (Configuration == "Complete") {
      DistanceMatrix <- matrix(Resistance,
                               nrow = Environments,
                               ncol = Environments)
      diag(DistanceMatrix) <- 0
    }
    return(DistanceMatrix)
  }
  )
}
