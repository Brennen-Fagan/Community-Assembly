CreateDispersalMatrix <- function(
  EnvironmentDistances, # Distances, which we will invert. Square Matrix >= 0.
  # Entries are zero if two environments are not connected. No self-loops.
  SpeciesSpeeds # Distances/Times.
  # Result is a matrix of frequencies (1/Times) which will be %*% Abundance.
) {
  NumEnvironments <- nrow(EnvironmentDistances)
  stopifnot(NumEnvironments == ncol(EnvironmentDistances))
  stopifnot(
    tryCatch(diag(EnvironmentDistances), error = function(e) {
      Matrix::diag(EnvironmentDistances)
    }) == 0
  )
  
  dispersalDiags <- NULL
  # Take i to be the (super/sub) diagonal index, aka "band".
  for (i in (NumEnvironments - 1):-(NumEnvironments - 1)) {
    if (i == 0) next
    # Start in top right band, move to bottom left.
    dispersalDiags <- c(
      dispersalDiags,
      list(c(
        unlist(lapply(
          1:NumEnvironments,
          function(index, offset, mat, vec) {
            if (offset != 0 &&
                nrow(mat) >= index + offset &&
                index + offset >= 1
            )
              if (mat[index, index + offset] != 0) { # To prevent Diffusive Infs.
                1/mat[index, index + offset] * vec
              } else {
                0 * vec
              }
          },
          offset = i, mat = EnvironmentDistances, vec = SpeciesSpeeds
        ))
      ))
    )
  }
  
  # Amount of travel from j to i is d[i,j]
  # Amount of gain to i from j is d[i,j]
  # Amount of travel from i is d[i,i]
  # This way, we can write the change in y from travel as d %*% y
  # Characterising this as proportions of a population, but that is assuming a
  # normalisation that I do not think is strictly necessary.
  # The matrix is sparse and has colsum = 0, diag < 0, offdiag >= 0.
  ns <- length(SpeciesSpeeds) * NumEnvironments
  ks <- length(SpeciesSpeeds) * c((NumEnvironments - 1):1,
                                  -(1:(NumEnvironments - 1)))
  ks <- ks[!unlist(lapply(dispersalDiags, is.null))]
  dispersalDiags <- dispersalDiags[!unlist(lapply(dispersalDiags, is.null))]
  
  dispersalMatrix <- Matrix::bandSparse(
    n = ns,
    k = ks,
    diagonals = dispersalDiags# c(
    # list(c(rep(0.0001, length(redCom)),   # Island 2 -> Island 1
    #        rep(0.0001, length(redCom)))), # Island 3 -> Island 2
    # list(rep(0, length(redCom))),         # Island 3 -> Island 1
    # list(c(rep(0.0001, length(redCom)),   # Island 1 -> Island 2
    #        rep(0.0001, length(redCom)))), # Island 2 -> Island 3
    # list(rep(0, length(redCom)))          # Island 1 -> Island 3
    # )
  )
  stopifnot(all(dispersalMatrix >= 0))
  diag(dispersalMatrix) <- -Matrix::colSums(dispersalMatrix)
  stopifnot(isTRUE(all.equal(Matrix::colSums(dispersalMatrix),
                             rep(0, ncol(dispersalMatrix)))))
  
  return(dispersalMatrix)
}