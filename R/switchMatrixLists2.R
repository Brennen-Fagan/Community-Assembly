switchMatrixLists2 <- function(matrixlist1, matrixlist2, switchtime) {
  # Note: only accepts 1 switchtime.
  #       multiple would require steps.
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(length(matrixlist1) == length(switchtime) ||
              length(switchtime) == 1)
  stopifnot(isTRUE(all.equal(lapply(matrixlist1, dim),
                             lapply(matrixlist2, dim))))

  matrix1 <- Matrix::bdiag(matrixlist1)
  matrix2 <- Matrix::bdiag(matrixlist2)

  switchMatrices(matrix1, matrix2, switchtime = switchtime)
}
