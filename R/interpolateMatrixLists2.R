interpolateMatrixLists2 <- function(matrixlist1, matrixlist2,
                                    timespan, switchtime) {
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(#length(matrixlist1) == length(switchtime) ||
    length(switchtime) == 1)
  stopifnot(#length(matrixlist1) == length(timespan) ||
    length(timespan) == 1)
  stopifnot(isTRUE(all.equal(lapply(matrixlist1, dim),
                             lapply(matrixlist2, dim))))

  matrix1 <- Matrix::bdiag(matrixlist1)
  matrix2 <- Matrix::bdiag(matrixlist2)

  interpolateMatrices(matrix1, matrix2,
                      timespan = timespan, switchtime = switchtime)
}
