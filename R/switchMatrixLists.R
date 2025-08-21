switchMatrixLists <- function(matrixlist1, matrixlist2, switchtimes) {
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(length(matrixlist1) == length(switchtimes) ||
              length(switchtimes) == 1)

  Matrix::bdiag(lapply(seq_along(matrixlist1), function(i, m1, m2, st) {
    switchMatrices(m1[[i]], m2[[i]], st[i])
  },
  m1 = matrixlist1, m2 = matrixlist2,
  st = if (length(switchtimes) == 1) {
    rep(switchtimes, length(matrixlist1))
  } else {switchtimes}
  ))
}
