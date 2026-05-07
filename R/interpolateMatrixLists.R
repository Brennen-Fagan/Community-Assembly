interpolateMatrixLists <- function(
  matrixlist1, matrixlist2, timespans, switchtimes
) {
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(length(matrixlist1) == length(timespans) ||
              length(timespans) == 1)
  stopifnot(length(matrixlist1) == length(switchtimes) ||
              length(switchtimes) == 1)

  Matrix::bdiag(lapply(
    seq_along(matrixlist1), function(i, m1, m2, ts, st) {
      interpolateMatrices(m1[[i]], m2[[i]], ts[i], st[i])
    },
    m1 = matrixlist1, m2 = matrixlist2,
    ts = if (length(timespans) == 1) {
      rep(timespans, length(matrixlist1))
    } else {timespans},
    st = if (length(switchtimes) == 1) {
      rep(switchtimes, length(matrixlist1))
    } else {switchtimes}
  ))
}
