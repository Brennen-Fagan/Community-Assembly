### Interventions: ############################################################
actTrivially <- function(elem) {force(elem); function(...){elem}}

switchMatrices <- function(matrix1, matrix2, switchtime) {
  stopifnot(dim(matrix1) == dim(matrix2))
  force(matrix1);force(matrix2);force(switchtime)
  function(t, ...) {
    if (t < switchtime) {
      matrix1
    } else if (t > switchtime) {
      matrix2
    } else {
      (matrix1 + matrix2)/2
    }
  }
}

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

# Choose a sigmoidal interpolation between the matrices.
# Problem: we want no deviation at 0 or t = timespan, but not too much in mid.
# Crit. Value = first numb. in tanh. 4 => 2% dev at edge, ~75% done in mid. 50%.
# Require: Output is a function.
# Require: Output: t = timespan / 2 => out = (mat1 + mat2) / 2.
interpolateMatrices <- function(matrix1, matrix2, timespan, switchtime = 0) {
  stopifnot(dim(matrix1) == dim(matrix2))
  force(matrix1);force(matrix2);force(timespan);force(switchtime)
  function(t, ...) {
    # if (t < switchtime) {
    # matrix1
    # } else {
    matrix1 + (matrix2 - matrix1) * (
      tanh(4 * ( (t - switchtime) / timespan - 0.5)) + 1
    ) / 2
    # }
  }
}

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