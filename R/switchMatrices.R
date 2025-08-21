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
