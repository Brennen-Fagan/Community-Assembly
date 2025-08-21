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
