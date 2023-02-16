print("Script stops early if a test fails.")

library(RMTRCode2)

if (!exists("testseed")) {
  testseed <- runif(1) * 1E8
  if (exists(".Random.seed")) {
    old.seed <- .Random.seed
  }
  set.seed(testseed)
  print(testseed)
}

stopifnot(
  # Test if is symmetric, as expected, that there are no 5's (off) or 6's (on).
  isSymmetric(
    temp <- GeneralSystems_CategoricalMatrix(
      100, c(rep(1, 4), 0, rep(1,5)), rep(1, 5), upper.tri.only = FALSE
    )
  ),
  !any(diag(temp) == 6),
  !any(temp[upper.tri(temp)] == 5 || temp[lower.tri(temp)] == 5),
  # Test if identity replacement is fine.
  isSymmetric(
    tempout <- GeneralSystems_CategoricalToProportioned(
      temp,
      list(
        function(n) {cbind(rep(1, n), rep(1, n))},
        function(n) {cbind(rep(2, n), rep(2, n))},
        function(n) {cbind(rep(3, n), rep(3, n))},
        function(n) {cbind(rep(4, n), rep(4, n))},
        function(n) {cbind(rep(5, n), rep(5, n))},
        function(n) {cbind(rep(6, n), rep(6, n))},
        function(n) {cbind(rep(7, n), rep(7, n))},
        function(n) {cbind(rep(8, n), rep(8, n))},
        function(n) {cbind(rep(9, n), rep(9, n))},
        function(n) {cbind(rep(10, n), rep(10, n))}
      ),
      list(
        function(n) {rep(1, n)},
        function(n) {rep(2, n)},
        function(n) {rep(3, n)},
        function(n) {rep(4, n)},
        function(n) {rep(5, n)},
        function(n) {rep(6, n)}
      ))
  ),
  all(tempout[upper.tri(tempout)] == temp[upper.tri(temp)])
)
