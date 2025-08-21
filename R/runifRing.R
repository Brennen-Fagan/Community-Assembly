# Run runif and organise in a smooth-ish ring.
runifRing <- function(n, ...) {
  indices <- if (n %% 2) {
    # Odd (why?)
    c(1, seq(from = 2, by = 2, to = n), seq(from = n, by = -2, to = 2))
  } else {
    # Even.
    c(1, seq(from = 2, by = 2, to = n), seq(from = n - 1, by = -2, to = 2))
  }
  sort(runif(n, ...))[indices]
}
