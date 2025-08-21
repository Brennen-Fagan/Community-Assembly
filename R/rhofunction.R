# reproductive rate r' = r * rho ^ (sign(r)), where rho is a function of
# distance between patch and species in some sort of trait way.
rhofunction <- function(
  base = 2, offset = 0, multiplier = 1, metric = "euclidean"
) {
  force(base);force(offset);force(multiplier)
  function(m, n) {
    base ^ (offset - multiplier * dist(
      matrix(c(m, n), byrow = TRUE, nrow = 2), method = metric)
    )
  }
}

rho.noop <- function(m, n) {1}
rho.2.0.1.euclidean <- rhofunction()
rho.2.1.2.euclidean <- rhofunction(2, 1, 2)
rho.5.0.1.euclidean <- rhofunction(5, 0, 1)
rho.5.1.2.euclidean <- rhofunction(5, 1, 2)
rho.10.0.1.euclidean <- rhofunction(10, 0, 1)
rho.10.1.2.euclidean <- rhofunction(10, 1, 2)
