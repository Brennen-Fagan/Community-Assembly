library(microbenchmark)
library(profvis)
library(Matrix)

RMTRCode2::PerCapitaDynamics_Mutualistic2(
  Pool$ReproductionRate,
  Matrix::bdiag(InteractionMatrices$Mats),
  NumEnvironments = 10,
  SpeciesTypes = Species,
  Saturations = Saturation
) -> dynamics

PerCapitaDynamics_Mutualistic3(
  Pool$ReproductionRate, 
  Matrix::bdiag(InteractionMatrices$Mats), 
  NumEnvironments = 10, 
  SpeciesTypes = Species,
  Saturations = Saturation
) -> dynamics3



dynamics_test <- function() {
  outval <- matrix(NA,nrow = ncol(result$Abundance) - 1, ncol = 10)
  index <- 1
  for (i in seq(1, nrow(result$Abundance), length.out = 10)) {
    target <- result$Abundance[i, ]
    outval[, index] <- dynamics(target[1], target[-1])[, 1]
    index <- index + 1
  }
  return(outval)
}
dynamics3_test <- function() {
  outval <- matrix(NA,nrow = ncol(result$Abundance) - 1, ncol = 10)
  index <- 1
  for (i in seq(1, nrow(result$Abundance), length.out = 10)) {
    target <- result$Abundance[i, ]
    outval[, index] <- dynamics3(target[1], target[-1])[, 1]
    index <- index + 1
  }
  return(outval)
}

dynamics_test_result <- dynamics_test()
dynamics3_test_result <- dynamics3_test()

stopifnot(all(dynamics_test_result == dynamics3_test_result))

microbenchmark(dynamics_test(), dynamics3_test(), times = 10)
profvis(dynamics3_test())
