# Coupon Collector's Problem
defaultEvents <- function(
  NumberOfEnvironments, NumberOfSpecies, constant = 3
) {
  ceiling(
    NumberOfEnvironments * NumberOfSpecies * (
      log(NumberOfEnvironments * NumberOfSpecies) + constant
    )
  )
}
