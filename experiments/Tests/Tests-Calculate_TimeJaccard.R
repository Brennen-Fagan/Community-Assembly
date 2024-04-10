# Tests-Calculate_TimeJaccard.R

# Standard Opening: ###########################################################
library(RMTRCode2)
library(dplyr)

print("Script stops early if a test fails.")

# Trivial Tests: ##############################################################
# Stochastic Tests: ###########################################################
# Deterministic Tests: ########################################################
environments <- 4; nspecies <- c(2, 1)
configurations <- list(
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    0:5, 5:0, rep(0, 6)
  )),
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    5:0, rep(0, 6), 0:5 
  )),
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    c(1, 1, 2, 2, 3, 3), c(1, 1, 2, 2, 3, 3), c(1, 1, 2, 2, 3, 3)
  )),
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    c(0, 0, 1, 2, 3, 3), rep(0, 6), rep(0, 6)
  ))
)
times <- matrix(1:nrow(configurations[[1]]), ncol = 1)
configurations[[1]] <- cbind(times, configurations[[1]])
ecosystem <- list(
  Abundance = do.call(cbind, configurations),
  Parameters = list(
    EliminationThreshold = 0.1
  ),
  NumEnvironments = environments
)

result <- Calculate_TimeJaccard(
  loaded = ecosystem, 
  nspecies = nspecies, 
  minTime = 1
)

# Reminder: only based on Species Identity, not abundance.
#           Addtionally breaks up by basal (1:2) and consumer (3).
#           No species: NaN.

expected <- tibble::tibble(dplyr::bind_cols(
  expand.grid(
    Measurement = paste0("JaccardTemporal", c("", "_Basal", "_Consumer")),
    Time = 1:5, Environment = c(as.character(1:4), "Mean"),
    stringsAsFactors = FALSE
  )[, c(2, 3, 1)],
  Value = c(
    # Env. 1
    0.5, 0.5, NaN, # 1 New Basal, no consumers
    rep(c(0, 0, NaN), 3), # No changes,
    0.5, 0.5, NaN, # 1 Lost Basal, no consumers
    # Env. 2
    0.5, 0, 1, # No Change in Basal, 1 New Consumer
    rep(c(0, 0, 0), 3), # No changes,
    0.5, 1, 0, # 1 Lost Basal, No Change in Consumer
    # Env. 3
    rep(c(0, 0, 0), 5), # No Gains or Losses of Species
    # Env. 4
    NaN, NaN, NaN, # Nothing Present
    1, 1, NaN, # Gain a Basal
    rep(c(0, 0, NaN), 3),
    # Env. Mean (NaN not removed)
    NaN, NaN, NaN,
    (0+0+0+1)/4, (0+0+0+1)/4, NaN,
    0, 0, NaN,
    0, 0, NaN,
    (0.5+0.5+0+0)/4, (0.5+1+0+0)/4, NaN
    )
))

stopifnot(isTRUE(all.equal(expected, result)))

# Standard Closing: ###########################################################
print("Success.")
