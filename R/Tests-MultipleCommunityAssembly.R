print("Script stops early if a test fails.")

numBasal <- 3; numConsum <- 2; numEnviron <- 2
print(paste("Settings", paste(numBasal, numConsum, numEnviron)))

print("Environment Seed Check")
egPool <- LawMorton1996_species(numBasal, numConsum)
CreateEnvironmentInteractions(egPool, numEnviron, LawMorton1996_CommunityMat,
                              Parameters = c(0.01, 10, 0.5, 0.2, 100, 10)
                              ) -> egInteractions
CreateEnvironmentInteractions(egPool, numEnviron, LawMorton1996_CommunityMat,
                              Parameters = c(0.01, 10, 0.5, 0.2, 100, 10),
                              EnvironmentSeeds = egInteractions$Seeds
                              ) -> egInteractions2
stopifnot(identical(egInteractions, egInteractions2))

egMatrix <- Matrix::bdiag(egInteractions$Mats)

print("History Seed Check")
CreateAssemblySequence(numBasal + numConsum, numEnviron) -> egEvents
CreateAssemblySequence(numBasal + numConsum, numEnviron,
                       HistorySeed = egEvents$Seed) -> egEvents2
stopifnot(identical(egEvents, egEvents2))

print("Dynamics check")
egDynamics <- PerCapitaDynamics_Type1(egPool$ReproductionRate, egMatrix,
                                      NumEnvironments = numEnviron)

egDynamics(0, rep(0, (numBasal + numConsum) * numEnviron), 0)

print("Events Check")
EliminationAndNeutralEvents(
  EventsAndSeed = egEvents, Species = nrow(egPool),
  PerCapitaDynamics = egDynamics,
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4
) -> egEventFunction

egCompareEventA <- which.max(egEvents$Events$Type == "Arrival")
egCompareEventE <- which.max(egEvents$Events$Type == "Extinct")

# Check below threshold removal.
print("Threshold Removal")
stopifnot(egEventFunction(0,
                          c(-1, 1E-7, 1, 1E3,
                            rep(0, (numBasal + numConsum) * numEnviron - 4)), 0)
           == c(0,0,1,1E3, rep(0, (numBasal + numConsum) * numEnviron - 4))
)

# Check that arrival works. Two modes: already present and new arrival.
print("Arrival, not present.")
with(egEvents$Events[egCompareEventA,], {
  expectedEvents <- egEventFunction(0,0,0, ReturnEvents = TRUE)
  expectedEvents[egCompareEventA, ]$Success <- TRUE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0.4
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(0,0,0, ReturnEvents = TRUE),
                      expectedEvents))
})

print("Arrival, already present.")
with(egEvents$Events[egCompareEventA,], {
  expectedEvents <- egEventFunction(0,0,0, ReturnEvents = TRUE)
  expectedEvents[egCompareEventA, ]$Success <- TRUE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  input[(Environment - 1) * nrow(egPool) + Species] <- 1E3
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0.4 +
    expected[(Environment - 1) * nrow(egPool) + Species]
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(0,0,0, ReturnEvents = TRUE),
                      expectedEvents))
})

# Check that elimination works.
print("Elimination, already present.")
with(egEvents$Events[egCompareEventE,], {
  expectedEvents <- egEventFunction(0,0,0, ReturnEvents = TRUE)
  expectedEvents[egCompareEventE, ]$Success <- TRUE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  input[(Environment - 1) * nrow(egPool) + Species] <- 1E3
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(0,0,0, ReturnEvents = TRUE),
                      expectedEvents))
})

print("Elimination, not present.")
with(egEvents$Events[egCompareEventE,], {
  expectedEvents <- egEventFunction(0,0,0, ReturnEvents = TRUE)
  expectedEvents[egCompareEventE, ]$Success <- FALSE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  input[(Environment - 1) * nrow(egPool) + Species] <- 0
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(0,0,0, ReturnEvents = TRUE),
                      expectedEvents))
})

print("Success.")
