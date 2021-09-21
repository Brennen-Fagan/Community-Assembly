print("Script stops early if a test fails.")

numBasal <- 3; numConsum <- 2; numEnviron <- 5
print(paste("Settings", paste(numBasal, numConsum, numEnviron)))
egArrivalRate <- 0.1

# CreateEnvironmentInteractions ################################################
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

# CreateAssemblySequence #######################################################
print("History Seed Check")
CreateAssemblySequence(numBasal + numConsum, numEnviron,
                       ArrivalRate = egArrivalRate) -> egEvents
CreateAssemblySequence(numBasal + numConsum, numEnviron,
                       ArrivalRate = egArrivalRate,
                       HistorySeed = egEvents$Seed) -> egEvents2
stopifnot(identical(egEvents, egEvents2))

# PerCapitaDynamics_Type1 ######################################################
print("Dynamics check")
egDynamics <- PerCapitaDynamics_Type1(egPool$ReproductionRate, egMatrix,
                                      NumEnvironments = numEnviron)

egDynamics(0, rep(0, (numBasal + numConsum) * numEnviron), 0)

# EliminationAndNeutralEvents ##################################################
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
  expectedEvents <- egEventFunction(ReturnEvents = TRUE)
  expectedEvents[egCompareEventA, ]$Success <- TRUE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0.4
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(ReturnEvents = TRUE),
                      expectedEvents))
})

print("Arrival, already present.")
with(egEvents$Events[egCompareEventA,], {
  expectedEvents <- egEventFunction(ReturnEvents = TRUE)
  expectedEvents[egCompareEventA, ]$Success <- TRUE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  input[(Environment - 1) * nrow(egPool) + Species] <- 1E3
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0.4 +
    expected[(Environment - 1) * nrow(egPool) + Species]
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(ReturnEvents = TRUE),
                      expectedEvents))
})

# Check that elimination works.
print("Elimination, already present.")
with(egEvents$Events[egCompareEventE,], {
  expectedEvents <- egEventFunction(ReturnEvents = TRUE)
  expectedEvents[egCompareEventE, ]$Success <- TRUE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  input[(Environment - 1) * nrow(egPool) + Species] <- 1E3
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(ReturnEvents = TRUE),
                      expectedEvents))
})

print("Elimination, not present.")
with(egEvents$Events[egCompareEventE,], {
  expectedEvents <- egEventFunction(ReturnEvents = TRUE)
  expectedEvents[egCompareEventE, ]$Success <- FALSE
  input <- rep(0, (numBasal + numConsum) * numEnviron)
  input[(Environment - 1) * nrow(egPool) + Species] <- 0
  expected <- input
  expected[(Environment - 1) * nrow(egPool) + Species] <- 0
  stopifnot(egEventFunction(Times, input, 0) == expected)
  stopifnot(identical(egEventFunction(ReturnEvents = TRUE),
                      expectedEvents))
})

# CreateDispersalMatrix ########################################################
# Test for failures.
print("Test Malformed Distance Matrix.")
Failure <- tryCatch({
  egDists <- matrix(1:6, nrow = 3)
  diag(egDists) <- 0
  egDispersal <- CreateDispersalMatrix(
    egDists, SpeciesSpeeds = 1:(numBasal + numConsum)
  )}, error = function(e) {
    TRUE
  })
stopifnot(Failure == TRUE)

print("Test Self-looping Distance Matrix.")
Failure <- tryCatch({
  egDists <- matrix(1:9, nrow = 3)
  egDispersal <- CreateDispersalMatrix(
    egDists, SpeciesSpeeds = 1:(numBasal + numConsum)
  )}, error = function(e) {
    TRUE
  })

print("Test regular behaviour of distance matrix.")
egDists <- matrix(1:numEnviron^2, nrow = numEnviron)
diag(egDists) <- 0
egDists[numEnviron - 1, numEnviron] <- Inf # Curiosity
egDispersal <- CreateDispersalMatrix(
  egDists, SpeciesSpeeds = 1:(numBasal + numConsum)
)

stopifnot(
  nrow(egDispersal) == numEnviron * (numBasal + numConsum),
  ncol(egDispersal) == nrow(egDispersal),
  isTRUE(all.equal(Matrix::colSums(egDispersal), rep(0, ncol(egDispersal)))),
  Matrix::diag(egDispersal) <= 0,
  # Check travel from last environment to second last is 0.
  isTRUE(Matrix::all.equal(egDispersal[
    (numEnviron - 2) * (numBasal + numConsum) + 1:(numBasal + numConsum),
    (numEnviron - 1) * (numBasal + numConsum) + 1:(numBasal + numConsum)
    ], as(as(Matrix::Diagonal(n = (numBasal + numConsum),
                              x = rep(0, (numBasal + numConsum))),
             "CsparseMatrix"), "dgCMatrix"),
    check.attributes = FALSE
    )),
  # Check only diagonal of the opposite block is non-zero
  Matrix::diag(egDispersal[
    (numEnviron - 1) * (numBasal + numConsum) + 1:(numBasal + numConsum),
    (numEnviron - 2) * (numBasal + numConsum) + 1:(numBasal + numConsum)
  ]) > 0,
  egDispersal[
    (numEnviron - 1) * (numBasal + numConsum) + 1:(numBasal + numConsum),
    (numEnviron - 2) * (numBasal + numConsum) + 1:(numBasal + numConsum)
  ][upper.tri(egDispersal[
    (numEnviron - 1) * (numBasal + numConsum) + 1:(numBasal + numConsum),
    (numEnviron - 2) * (numBasal + numConsum) + 1:(numBasal + numConsum)
  ])] == 0,
  egDispersal[
    (numEnviron - 1) * (numBasal + numConsum) + 1:(numBasal + numConsum),
    (numEnviron - 2) * (numBasal + numConsum) + 1:(numBasal + numConsum)
  ][lower.tri(egDispersal[
    (numEnviron - 1) * (numBasal + numConsum) + 1:(numBasal + numConsum),
    (numEnviron - 2) * (numBasal + numConsum) + 1:(numBasal + numConsum)
  ])] == 0
  )

# MultipleNumericalAssembly_Dispersal ##########################################
print("Testing Main Function: Dispersal")

egResults_Dispersal <- MultipleNumericalAssembly_Dispersal(
  Pool = egPool,
  NumEnvironments = numEnviron,
  ComputeInteractionMatrix = LawMorton1996_CommunityMat,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 10),
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = egDispersal,
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  ArrivalEvents = 10, ArrivalRate = egArrivalRate, ArrivalFUN = NULL,
  ExtinctEvents = 10, ExtinctRate = NULL, ExtinctFUN = NULL,
  EnvironmentSeeds = egInteractions$Seeds,
  HistorySeed = egEvents$Seed
)

egResults_Dispersal2 <- MultipleNumericalAssembly_Dispersal(
  Pool = egPool,
  NumEnvironments = numEnviron,
  InteractionMatrices = egInteractions,
  Events = egEvents,
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = egDispersal,
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4
)

print("Success.")
