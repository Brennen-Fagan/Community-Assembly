print("Script stops early if a test fails.")

library(RMTRCode2)

if (!exists("testseed")) {
  testseed <- runif(1) * 1E8
  if (exists(".Random.seed")) {
    old.seed <- .Random.seed
  }
  set.seed(testseed)
}

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

print("Test 0 behaviour of distance matrix.")
egDispersal0 <- CreateDispersalMatrix(
  Matrix::sparseMatrix(
    i = numEnviron, j = numEnviron, x = 0),
  SpeciesSpeeds = 1:(numBasal + numConsum)
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

print("Test less than fully connected dispersal.")
egDists <- matrix(1:numEnviron^2, nrow = numEnviron)
diag(egDists) <- 0
egDists[numEnviron - 1, numEnviron] <- Inf # Curiosity
egDists[1:2,4:5] <- 0
egDists[4:5,1:2] <- 0
egDispersal <- CreateDispersalMatrix(
  egDists, SpeciesSpeeds = 1:(numBasal + numConsum)
)

stopifnot(
  !any(egDispersal0 != 0),
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

print("Test reproducibility.")
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
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  Verbose = TRUE
)

stopifnot(isTRUE(all.equal(egResults_Dispersal$Abundance,
                           egResults_Dispersal2$Abundance)))

print("Test that no interaction is same as separate runs.")

egResults_Dispersal3 <- MultipleNumericalAssembly_Dispersal(
  Pool = egPool,
  NumEnvironments = numEnviron,
  InteractionMatrices = egInteractions,
  Events = egEvents,
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = egDispersal0,
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  Verbose = TRUE
)


print("Subtest 1: 1 Island at a time.")

egEventsSubsets <- lapply(
  1:numEnviron, function(i, events, seed) {
    temp <- events
    temp$Type[temp$Environment != i] <- "Dummy" # Should not run.
    temp$Environment <- 1
    list(Events = temp,
         Seed = seed)
  }, events = egEvents$Events, seed = egEvents$Seed
)
egEventsSubsets2 <- lapply(
  1:numEnviron, function(i, events, seed) {
    temp <- events
    temp$Type[temp$Environment != i] <- "Dummy" # Should not run.
    #temp$Environment <- 1
    list(Events = temp,
         Seed = seed)
  }, events = egEvents$Events, seed = egEvents$Seed
)
egMatsSubsets <- lapply(
  1:numEnviron, function(i, matrices, seed) {
    list(Mats = list(matrices[[i]]),
         Seeds = seed[i])
  }, matrices = egInteractions$Mats, seed = egInteractions$Seeds
)
egDynamicsSubsets <- lapply(
  egMatsSubsets, function(mat, reprate) {
    PerCapitaDynamics_Type1(reprate, mat$Mats[[1]],
                            NumEnvironments = 1)
  }, reprate = egPool$ReproductionRate
)

egResults_Dispersal4Calc <- lapply(
  1:numEnviron, function(i, pl, Di, mats, events, dynamics) {
    library(RMTRCode2)
    #print(i)
    #print("###################################################")
    MultipleNumericalAssembly_Dispersal(
      Pool = pl,
      NumEnvironments = 1,
      InteractionMatrices = mats,
      Events = events[[i]],
      PerCapitaDynamics = dynamics[[i]],
      DispersalMatrix = Di,
      EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
      Verbose = FALSE
    )
  },
  pl = egPool,
  Di = egDispersal0[1:5, 1:5],
  mats = egInteractions, #egMatsSubsets, # Bad implementation, we don't use the
                                         # full interaction matrix but for rate.
  events = egEventsSubsets, dynamics = egDynamicsSubsets
)

print("Subtest 2: All Islands, but only events on 1 at a time.")
egResults_Dispersal5Calc <- lapply(
  egEventsSubsets, function(
    events, pl, nE, In, Dy, Di) {
    library(RMTRCode2)
    MultipleNumericalAssembly_Dispersal(
      Pool = pl,
      NumEnvironments = nE,
      InteractionMatrices = In,
      Events = events,
      PerCapitaDynamics = Dy,
      DispersalMatrix = Di,
      EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
      Verbose = FALSE
    )
  },
  pl = egPool,
  nE = numEnviron,
  In = egInteractions,
  Dy = egDynamics,
  Di = egDispersal0
)

egResults_Dispersal4 <- list()
egResults_Dispersal4$Events <- do.call(
  rbind, lapply(1:length(egResults_Dispersal4Calc), function(i, dat) {
    temp <- dat[[i]]$Events[!is.na(dat[[i]]$Events$Success), ]
    temp$Environment <- i
    temp
  }, dat = egResults_Dispersal4Calc)
)
egResults_Dispersal4$Events <- egResults_Dispersal4$Events[
  order(egResults_Dispersal4$Events$Times),
]
egResults_Dispersal4$Abundance <- cbind(
  time = egResults_Dispersal4Calc[[1]]$Abundance[, 1],
  do.call(
    cbind,
    lapply(1:length(egResults_Dispersal4Calc), function(i, x) {
      x[[i]]$Abundance[
        , 1:(numBasal+numConsum) + 1
      ]
    }, x = egResults_Dispersal4Calc)
  )
)
egResults_Dispersal4$NumEnvironments <- length(egResults_Dispersal4Calc)
egResults_Dispersal4$EnvironmentSeeds <- unique(unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$EnvironmentSeeds
)))
egResults_Dispersal4$HistorySeed <- unique(unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$HistorySeed
)))
egResults_Dispersal4$Parameters <- unique(unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$Parameters
), recursive = FALSE))
names(egResults_Dispersal4$Parameters) <- names(egResults_Dispersal4Calc[[1]]$Parameters)
egResults_Dispersal4$Ellipsis <- unique(unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$Ellipsis
), recursive = FALSE))

stopifnot(
  unlist(lapply(
    1:length(egResults_Dispersal3),
    function(i) all.equal(egResults_Dispersal3[[i]],
                          egResults_Dispersal4[[i]],
                          use.names = FALSE, check.attributes = FALSE)
    )) == c(TRUE, "target is deSolve, current is matrix", rep(TRUE, 5))
)

egResults_Dispersal5 <- list()
egResults_Dispersal5$Events <- do.call(
  rbind, lapply(1:length(egResults_Dispersal5Calc), function(i, dat) {
    temp <- dat[[i]]$Events[!is.na(dat[[i]]$Events$Success), ]
    temp$Environment <- i
    temp
  }, dat = egResults_Dispersal5Calc)
)
egResults_Dispersal5$Events <- egResults_Dispersal5$Events[
  order(egResults_Dispersal5$Events$Times),
]
egResults_Dispersal5$Abundance <- cbind(
  time = egResults_Dispersal5Calc[[1]]$Abundance[, 1],
  do.call(
    cbind,
    lapply(1:length(egResults_Dispersal5Calc), function(i, x) {
      x[[i]]$Abundance[
        , (i - 1) * (numBasal + numConsum) + 1:(numBasal+numConsum) + 1
      ]
    }, x = egResults_Dispersal5Calc)
  )
)
egResults_Dispersal5$NumEnvironments <- length(egResults_Dispersal5Calc)
egResults_Dispersal5$EnvironmentSeeds <- unique(unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$EnvironmentSeeds
)))
egResults_Dispersal5$HistorySeed <- unique(unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$HistorySeed
)))
egResults_Dispersal5$Parameters <- unique(unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$Parameters
), recursive = FALSE))
names(egResults_Dispersal5$Parameters) <- names(egResults_Dispersal5Calc[[1]]$Parameters)
egResults_Dispersal5$Ellipsis <- unique(unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$Ellipsis
), recursive = FALSE))

stopifnot(
  unlist(lapply(
    1:length(egResults_Dispersal3),
    function(i) all.equal(egResults_Dispersal3[[i]],
                          egResults_Dispersal5[[i]],
                          use.names = FALSE, check.attributes = FALSE)
  )) == c(TRUE, "target is deSolve, current is matrix", rep(TRUE, 5))
)

# MultipleNumericalAssembly_Dispersal, Trophics ################################

egTrophicFunction <- CalculateTrophicStructure(
  Pool = egPool,
  NumEnvironments = numEnviron,
  InteractionMatrices = egInteractions,
  EliminationThreshold = 10^-4
)

egTrophicLevels <- apply(
  egResults_Dispersal2$Abundance[, -1], # No Time Column
  MARGIN = 1, # Rows
  FUN = egTrophicFunction
)

stopifnot(
  sapply(egTrophicLevels, function(lst) {
    (length(lst$EdgeVertexList) == numEnviron) && # Length is correct
      all(sapply(lst$EdgeVertexList, function(lst2) # Edgelist is correct
        (is.na(lst2$Edges) && is.na(lst2$Vertices)) || # No species case
          (all(lst2$Edges$effectNormalised >= 0) && # Normalisation
             all(lst2$Edges$effectNormalised <= 1) &&
             all(lst2$Edges %>% dplyr::group_by( # Check that norm'ed to 1.
               to, effectSign
             ) %>% dplyr::summarise(
               .groups = "drop",
               test = isTRUE(all.equal(sum(effectNormalised), 1))
               # Note identical returns FALSE.
             ) %>% dplyr::pull(test))
          )
      ))
  })
)

print("Success.")

if (exists("old.seed")) set.seed(old.seed)
