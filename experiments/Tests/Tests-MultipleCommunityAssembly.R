print("Script stops early if a test fails.")

# A useful construction:
# for (i in 1:100) {
#   testseed <- runif(i)[i]*1E8; print(testseed)
#   source("Tests-MultipleCommunityAssembly.R")
# }

# E.g.
# invadeCounter <- NULL
# for (i in 1:100) {
#   testseed <- runif(i)[i]*1E8; print(testseed)
#   source("Tests-MultipleCommunityAssembly.R")
#   if(identical(invadeCounter, NULL))
#     invadeCounter <-
#       table(factor(Invadables$values, levels = c(FALSE, TRUE)),
#             factor(Invadables$simulation, levels = c(FALSE, TRUE))
#     )
#   else
#     invadeCounter <- invadeCounter +
#       table(factor(Invadables$values, levels = c(FALSE, TRUE)),
#             factor(Invadables$simulation, levels = c(FALSE, TRUE)))
# }

library(RMTRCode2)

if (!exists("testseed")) {
  testseed <- runif(1) * 1E8
  if (exists(".Random.seed")) {
    old.seed <- .Random.seed
  }
  set.seed(testseed)
  print(testseed)
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
stopifnot(identical(egEvents, egEvents2),
          nrow(egEvents$Events) == 2 * 10)

CreateAssemblySequence(
  Species = numBasal + numConsum, NumEnvironments = numEnviron,
  ArrivalRate = egArrivalRate,
  ArrivalFUN = ArrivalFUN_Example2,
  ExtinctFUN = ExtinctFUN_Example2
) -> egEvents3

CreateAssemblySequence(
  Species = numBasal + numConsum, NumEnvironments = numEnviron,
  ArrivalRate = egArrivalRate,
  ArrivalFUN = ArrivalFUN_Example2,
  ExtinctFUN = ExtinctFUN_Example2,
  HistorySeed = egEvents3$Seed
) -> egEvents4
stopifnot(identical(egEvents3, egEvents4),
          # Default argument divided by rate.
          max(egEvents3$Events$Times) < 10 / egArrivalRate)

stopifnot(rapply(object = egEvents, f = typeof) ==
            rapply(object = egEvents3, f = typeof))

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

egEventsNo1 <- egEvents
egEventsNo1$Events <- egEventsNo1$Events[
  !(egEventsNo1$Events$Species == 1 |
      egEventsNo1$Events$Environment == 1
  ),
]

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
  # ComputeInteractionMatrix = LawMorton1996_CommunityMat,
  CharacteristicRate = max(
    unlist(lapply(egInteractions$Mats, function(mat) {
      abs(eigen(mat, only.values = TRUE)$values)
    }))
  ),
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 10),
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = egDispersal,
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  ArrivalEvents = 10, ArrivalRate = egArrivalRate, ArrivalFUN = NULL,
  ExtinctEvents = 10, ExtinctRate = NULL, ExtinctFUN = NULL,
  # EnvironmentSeeds = egInteractions$Seeds,
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
    temp$Environment <- 1 # 1 # Each environment has different matrices!
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
    PerCapitaDynamics_Type1(reprate, mat$Mats[[1]], # only 1 mat in each subset.
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
  Di = egDispersal0[1:nrow(egPool), 1:nrow(egPool)],
  mats = egInteractions, #egMatsSubsets, # Bad implementation, we don't use the
  # full interaction matrix but for rate.
  events = egEventsSubsets, dynamics = egDynamicsSubsets
)

print("Subtest 2: All Islands, but only events on 1 at a time.")
egResults_Dispersal5Calc <- lapply(
  egEventsSubsets2, function(
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

print("Subtest 3: 1 Island at a time using characteristic rate.")
egResults_Dispersal6Calc <- lapply(
  1:numEnviron, function(i, pl, Di, mats, events, dynamics) {
    #library(RMTRCode2)
    MultipleNumericalAssembly_Dispersal(
      Pool = pl,
      NumEnvironments = 1,
      #InteractionMatrices = mats,
      CharacteristicRate = max(
        unlist(lapply(mats$Mats, function(mat) {
          abs(eigen(mat, only.values = TRUE)$values)
        }))
      ),
      Events = events[[i]],
      PerCapitaDynamics = dynamics[[i]],
      DispersalMatrix = Di,
      EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
      Verbose = FALSE
    )
  },
  pl = egPool,
  Di = egDispersal0[1:nrow(egPool), 1:nrow(egPool)],
  mats = egInteractions,
  events = egEventsSubsets, dynamics = egDynamicsSubsets
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
egResults_Dispersal4$ReactionTime <- unique(unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$ReactionTime
)))
egResults_Dispersal4$HistorySeed <- unique(unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$HistorySeed
)))
egResults_Dispersal4$Parameters <- unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$Parameters
), recursive = FALSE)
egResults_Dispersal4$Parameters <- egResults_Dispersal4$Parameters[
  unique(names(egResults_Dispersal4$Parameters))
  ]
egResults_Dispersal4$Ellipsis <- unique(unlist(lapply(
  egResults_Dispersal4Calc, function(x) x$Ellipsis
), recursive = FALSE))

stopifnot(
  unlist(lapply(
    names(egResults_Dispersal3)[-2], # Test abundances separately.
    function(i) all.equal(as.matrix(egResults_Dispersal3[[i]])[,],
                          as.matrix(egResults_Dispersal4[[i]])[,],
                          use.names = FALSE, check.attributes = FALSE)
  ))  == "TRUE"
)

stopifnot(
  max(ifelse(egResults_Dispersal3$Abundance <= egResults_Dispersal3$Parameters$EliminationThreshold &
           egResults_Dispersal4$Abundance <= egResults_Dispersal4$Parameters$EliminationThreshold,
         0,
         abs(egResults_Dispersal3$Abundance - egResults_Dispersal4$Abundance) / egResults_Dispersal4$Abundance
  )) < 0.01 # Do we break 1% deviation?
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
egResults_Dispersal5$ReactionTime <- unique(unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$ReactionTime
)))
egResults_Dispersal5$HistorySeed <- unique(unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$HistorySeed
)))
egResults_Dispersal5$Parameters <- unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$Parameters
), recursive = FALSE)
egResults_Dispersal5$Parameters <- egResults_Dispersal5$Parameters[
  unique(names(egResults_Dispersal5$Parameters))
  ]
egResults_Dispersal5$Ellipsis <- unique(unlist(lapply(
  egResults_Dispersal5Calc, function(x) x$Ellipsis
), recursive = FALSE))

stopifnot(
  unlist(lapply(
    names(egResults_Dispersal4),
    function(i) all.equal(as.matrix(egResults_Dispersal4[[i]])[,],
                          as.matrix(egResults_Dispersal5[[i]])[,],
                          use.names = FALSE, check.attributes = FALSE)
  ))  == "TRUE"
)

egResults_Dispersal6 <- list()
egResults_Dispersal6$Events <- do.call(
  rbind, lapply(1:length(egResults_Dispersal6Calc), function(i, dat) {
    temp <- dat[[i]]$Events[!is.na(dat[[i]]$Events$Success), ]
    temp$Environment <- i
    temp
  }, dat = egResults_Dispersal6Calc)
)
egResults_Dispersal6$Events <- egResults_Dispersal6$Events[
  order(egResults_Dispersal6$Events$Times),
]
egResults_Dispersal6$Abundance <- cbind(
  time = egResults_Dispersal6Calc[[1]]$Abundance[, 1],
  do.call(
    cbind,
    lapply(1:length(egResults_Dispersal6Calc), function(i, x) {
      x[[i]]$Abundance[
        , 1:(numBasal+numConsum) + 1
      ]
    }, x = egResults_Dispersal6Calc)
  )
)
egResults_Dispersal6$NumEnvironments <- length(egResults_Dispersal6Calc)
egResults_Dispersal6$ReactionTime <- unique(unlist(lapply(
  egResults_Dispersal6Calc, function(x) x$ReactionTime
)))
egResults_Dispersal6$HistorySeed <- unique(unlist(lapply(
  egResults_Dispersal6Calc, function(x) x$HistorySeed
)))
egResults_Dispersal6$Parameters <- unlist(lapply(
  egResults_Dispersal6Calc, function(x) x$Parameters
), recursive = FALSE)
egResults_Dispersal6$Parameters <- egResults_Dispersal6$Parameters[
  unique(names(egResults_Dispersal6$Parameters))
  ]
egResults_Dispersal6$Ellipsis <- unique(unlist(lapply(
  egResults_Dispersal6Calc, function(x) x$Ellipsis
), recursive = FALSE))

stopifnot(
  unlist(lapply(
    names(egResults_Dispersal4),
    function(i) if (i == "Parameters") { # We didn't include seeds in this case
      all.equal(as.matrix(egResults_Dispersal4[[i]][1:5])[,],
                as.matrix(egResults_Dispersal6[[i]])[,],
                use.names = FALSE, check.attributes = FALSE)
    } else {
      all.equal(as.matrix(egResults_Dispersal4[[i]])[,],
                as.matrix(egResults_Dispersal6[[i]])[,],
                use.names = FALSE, check.attributes = FALSE)
    }
  ))  == "TRUE"
)

print("Subtest 4: No dispersal if no connections")
egResults_Dispersal7 <- MultipleNumericalAssembly_Dispersal(
  Pool = egPool,
  NumEnvironments = numEnviron,
  InteractionMatrices = egInteractions,
  Events = egEventsNo1,
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = CreateDispersalMatrix(
    matrix(0, nrow = numEnviron, ncol = numEnviron),
    SpeciesSpeeds = 1:(numBasal + numConsum)
  ),
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  Verbose = FALSE
)

stopifnot(
  !any(egResults_Dispersal7$Abundance[, 1 + 1:(numBasal + numConsum)] > 0),
  !any(egResults_Dispersal7$Abundance[
    , 1 + 1 + (0:(numEnviron - 1)) * (numBasal + numConsum)] > 0)
)


print("Subtest 5: Interrupted Evaluation")
TimeNoLate <- 40
egEventsNo1NoLate <- egEventsNo1
egEventsNo1NoLate$Events <- egEventsNo1NoLate$Events[
  egEventsNo1NoLate$Events$Times < TimeNoLate,
]
# NOTE: Need to include the last event of the previous, since we use that as the
# (post-)ending state.
egEventsNo1Late <- egEventsNo1
egEventsNo1Late$Events <- egEventsNo1Late$Events[
  egEventsNo1Late$Events$Times >= max(egEventsNo1NoLate$Events$Times),
]

egResults_Dispersal8 <- MultipleNumericalAssembly_Dispersal(
  Pool = egPool,
  NumEnvironments = numEnviron,
  InteractionMatrices = egInteractions,
  Events = egEventsNo1NoLate,
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = CreateDispersalMatrix(
    matrix(0, nrow = numEnviron, ncol = numEnviron),
    SpeciesSpeeds = 1:(numBasal + numConsum)
  ),
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  Verbose = FALSE
)

TimeLastEvent <- max(egResults_Dispersal8$Events$Times)
TimeNoLateRow <- which.max(egResults_Dispersal8$Abundance[, 1] >= TimeLastEvent) - 1

# Together, should be same as egResults_Dispersal7
egResults_Dispersal9 <- MultipleNumericalAssembly_Dispersal(
  PopulationInitial = egResults_Dispersal8$Abundance[TimeNoLateRow, -1],
  TimeInitial = egResults_Dispersal8$Abundance[TimeNoLateRow, 1],
  Pool = egPool,
  NumEnvironments = numEnviron,
  InteractionMatrices = egInteractions,
  Events = egEventsNo1Late,
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = CreateDispersalMatrix(
    matrix(0, nrow = numEnviron, ncol = numEnviron),
    SpeciesSpeeds = 1:(numBasal + numConsum)
  ),
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  Verbose = FALSE
)

# Should be the same as above. Diff is egEventsNo1 vs egEventsNo1 without prior
# events in the events listing.
egResults_Dispersal10 <- MultipleNumericalAssembly_Dispersal(
  PopulationInitial = egResults_Dispersal8$Abundance[TimeNoLateRow, -1],
  TimeInitial = egResults_Dispersal8$Abundance[TimeNoLateRow, 1],
  Pool = egPool,
  NumEnvironments = numEnviron,
  InteractionMatrices = egInteractions,
  Events = egEventsNo1,
  PerCapitaDynamics = egDynamics,
  DispersalMatrix = CreateDispersalMatrix(
    matrix(0, nrow = numEnviron, ncol = numEnviron),
    SpeciesSpeeds = 1:(numBasal + numConsum)
  ),
  EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
  Verbose = FALSE
)

stopifnot(
  isTRUE(all.equal(egResults_Dispersal9, egResults_Dispersal10)),
  isTRUE(all.equal(
    rbind(egResults_Dispersal8$Events,
          # Remove the duplicated event.
          egResults_Dispersal9$Events[-1, ]),
    egResults_Dispersal7$Events
  )),
  # More complicated: Abundances except for the gap?
  # There is some carry through. How about a tolerance of 1E-4 for now.
  # We'll run this a few times to see if that tolerance is exceeded, and,
  # if it is, see if anything changes about the resulting story.
  unlist(lapply(
    1:ncol(egResults_Dispersal7$Abundance),
    function(col) {
      cbind({
        tempTimeOfLastEvent <- tail(egResults_Dispersal8$Events$Times, 1)
        tempRowOfLastEvent <- which(
          egResults_Dispersal8$Abundance[, 1] == tempTimeOfLastEvent
        )
        tempRowOfFirstEvent <- which(
          egResults_Dispersal9$Abundance[, 1] ==
            egResults_Dispersal9$Events$Times[1]
        )
        rbind(
          egResults_Dispersal8$Abundance[1:tempRowOfLastEvent, ],
          egResults_Dispersal9$Abundance[
            tempRowOfFirstEvent:nrow(egResults_Dispersal9$Abundance),]
        )
      }[, col], {
        tempTimeOfLastEvent <- tail(egResults_Dispersal8$Events$Times, 1)
        tempRowOfLastEvent <- which(
          egResults_Dispersal7$Abundance[, 1] == tempTimeOfLastEvent
        )
        tempRowOfFirstEvent <- which(
          egResults_Dispersal7$Abundance[, 1] ==
            egResults_Dispersal9$Events$Times[1]
        )
        rbind(
          egResults_Dispersal7$Abundance[1:tempRowOfLastEvent, ],
          egResults_Dispersal7$Abundance[
            tempRowOfFirstEvent:nrow(egResults_Dispersal7$Abundance),]
        )
      }[, col]) -> temp

      if(isTRUE(all.equal(temp[, 1], temp[, 2]))) return(TRUE)

      # If the difference is small or they are both to be eliminated, TRUE.
      if(all(abs((temp[, 1] - temp[, 2]) / temp[, 1]) < 1E-4 |
             (temp[, 1] < egResults_Dispersal7$Parameters$EliminationThreshold &
              temp[, 2] < egResults_Dispersal7$Parameters$EliminationThreshold)))
        return(TRUE)
      else
        return(FALSE)
    }))
)

# Check Invadability: ##########################################################
# Recall:
# egResults_Dispersal2 <- MultipleNumericalAssembly_Dispersal(
#   Pool = egPool,
#   NumEnvironments = numEnviron,
#   InteractionMatrices = egInteractions,
#   Events = egEvents,
#   PerCapitaDynamics = egDynamics,
#   DispersalMatrix = egDispersal,
#   EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
#   Verbose = TRUE
# )
print("Invasibility")

# Calculate:
Invadables <- CalculateLocalInvadables_BareBones(
  with(egResults_Dispersal2,
       Abundance[nrow(Abundance), -1]),
  egDynamics,
  Environments = numEnviron,
  ArrivalDensity = 0.4,
  DispersalMatrix = egDispersal
)

Invadables2 <- CalculateLocalInvadables_KnockOn(
  with(egResults_Dispersal2,
       Abundance[nrow(Abundance), -1]),
  egDynamics,
  Environments = numEnviron,
  ArrivalDensity = 0.4,
  DispersalMatrix = egDispersal,
  TimeScale = 1/max(
    unlist(lapply(egInteractions$Mats, function(mat) {
      abs(eigen(mat, only.values = TRUE)$values)
    }))
  )
)

Invadables3 <- CalculateLocalInvadables_KnockOnRunSteady(
  with(egResults_Dispersal2,
       Abundance[nrow(Abundance), -1]),
  egDynamics,
  Environments = numEnviron,
  ArrivalDensity = 0.4,
  DispersalMatrix = egDispersal
)

Invadables <- data.frame(
  names = as.numeric(names(Invadables)),
  values1 = Invadables,
  values2 = Invadables2$invadable,
  values3 = Invadables3$invadable,
  simulation = NA
)

Invadables$species <- ((Invadables$names - 1) %% nrow(egPool) + 1)
Invadables$environment <- ((Invadables$names - 1) %/% nrow(egPool)) + 1

# Verify:
for (i in 1:nrow(Invadables)) {
  testEvent <- list(Events = data.frame(
    Times = with(egResults_Dispersal2,
                 Abundance[nrow(Abundance), 1]) + 0.01,
    Species = Invadables$species[i],
    Environment = Invadables$environment[i],
    Type = "Arrival",
    Success = NA
  ), Seed = NA)

  test <- MultipleNumericalAssembly_Dispersal(
    Pool = egPool,
    NumEnvironments = numEnviron,
    InteractionMatrices = egInteractions,
    Events = testEvent,
    PerCapitaDynamics = egDynamics,
    DispersalMatrix = egDispersal,
    EliminationThreshold = 10^-4, ArrivalDensity = 0.4,
    Verbose = FALSE,
    PopulationInitial = with(egResults_Dispersal2,
                             Abundance[nrow(Abundance), -1]),
    TimeInitial = with(egResults_Dispersal2,
                       Abundance[nrow(Abundance), 1])
  )

  # TRUE if species establishes and grows
  testAbundance <- tail(test$Abundance[
    , 1 + (Invadables$environment[i] - 1) * numEnviron + Invadables$species[i]
  ])

  Invadables$simulation[i] <-
    test$Events$Success && testAbundance[6] > 0
}

# Fundamentally, we can't really predict invasion success, because
# it isn't terribly well-defined here.
# What is an appropriate time scale?
# What about if it dies out but irreparably damages the system?
# Instead of stopping if the long term dynamics "simulation" disagree
# with our invadability metric "values", we merely compare and contrast.
# Ideally, they shouldn't diverge much: if a species is identified as
# uninvadable, we should be reasonably confident it won't stay in the system.
# In the case where it manages to succeed despite appearing to be about to
# drop out, we will end up miss-classifying it as a rescue event.
print(Invadables)

print("Success.")
print("Testseed:")
print(testseed)

if (exists("old.seed")) set.seed(old.seed)
