# Islands and Interactions
# So for our initial interaction, we need to have a graph of patches.
# Each patch contains an ecosystem.
# Patches are effectively homogeneous at this point.
# Heterogeneity is hard to decide how to introduce, given that the obvious
# answer would be to modify the effective species carrying capacities, but
# those are implicit in the model, appearing due to reproduction rates and
# intra-species competition.
# Patches also need their dynamics under which the community is run.
# Then there are dispersal rates which are some function of species and distance
# from the patch perspective.
# From this perspective, islands are closely related sets of patches in that
# their adjacent distances are 0 (which does not mean instantaneous travel!).
#       Patch                        Patch
#   |-----------|                |-----------|
#   | Community |    Dispersal   | Community |
#   | Dynamics  | -------------> | Dynamics  |
#   |           | <------------- |           |
#   |           |                |           |
#   |-----------|                |-----------|
#
#
# Inputs: Species Pool, Species Interaction Matrix, which imply dynamics
# Species present on each patch, Patch distances

CsvRowSplit <- function(csv) {
  return(as.numeric(unlist(strsplit(csv, split = ", "))))
}

# https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-tp859876p859882.html
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

FindSteadyStateFromEstimate <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
  Tolerance = 1E-1,
  MaxAttempts = 1E2,
  maxRandVal = 1E6,
  Verbose = FALSE
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }

  if (is.character(Populations)) {
    init <- CsvRowSplit(Populations)
  } else {
    init <- Populations
  }
  init_min <- min(c(init, 0))
  init_max <- max(c(init, maxRandVal))

  anyZeroOrNotSame <- TRUE # Set T to make at least one look.
  epsilon <- Tolerance
  attempt <- 1

  # Run and check to see if anyone dies (bad) or changes (not at steady).
  while (anyZeroOrNotSame && MaxAttempts > attempt) {
    if (Verbose) print(paste0(attempt, ":"))
    init_old <- init

    init <- rootSolve::steady(
      y = init,
      func = Dynamics,
      parms = list(a = InteractionMatrix[com, com],
                   r = Pool$ReproductionRate[com],
                   epsilon = epsilon),
      positive = TRUE
    )$y
    if (Verbose) print(init)

    if (any(init < epsilon)) {
      # Someone died, reset to random location.
      if (Verbose) print("Died")
      anyZeroOrNotSame <- TRUE
      init[init < epsilon] <- runif(
        n = sum(init < epsilon),
        min = init_min,
        max = init_max
      )
    } else if (any(round(init / init_old, 1) != 1)) {
      # Not done, keep going.
      if (Verbose) print("Changed")
      anyZeroOrNotSame <- TRUE
    } else {
      anyZeroOrNotSame <- FALSE
    }

    attempt <- attempt + 1
  }

  if (attempt >= MaxAttempts) {
    warning(paste("Failed to converge after", attempt, "attempts."))
  }

  return(init)
}

FindSteadyStateFromBasal <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
  Tolerance = 1E-1,
  MaxAttempts = 1E2,
  maxRandVal = 1E6,
  Verbose = FALSE
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }

  if (is.character(Populations)) {
    init <- CsvRowSplit(Populations)
  } else {
    init <- Populations
  }
  init_min <- min(c(init, 0))
  init_max <- max(c(init, maxRandVal))

  basal <- (diag(InteractionMatrix) != 0)[com]

  epsilon <- Tolerance

  # Allow basal species to establish themselves.
  abund <- init; abund[!basal] <- 0
  abund <- rootSolve::steady(
    y = abund,
    func = Dynamics,
    parms = list(a = InteractionMatrix[com, com],
                 r = Pool$ReproductionRate[com],
                 epsilon = epsilon),
    positive = TRUE
  )$y
  init[basal] <- abund[basal]

  anyZeroOrNotSame <- TRUE # Set T to make at least one look.
  attempt <- 1

  # Run and check to see if anyone dies (bad) or changes (not at steady).
  while (anyZeroOrNotSame && MaxAttempts > attempt) {
    if (Verbose) print(paste0(attempt, ":"))
    init_old <- init

    init <- rootSolve::steady(
      y = init,
      func = Dynamics,
      parms = list(a = InteractionMatrix[com, com],
                   r = Pool$ReproductionRate[com],
                   epsilon = epsilon),
      positive = TRUE
    )$y
    if (Verbose) print(init)

    if (any(init < epsilon)) {
      # Someone died, reset to random location.
      if (Verbose) print("Died")
      anyZeroOrNotSame <- TRUE
      init[init < epsilon] <- runif(
        n = sum(init < epsilon),
        min = init_min,
        max = init_max
      )
    } else if (any(round(init / init_old, 1) != 1)) {
      # Not done, keep going.
      if (Verbose) print("Changed")
      anyZeroOrNotSame <- TRUE
    } else {
      anyZeroOrNotSame <- FALSE
    }

    attempt <- attempt + 1
  }

  if (attempt >= MaxAttempts) {
    warning(paste("Failed to converge after", attempt, "attempts."))
  }

  return(init)
}

FindSteadyStateFromSize <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
  Tolerance = 1E-1,
  MaxAttempts = 1E2,
  maxRandVal = 1E6,
  Verbose = FALSE
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }

  if (is.character(Populations)) {
    init <- CsvRowSplit(Populations)
  } else {
    init <- Populations
  }
  init_min <- min(c(init, 0))
  init_max <- max(c(init, maxRandVal))

  basal <- (diag(InteractionMatrix) != 0)[com]

  epsilon <- Tolerance

  # Allow basal species to establish themselves.
  abund <- init; abund[!basal] <- 0
  abund <- rootSolve::steady(
    y = abund,
    func = Dynamics,
    parms = list(a = InteractionMatrix[com, com],
                 r = Pool$ReproductionRate[com],
                 epsilon = epsilon),
    # positive = TRUE
    method = "runsteady"
  )$y
  if (Verbose) print(abund)

  # Order consumer's by size (correlated with trophic level).
  consumerOrder <- order(Pool[com,][!basal,]$Size)
  for (i in consumerOrder) {
    if (Verbose) print(paste("Consumer:", i))
    # Save,
    abund_old <- abund
    # Invade,
    abund[!basal][i] <- 1
    # Resolve, and
    abund <- tryCatch(quiet(rootSolve::steady(
      y = abund,
      func = Dynamics,
      parms = list(a = InteractionMatrix[com, com],
                   r = Pool$ReproductionRate[com],
                   epsilon = epsilon),
      method = "runsteady" # positive = TRUE
    )$y),
    error = function(e) {
      if (Verbose) warning(e)

      return(abund_old)
    })
    if (Verbose) print(abund)
    # Reject if something died.
    if (any(abund_old[abund_old > 0] < epsilon)) {
      abund <- abund_old
    }
  }
  # Use this as the new starting place.
  init[abund > 0] <- abund[abund > 0]
  if (Verbose) print(init)

  anyZeroOrNotSame <- TRUE # Set T to make at least one look.
  attempt <- 1

  # Run and check to see if anyone dies (bad) or changes (not at steady).
  while (anyZeroOrNotSame && MaxAttempts > attempt) {
    if (Verbose) print(paste0(attempt, ":"))
    init_old <- init

    init <- rootSolve::steady(
      y = init,
      func = Dynamics,
      parms = list(a = InteractionMatrix[com, com],
                   r = Pool$ReproductionRate[com],
                   epsilon = epsilon),
      positive = TRUE
    )$y
    if (Verbose) print(init)

    if (any(init < epsilon)) {
      # Someone died, reset to random location.
      if (Verbose) print("Died")
      anyZeroOrNotSame <- TRUE
      init[init < epsilon] <- runif(
        n = sum(init < epsilon),
        min = init_min,
        max = init_max
      )
    } else if (any(round(init / init_old, 1) != 1)) {
      # Not done, keep going.
      if (Verbose) print("Changed")
      anyZeroOrNotSame <- TRUE
    } else {
      anyZeroOrNotSame <- FALSE
    }

    attempt <- attempt + 1
  }

  if (attempt >= MaxAttempts) {
    warning(paste("Failed to converge after", attempt, "attempts."))
  }

  return(init)
}

Productivity <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }

  if (is.character(Populations)) {
    pop <- CsvRowSplit(Populations)
  } else {
    pop <- Populations
  }

  if (length(com) == 0) {
    return(NA)
  }

  comMatPos <- InteractionMatrix[com, com]; comMatPos[comMatPos < 0] <- 0
  poolRepPos <- Pool$ReproductionRate[com]; poolRepPos[poolRepPos < 0] <- 0

  parameters <- list(
    a = comMatPos,
    r = poolRepPos
  )

  return(
    sum(Dynamics(0, pop, parameters)[[1]] * pop / sum(pop))
  )
}

IslandLotkaVolterra <- function(t, y, parms) {
  with(as.list(parms), {
    list(as.numeric(y * (r + a %*% y) + d %*% y))
    # as.numeric since the solver doesn't know what Matrix::Matrices are.
  })
}

IslandDynamics <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  DispersalPool, # Species related dispersal rates
                 # Should have length == nrow(Pool). Multiplied by entries of
  DispersalIsland, # Island related dispersal rates. Is a matrix, row = to.
  Dynamics = IslandLotkaVolterra,
  Tolerance = 1E-1,
  Times = seq(from = 0,
              to = 2E4,
              by = 5E2),
  #DispersalRates, # For now, a matrix, to-from notation, for each island.
  #(NOT species on each island, since that requires the user knowing things in advance...)
  #DispersalRates, # List of matrices: column = species, row = (TO other) island, entry = travel rate
  Method = NULL, # Passed through to deSolve.
  TimesEvents = NULL
) {
  preprocessed <- IslandPreprocess(
    Pool = Pool,
    InteractionMatrix = InteractionMatrix,
    Communities = Communities,
    Populations = Populations,
    DispersalPool = DispersalPool,
    DispersalIsland = DispersalIsland,
    Tolerance = Tolerance
  )

  abundance <- with(
    preprocessed,
    deSolve::ode(
      abundance_init,
      times = Times,
      func = Dynamics,
      parms = parameters,
      events = list(func = function(t, y, parms) {
        y[y < parms$epsilon] <- 0
        y
      }, time = if(is.null(TimesEvents)) Times else TimesEvents),
      method = Method
    )
  )

  return(abundance)
}

IslandDynamicsSTODE <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  DispersalPool, # Species related dispersal rates
  # Should have length == nrow(Pool). Multiplied by entries of
  DispersalIsland, # Island related dispersal rates. Is a matrix, row = to.
  Dynamics = IslandLotkaVolterra,
  Tolerance = 1E-1,
  Times = seq(from = 0,
              to = 2E4,
              by = 5E2)
  #DispersalRates, # For now, a matrix, to-from notation, for each island.
  #(NOT species on each island, since that requires the user knowing things in advance...)
  #DispersalRates, # List of matrices: column = species, row = (TO other) island, entry = travel rate
) {
  preprocessed <- IslandPreprocess(
    Pool = Pool,
    InteractionMatrix = InteractionMatrix,
    Communities = Communities,
    Populations = Populations,
    DispersalPool = DispersalPool,
    DispersalIsland = DispersalIsland,
    Tolerance = Tolerance
  )

  # Allocate Time
  abundance <- matrix(with(preprocessed, c(0, abundance_init)), nrow = 1)

  times <- Times
  if (times[1] != 0) times <- c(0, times)

  for (i in 1:length(times)) {
    if (times[i] == 0) next()

    stodeRun <- with(
      preprocessed,
      rootSolve::stode(
        abundance[i - 1, -1], # remove time, running to i'th entry.
        times = times[i] - times[i - 1],
        func = Dynamics,
        parms = parameters,
        positive = TRUE
      )
    )

    abundance <- rbind(
      abundance,
      c(times[i], stodeRun$y)
      )

    if (stodeRun$steady) break
  }

  return(abundance)
}

IslandNumericalAssembly <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  # DispersalPool, # Species related dispersal rates
  # Should have length == nrow(Pool). Multiplied by entries of
  DispersalIsland, # Island related dispersal rates. Is a matrix, row = to.
  Dynamics = IslandLotkaVolterra,
  Tolerance = 1E-1,
  ArrivalDensity = 0.1,
  ArrivalEvents = 10,
  ArrivalSampler = c("rearrange", "iid"),
  IntegratorTimeStep = 100,
  ExtinctionTimeSteps = NULL,
  seed = NULL,
  ReturnValues = c("Abundance", "Sequence", "Pool", "Matrix",
                   "Uninvadable", "Steady", "Events")
) {
  if (!is.null(seed)) {
    if (exists(".Random.seed"))
      oldSeed <- .Random.seed
    set.seed(seed)
  }

  preprocessed <- IslandPreprocess(
    Pool = Pool,
    InteractionMatrix = InteractionMatrix,
    Communities = Communities,
    Populations = Populations,
    DispersalPool = 0,
    DispersalIsland = DispersalIsland,
    Tolerance = Tolerance
  )

  arrivalSamplerOptions <- c("rearrange", "iid")
  arrivalSampler <- match.arg(ArrivalSampler,
                              arrivalSamplerOptions)

  arrivalIDs <- # Goes down columns first.
    with(
      preprocessed,
      if (arrivalSampler == "rearrange") {
        lapply(
          DispersalIsland,
          function(i, len, events) {
            if (i > 0) {
              replicate(
                n = ceiling(events / len),
                sample.int(len, replace = FALSE)
              )[1:events]
            } else NULL
          },
          len = length(redCom),
          events = ArrivalEvents
        )
      } else if (arrivalSampler == "iid") {
        lapply(
          DispersalIsland,
          function(i, len, events) {
            if (i > 0) {
              sample.int(length(redCom),
                         size = ArrivalEvents,
                         replace = TRUE)
            } else NULL
          },
          len = length(redCom),
          events = ArrivalEvents
        )
      } else {
        stop(paste("ArrivalSampler not recognised. Should be one of",
                   arrivalSamplerOptions, collapse = ", "))
      }
    )


  # Setup
  abundance <- preprocessed$abundance_init
  speciesNum <- length(abundance) / nrow(DispersalIsland)

  allSpeciesOld <- lapply(1:nrow(DispersalIsland), function(i, s, ab, eps) {
    which(ab[((i - 1) * s + 1) : (i * s)] > eps)
  }, s = speciesNum, ab = abundance, eps = Tolerance)

  if ("Abundance" %in% ReturnValues) {
    abundanceHistory <- c("Time" = 0, abundance)
  }

  if ("Sequence" %in% ReturnValues) {
    #TODO DOUBLECHECK
    SequenceRetVal <- data.frame(
      Events = c(0, (IntegratorTimeStep) * 0:(ArrivalEvents - 1))
      # Addition = c(NA, ArrivalIDs),
      # Outcome = factor(NA, levels = c(
      #   "Present",
      #   "Type 1 (Failure)",
      #   "Type 2 (Invade)",
      #   "Type 3 (Contract)")),
      # Community = NA
    )

    # Links
    SRVLinks <- lapply(
      seq_along(DispersalIsland),
      function(i, nr, nc, invaders) {
        if (DispersalIsland[i] > 0) {
          retval <- data.frame(
            c(NA, invaders[[i]]),
            factor(NA, levels = c(
              "Present",
              "Type 1 (Failure)",
              "Type 2 (Invade)",
              "Type 3 (Contract)",
              "DNE")) # Better suggestion than DNE?
          )
          colnames(retval) <- paste(
            c("Addition", "Outcome"), nc[i], "->", nr[i]
          )
          retval
        }
      }, invaders = arrivalIDs,
      nr = row(DispersalIsland),
      nc = col(DispersalIsland))

    SRVLinks <- do.call(
      cbind,
      SRVLinks[!unlist(lapply(SRVLinks, is.null))]
    )

    SequenceRetVal <- cbind(
      SequenceRetVal,
      SRVLinks,
      # Islands
      do.call(
        cbind,
        lapply(
          1:nrow(DispersalIsland),
          function(i, events) {
            retval <- data.frame(
              rep(NA, events + 1)
            )
            colnames(retval) <- paste(
              "Community", i
            )
            retval
          }, events = ArrivalEvents)
      )
    )

    SequenceRetVal[1,
                   # Last columns correspond to Communities
                   (ncol(SequenceRetVal) - nrow(DispersalIsland) + 1):
                     ncol(SequenceRetVal)
    ] <- # Divide up all Species into islands and convert.
      unlist(lapply(allSpeciesOld, toString))
  }

  for (event in 1:ArrivalEvents) {
    abundanceChange <- rep(0, length(abundance))

    outcomes <- matrix(nrow = nrow(DispersalIsland),
                       ncol = ncol(DispersalIsland))

    # The elements in the copies tell us who is
    # attempting to invade.
    for (island in 1:nrow(DispersalIsland)) {# Rows are TO
      # For each island, check all invading links/nonzero col(DispersalIsland)
      # Retrieve the indices acting on and the species present.
      islandIndices <- ((island - 1)  * speciesNum + 1) : (island * speciesNum)
      islandIDs <- islandIndices[abundance[islandIndices] > Tolerance]
      islandSpecies <- allSpeciesOld[[island]] #((islandIDs - 1) %% speciesNum) + 1

      for (i in 1:ncol(DispersalIsland)) {# Columns are FROM
        # For each invading link,
        if (DispersalIsland[island, i] <= 0) next

        # entries go down each row in a column before moving to next column.
        candidateIndex <- (i - 1) * nrow(DispersalIsland) + island

        candidate <- arrivalIDs[[candidateIndex]][event]

        # check that the invasion is valid (present on the link) (DNE)

        if (abundance[candidate + (i - 1) * speciesNum] < ArrivalDensity) {
          outcomes[island, i] <- "DNE"
          next
        }

        # Check that additional species is not in islandIDs (Present)

        if (candidate %in% islandSpecies) {
          outcomes[island, i] <- "Present"
          next
        }

        # and check that the island is invadable (by that species). (Type 1)
        if (
          #length(islandSpecies) & # Is there anyone there?
          (preprocessed$redPool$ReproductionRate[candidate] +
           sum(preprocessed$redComMat[candidate,
                                      islandSpecies] *
               abundance[islandIDs],
               na.rm = TRUE) < 0) # Can you not reproduce from infinitesimal?
        ) {
          outcomes[island, i] <- "Type 1 (Failure)"
          next
        }

        # Remove the amount moving from its home
        abundanceChange[
          candidate + (i - 1) * speciesNum
        ] <- abundanceChange[
          candidate + (i - 1) * speciesNum
        ] - ArrivalDensity
        # Add the amount to its new home.
        abundanceChange[
          candidate + (island - 1)  * speciesNum
        ] <- abundanceChange[
          candidate + (island - 1)  * speciesNum
        ] + ArrivalDensity

        # print(paste(island, i, candidate))

        #TODO consider the rare case in which very little abundance remains,
        # and more population is selected to move than remains.
        # (Not too concerned: the abundance will be zeroed anyways if this is
        #  the case, and this would be very rare: (1/species)^2 at least for
        #  getting the right selections, but then also to have that low
        #  abundance?)

      }
    }

    #if (!all(abundanceChange == 0)) {
    # Apply the density record
    abundance <- abundance + abundanceChange

    # Run the system.
    abundance <- with(
      preprocessed,
      RMTRCode2::quiet(deSolve::ode(
        abundance,
        times = c(0, #IntegratorTimeStep/2,
                  IntegratorTimeStep),
        func = Dynamics,
        parms = parameters,
        events = list(func = function(t, y, parms) {
          y[y < parms$epsilon] <- 0
          y
        }, time = seq(0, to = IntegratorTimeStep,
                      by = if (!is.null(ExtinctionTimeSteps))
                        ExtinctionTimeSteps
                      else IntegratorTimeStep))
      ))
    )

    if ("Abundance" %in% ReturnValues) {
      abundanceHistory <- rbind(abundanceHistory,
                                abundance[c(1, nrow(abundance)), ])
    }

    abundance <- abundance[nrow(abundance), -1]

    # 0 anything that should be extinct.
    abundance[abundance < Tolerance] <- 0
    #}

    # Record results.

    allSpecies <- lapply(1:nrow(DispersalIsland), function(i, s, ab, eps) {
      which(ab[((i - 1) * s + 1) : (i * s)] > eps)
    }, s = speciesNum, ab = abundance, eps = Tolerance)

    if ("Sequence" %in% ReturnValues) {

      # Record Communities.
      SequenceRetVal[event + 1,
                     # Last columns correspond to Communities
                     (ncol(SequenceRetVal) - nrow(DispersalIsland) + 1):
                       ncol(SequenceRetVal)
      ] <- # Divide up all Species into islands and convert.
        unlist(lapply(allSpecies, toString))

      # Record Outcomes.
      for (clmn in seq(from = 2, by = 2, to = ncol(SequenceRetVal) - 3)) {
        # clmn is the Addition, clmn + 1 is the Outcome.
        # Source is the first number in the name.
        # Destination is the last number in the name.
        # Checks proceed in sequence.
        # Outcome is NA only if unevaluated.
        # Outcome is DNE if Addition is not present in source.
        # Outcome is Present if Addition is present in allSpeciesOld.
        # Outcome is Type 1 (Failure) if Addition could not invade.
        # Outcome is Type 2 (Invade) if Addition is present in allSpecies and
        #   all species in allSpeciesOld are present in allSpecies.
        # Outcome is Type 3 (Contract) if species in allSpeciesOld are not
        #   present in allSpecies.
        # The first 4 are the default or recorded in the outcomes string vector.
        # Only the last 2 need to be checked explicitly here.

        srce <- strsplit(colnames(SequenceRetVal[event, ])[clmn], split = " ")
        dest <- as.numeric(srce[[1]][4])
        srce <- as.numeric(srce[[1]][2])

        if (!is.na(outcomes[dest, srce])) {
          SequenceRetVal[event + 1, clmn + 1] <- outcomes[dest, srce]
        } else {
          # Check if allSpeciesOld are in allSpecies
          SequenceRetVal[event + 1, clmn + 1] <-
            if (length(allSpeciesOld[[dest]]) == 0) {
              if (SequenceRetVal[event + 1, clmn] %in% allSpecies[[dest]]) {
                "Type 2 (Invade)"
              } else {
                "Type 1 (Failure)"
              }
            }
          else if (all(allSpeciesOld[[dest]] %in% allSpecies[[dest]])) {
            if (SequenceRetVal[event + 1, clmn] %in% allSpecies[[dest]]) {
              "Type 2 (Invade)"
            } else {
              "Type 1 (Failure)"
            }
          } else {
            "Type 3 (Contract)"
          }
        }
      }
    }

    allSpeciesOld <- allSpecies

    # Check if uninvadable steady-state for each island relative to links.
    # If we are uninvadable AND (approximately) in steady-state, then we are done.
    # If so, we are done, otherwise continue running.

    completelyUninvadable <- TRUE
    for (island in 1:nrow(DispersalIsland)) {
      islandIndices <- ((island - 1)  * speciesNum + 1) : (island * speciesNum)
      islandIDs <- islandIndices[abundance[islandIndices] > Tolerance]
      islandSpecies <- allSpeciesOld[[island]]
      for (i in 1:ncol(DispersalIsland)) {
        # For each invading link,
        if (DispersalIsland[island, i] <= 0) next

        iSpecies <- allSpeciesOld[[i]]

        # Only concerned with species that are not already present.
        iSpecies <- iSpecies[!(iSpecies %in% islandSpecies)]

        # Check if no-one can invade from i to island.
        if (
          # length(islandSpecies) &
          any(preprocessed$redPool$ReproductionRate[iSpecies] +
              sum(preprocessed$redComMat[iSpecies,
                                         islandSpecies] *
                  abundance[islandIDs],
                  na.rm = TRUE) > 0)
        ) {
          completelyUninvadable <- FALSE
          break
        }
      }
    }

    if (completelyUninvadable) {
      # Check steady-state
      steady <- with(
        preprocessed,
        rootSolve::steady(
          abundance,
          func = Dynamics,
          parms = parameters
        ),
        time = c(0, 1), # If in steady-state, the duration shouldn't matter.
        method = "runsteady"
      )

      if (attr(steady, "steady") == TRUE){
        steadycheck <- all(# Double check that the abundances do match...
          (steady$y == 0 & abundance == 0) |
            (round((abundance - steady$y)/abundance,
                   digits = -log10(Tolerance)) == 0)
        )
        if (steadycheck) break
      }
    }
  }

  if (!is.null(seed)) {
    if (exists("oldSeed"))
      set.seed(oldSeed)
  }
  retval <- list()

  if ("Abundance" %in% ReturnValues) {
    abundanceHistory[, 1] <- cumsum(abundanceHistory[, 1])
    colnames(abundanceHistory)[-1] <- as.character(1:(ncol(abundanceHistory) - 1))
    retval$Abundance <- abundanceHistory
  }
  if ("Sequence" %in% ReturnValues) {
    retval$Sequence <- SequenceRetVal
  }
  if ("Pool" %in% ReturnValues) {
    retval$Pool <- Pool
  }
  if ("Matrix" %in% ReturnValues) {
    retval$Matrix <- InteractionMatrix
  }
  if ("Uninvadable" %in% ReturnValues) {
    retval$Uninvadable <- completelyUninvadable
  }
  if ("Steady" %in% ReturnValues) {
    retval$Steady <- steady
    retval$SteadyCheck <- steadycheck
  }
  if ("Events" %in% ReturnValues) {
    retval$Events <- event
  }

  return(retval)
}

IslandPreprocess <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  DispersalPool, # Species related dispersal rates
  # Should have length == nrow(Pool). Multiplied by entries of
  DispersalIsland,
  Tolerance
) {
  # Sanity check 1. ############################################################
  stopifnot(length(Communities) == length(Populations))

  if (length(DispersalPool) == 1) {
    DispersalPool <- rep(DispersalPool, nrow(Pool))
  }

  # Reduce to necessary information. ###########################################
  # Total list of species present.
  # Make sure that formatting is handled.
  CommunitiesNumeric <- lapply(Communities, function(com) {
    if (is.character(com)) {
      CsvRowSplit(com)
    } else if (is.numeric(com)) {
      com
    } else {
      stop("Community is not numeric or string.")
    }
  }
  )

  PopulationsNumeric <- lapply(Populations, function(pop) {
    if (is.character(pop)) {
      CsvRowSplit(pop)
    } else if (is.numeric(pop)) {
      pop
    } else {
      stop("Population is not numeric or string.")
    }
  }
  )

  # Sanity Check 2. ###########################################################
  if (!((all(
    unlist(lapply(CommunitiesNumeric, length)) ==
    unlist(lapply(PopulationsNumeric, length))
  )))) {
    stop("Entries in Communities and Populations differ.")
  }

  redCom <- sort(unique(unlist(CommunitiesNumeric))) # Uses recursive unlisting

  redPool <- Pool[redCom,]
  redComMat <- InteractionMatrix[redCom, redCom]

  redComs <- lapply(Communities, function(strCom) {
    which(redCom %in% CsvRowSplit(strCom))
  })

  redPops <- PopulationsNumeric # Unify language, all formatting done.

  # Sanity Check: Lengths should match.
  stopifnot(all(
    unlist(lapply(redComs, length)) == unlist(lapply(redPops, length))
  ))

  redDisPool <- DispersalPool[redCom]

  # Create the dispersal matrix. ###############################################

  # Alternative: create multiple sparse diagonal matrices, then r and c bind.
  # Loop through, keep track of which island we are on to add in an empty block.

  # Alternative for simple case of just an island, not species, travel matrix.
  # redDisMat <- DispersalRates
  # diag(redDisMat) <- -Matrix::colSums(dispersalMatrix)

  # So the question is how to format the diagonals then.
  dispersalDiags <- NULL
  for (i in (length(Communities) - 1):-(length(Communities) - 1)) {
    if (i == 0) next
    # Start in top right band, move to bottom left.
    dispersalDiags <- c(
      dispersalDiags,
      list(c(
        unlist(lapply(
          1:length(Communities),
          function(index, offset, mat, vec) {
            if (offset != 0 &&
                nrow(mat) >= index + offset &&
                index + offset >= 1)
              mat[index, index + offset] * vec
          },
          offset = i, mat = DispersalIsland, vec = redDisPool
        ))
      ))
    )
  }

  # Amount of travel from j to i is d[i,j]
  # Amount of gain to i from j is d[i,j]
  # Amount of travel from i is d[i,i]
  # This way, we can write the change in y from travel as d %*% y
  # Characterising this as proportions of a population, but that is assuming a
  # normalisation that I do not think is strictly necessary.
  # The matrix is sparse and has colsum = 0, diag < 0, offdiag >= 0.
  dispersalMatrix <- Matrix::bandSparse(
    n = length(redCom) * length(Communities),
    k = length(redCom) * c((length(Communities) - 1):1,
                           -(1:(length(Communities) - 1))),
    diagonals = dispersalDiags# c(
    # list(c(rep(0.0001, length(redCom)),   # Island 2 -> Island 1
    #        rep(0.0001, length(redCom)))), # Island 3 -> Island 2
    # list(rep(0, length(redCom))),         # Island 3 -> Island 1
    # list(c(rep(0.0001, length(redCom)),   # Island 1 -> Island 2
    #        rep(0.0001, length(redCom)))), # Island 2 -> Island 3
    # list(rep(0, length(redCom)))          # Island 1 -> Island 3
    # )
  )
  stopifnot(all(dispersalMatrix >= 0))
  diag(dispersalMatrix) <- -Matrix::colSums(dispersalMatrix)
  stopifnot(all(Matrix::colSums(dispersalMatrix) == 0))

  parameters <- list(
    r = rep(redPool$ReproductionRate, length(Communities)),
    a = Matrix::bdiag(rep(list(redComMat), length(Communities))),
    d = dispersalMatrix,
    epsilon = Tolerance
  )

  # Technically, a bit of extra work being done here since we already copied the
  # populations per island. The result is somewhat more readable though.
  abundance_init <- unlist(lapply(
    seq_along(redComs),
    function(i, com, pop, numPops) {
      # Interlace 0's with population values
      k <- 1
      retval <- rep(0, numPops)
      for (j in 1:numPops) {
        if (j %in% com[[i]]) {
          retval[j] <- pop[[i]][k]
          k <- k + 1
        }
      }
      return(retval)
    },
    com = redComs,
    pop = redPops,
    numPops = length(redCom)
  ))

  return(list(
    abundance_init = abundance_init,
    redCom = redCom,
    redComs = redComs,
    redPops = redPops,
    redPool = redPool,
    redDisPool = redDisPool,
    redComMat = redComMat,
    parameters = parameters,
    dispersalMatrix = dispersalMatrix
  ))
}
