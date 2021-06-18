IslandNumericalAssembly <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  DispersalPool, # Species related dispersal rates
  # Should have length == nrow(Pool). Multiplied by entries of
  DispersalIsland, # Island related dispersal rates. Is a matrix, row = to.
  Dynamics = GeneralisedLotkaVolterra,
  Tolerance = 1E-1,
  ArrivalDensity = 0.1,
  ArrivalEvents = 10,
  ArrivalSampler = c("rearrange", "iid"),
  IntegratorTimeStep = 100,
  ExtinctionTimeSteps = NULL,
  seed = NULL,
  ReturnValues = c("Abundance", "Sequence", "Pool", "Matrix")
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

  arrivalIDs <-
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
            c("Addition", "Outcome"), nr[i], nc[i]
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

    allSpeciesOld <- lapply(1:nrow(DispersalIsland), function(i, s, ab, eps) {
      which(ab[((i - 1) * s + 1) : (i * s)] > eps)
    }, s = speciesNum, ab = abundance, eps = Tolerance)

    SequenceRetVal[1,
                   # Last columns correspond to Communities
                   (ncol(SequenceRetVal) - nrow(DispersalIsland) + 1):
                     ncol(SequenceRetVal)
    ] <- # Divide up all Species into islands and convert.
      unlist(lapply(allSpeciesOld, toString))
  }

  for (event in 1:ArrivalEvents) {
    abundanceChange <- rep(0, length(abundance))

    outcomes <- NULL

    # The elements in the copies tell us who is
    # attempting to invade.
    for (island in 1:nrow(DispersalIsland)) {
      # For each island, check all invading links/nonzero col(DispersalIsland)
      # Retrieve the indices acting on and the species present.
      islandIndices <- ((island - 1)  * speciesNum + 1) : (island * speciesNum)
      islandIDs <- which(abundance[islandIndices] > Tolerance)
      islandSpecies <- allSpeciesOld[[island]] #((islandIDs - 1) %% speciesNum) + 1

      for (i in 1:ncol(DispersalIsland)) {
        # For each invading link,
        if (DispersalIsland[island, i] <= 0) next

        # entries go down each row in a column before moving to next column.
        candidateIndex <- (i - 1) * nrow(DispersalIsland) + island

        candidate <- arrivalIDs[[candidateIndex]][event]

        # check that the invasion is valid (present on the link) (DNE)

        if (abundance[candidate + (i - 1) * speciesNum] < ArrivalDensity) {
          outcomes <- c(outcomes, "DNE")
          next
        }

        # Check that additional species is not in islandIDs (Present)

        if (candidate %in% islandSpecies) {
          outcomes <- c(outcomes, "Present")
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
          outcomes <- c(outcomes, "Type 1 (Failure)")
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

        outcomes <- c(outcomes, NA)
      }
    }

    if (!all(abundanceChange == 0)) {
      # Apply the density record
      abundance <- abundance + abundanceChange

      # Run the system.
      abundance <- with(
        preprocessed,
        deSolve::ode(
          abundance,
          times = c(0, #IntegratorTimeStep/2,
                    IntegratorTimeStep),
          func = IslandLotkaVolterra,
          parms = parameters,
          events = list(func = function(t, y, parms) {
            y[y < parms$epsilon] <- 0
            y
          }, time = seq(0, to = IntegratorTimeStep,
                        by = if (!is.null(ExtinctionTimeSteps))
                          ExtinctionTimeSteps
                        else IntegratorTimeStep))
        )
      )

      if ("Abundance" %in% ReturnValues) {
        abundanceHistory <- rbind(abundanceHistory,
                                  abundance[c(1, nrow(abundance)), ])
      }

      abundance <- abundance[nrow(abundance), -1]

      # 0 anything that should be extinct.
      abundance[abundance < Tolerance] <- 0


    }

    # Record results.

    if ("Sequence" %in% ReturnValues) {
      allSpecies <- lapply(1:nrow(DispersalIsland), function(i, s, ab, eps) {
        which(ab[((i - 1) * s + 1) : (i * s)] > eps)
      }, s = speciesNum, ab = abundance, eps = Tolerance)

      # Record Communities.
      SequenceRetVal[event + 1,
                     # Last columns correspond to Communities
                     (ncol(SequenceRetVal) - nrow(DispersalIsland) + 1):
                       ncol(SequenceRetVal)
                     ] <- # Divide up all Species into islands and convert.
        unlist(lapply(allSpecies, toString))

      # Record Outcomes.
      outcomeNum <- 1
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

        if (!is.na(outcomes[outcomeNum])) {
          SequenceRetVal[event + 1, clmn + 1] <- outcomes[outcomeNum]
        } else {
          srce <- strsplit(colnames(SequenceRetVal[event, ])[clmn], split = " ")
          dest <- as.numeric(srce[[1]][3])
          srce <- as.numeric(srce[[1]][2])

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
        outcomeNum <- outcomeNum + 1
      }

      allSpeciesOld <- allSpecies
    }

    # Check if uninvadable steady-state for each island relative to links.
    # If so, we are done, otherwise continue running.

    completelyUninvadable <- TRUE
    for (island in 1:nrow(DispersalIsland)) {
      islandIndices <- ((island - 1)  * speciesNum + 1) : (island * speciesNum)
      islandIDs <- which(abundance[islandIndices] > Tolerance)
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
      break
    }
  }

  if (!is.null(seed)) {
    if (exists("oldSeed"))
      set.seed(oldSeed)
  }
  retval <- list()

  if ("Abundance" %in% ReturnValues) {
    abundanceHistory[, 1] <- cumsum(abundanceHistory[, 1])
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

  return(retval)
}

# debugonce(IslandNumericalAssembly)

dmat <- matrix(c(
  0, 1, 0, # Island 2 -> 1
  1, 0, 1, # Island 1 -> 2, Island 3 -> 2
  0, 1, 0  # Island 2 -> 3
), nrow = 3, ncol = 3, byrow = TRUE)

with(communitiesEX[[5]],
IslandNumericalAssembly(
  Pool = pools[[DatasetID[1]]][[CombnNum[1]]],
  InteractionMatrix = mats[[DatasetID[1]]][[CombnNum[1]]],
  DispersalPool = 0.001, DispersalIsland = dmat,
  Communities = c(
    list(Communities[1]),
    rep("", nrow(dmat) - 2),
    Communities[2]
  ),
  Populations = c(
    list(CommunityAbund[1]),
    rep("", nrow(dmat) - 2),
    list(CommunityAbund[2])
  ),
  ArrivalEvents = 100,
  ReturnValues = c("Sequence", "Abundance")
))
