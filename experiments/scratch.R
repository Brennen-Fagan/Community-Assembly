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
    #}

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

      allSpeciesOld <- allSpecies
    }

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
          func = IslandLotkaVolterra,
          parms = parameters
        ),
        time = c(0, 1), # If in steady-state, the duration shouldn't matter.
        method = "runsteady"
      )

      if (attr(steady, "steady") == TRUE &&
          all(# Double check that the abundances do match...
            (steady$y == 0 & abundance == 0) |
            (round((abundance - steady$y)/abundance,
                   digits = -log10(Tolerance)) == 0)
          )){
        break
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
  }
  if ("Events" %in% ReturnValues) {
    retval$Events <- event
  }

  return(retval)
}

# debugonce(IslandNumericalAssembly)

dmat <- matrix(c( # row = to, column = from
  0, 1, 0, # Island 2 -> 1, Try 0, 0, 0 if need to test asymmetry.
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
  ReturnValues = c("Sequence", "Abundance",
                   "Uninvadable", "Steady", "Events"),
  seed = 1
)) -> temp

tempabund <- lapply(
  list(temp$Abundance),
  function(listtib, islandsNum) {
    speciesNum <- (ncol(listtib) - 1) / islandsNum

    retval <- tidyr::pivot_longer(
      listtib %>% as.data.frame,
      cols = !"Time",
      names_to = "Species",
      values_to = "Abundance"
    ) %>% dplyr::mutate(
      Species = as.numeric(Species),
      Island = floor((Species - 1) / speciesNum) + 1,
      Species = ((Species - 1) %% speciesNum) + 1 # maps 1,2,3,4 %% 4 -> 1:4
    ) %>% dplyr::group_by(
      Species, Island
    ) %>% dplyr::mutate(
      Native = dplyr::first(Abundance > 0),
      Invasive = !Native
    ) %>% dplyr::ungroup() %>% dplyr::group_by(
      Species, Time
    ) %>% dplyr::mutate(
      IslandsOccupied = sum(Abundance > 0) / 2
    ) %>% dplyr::ungroup() %>% dplyr::group_by(
      Species, Island
    ) %>% dplyr::mutate(
      Endemic = Native & dplyr::first(IslandsOccupied) == 1,
      Type = dplyr::case_when(
        Endemic ~ "Endemic",
        Native ~ "Native",
        Invasive ~ "Invasive",
        TRUE ~ "Oops"
      )
    )

    return(retval)
  }, islandsNum = 3
)

tempplot <- lapply(
  tempabund,
  function(tib) {
    ggplot2::ggplot(
      tib,
      ggplot2::aes(
        x = Time,
        y = Abundance,
        color = as.factor(Species),
        # linetype = Type,
        group = interaction(Species, Island)
      )
    ) + ggplot2::geom_point(
    ) + ggplot2::facet_grid(
      Island ~ Type
    ) + ggplot2::scale_y_log10(
      limits = c(10^-12, 10^4.5)
    )
  }
)
