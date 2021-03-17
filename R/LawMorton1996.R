# Law and Morton, 1996 #########################################################
LawMorton1996_NumericalAssembly <- function(
  Basal = NULL,
  Consumer = NULL,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  LogBodySize = c(-2, -1, -1, 0), # c(-2, -1) for Basal, c(-1, 0) for Consumer
  EliminationThreshold = 10^-4,
  IntegratorTimeStep = 1000,
  ArrivalDensity = 0.1,
  ArrivalEvents = 10,
  ArrivalSampler = c("rearrange", "iid"),
  InnerTimeStepSize = 100,
  ReturnValues = c("Abundance", "Sequence", "Pool", "Matrix"),
  Pool = NULL,
  CommunityMat = NULL,
  seed = NULL
) {
  # EliminationThreshold = (note: actual threshold is X * Threshold.)

  if (is.null(Pool) & !is.null(CommunityMat)) {
    stop("CommunityMat should not be specified if Pool is not specified.")
  }

  if (is.null(Basal) & is.null(Consumer) & is.null(Pool)) {
    stop("Either number of Basal and Consumer species must be specified or a Pool must be provided.")
  }

  if (is.null(Pool)) {
    # Create species pool.
    Pool <- LawMorton1996_species(Basal, Consumer, Parameters, LogBodySize)
    speciesNum <- Basal + Consumer
  } else {
    speciesNum <- nrow(Pool)
  }

  if (is.null(CommunityMat)) {
    # Create interaction matrix.
    CommunityMat <- LawMorton1996_CommunityMat(Pool, Parameters)
    # Constructed so that i is row, j is column.
  }

  if (!is.null(seed)) {
    oldSeed <- .Random.seed
    set.seed(seed)
  }

  # Setup.
  CurrentAbundance <- rep(0, speciesNum)
  SpeciesPresent <- NULL
  if (!is.null(InnerTimeStepSize)) {
    stopifnot(InnerTimeStepSize < IntegratorTimeStep)
    # Store Time and Abundances for plotting.
    TimeAbundances <- rbind(
      matrix(
        data = 0, nrow = 1, ncol = 1 + length(CurrentAbundance)
      ),
      matrix(
        nrow = (IntegratorTimeStep / InnerTimeStepSize + 1) * ArrivalEvents,
        ncol = 1 + length(CurrentAbundance)
      )
    )
    TimeAbundances_row <- 2
    TimeAbundances_time <- 0
  }

  # Generate arrival events.
  arrivalSampler <- match.arg(ArrivalSampler, c("rearrange", "iid"))
  if (arrivalSampler == "iid") {
    ArrivalIDs <- sample.int(nrow(Pool), size = ArrivalEvents, replace = TRUE)
  } else if (arrivalSampler == "rearrange") {
    ArrivalIDs <- replicate(
      n = ceiling(ArrivalEvents / nrow(Pool)),
      sample.int(nrow(Pool), replace = FALSE)
    )[
      1:ArrivalEvents
    ]
  }


  if ("Sequence" %in% ReturnValues) {
    SequenceRetVal <- data.frame(
      Events = c(0, (IntegratorTimeStep + 1) * 0:(ArrivalEvents - 1)),
      Addition = c(NA, ArrivalIDs),
      Outcome = factor(NA, levels = c(
        "Present",
        "Type 1 (Failure)",
        "Type 2 (Invade)",
        "Type 3 (Contract)")),
      Community = NA
    )
  }

  eventNumber <- 1

  # Resolve each arrival.
  for (ID in ArrivalIDs) {
    SpeciesPresent_Old <- SpeciesPresent
    CurrentAbundance[ID] <- CurrentAbundance[ID] + ArrivalDensity

    if (!(ID %in% SpeciesPresent)) {
      SpeciesPresent[length(SpeciesPresent) + 1] <- ID
    }

    # Run_GLV <- LawMorton1996_NumIntegration(
    #   CommunityMat[SpeciesPresent, SpeciesPresent],
    #   Pool$ReproductionRate[SpeciesPresent],
    #   CurrentAbundance[SpeciesPresent],
    #   OuterTimeStepSize = IntegratorTimeStep,
    #   InnerTimeStepSize = if(is.null(InnerTimeStepSize)) {
    #     IntegratorTimeStep
    #   } else {InnerTimeStepSize}
    # )

    rootSolveCounter <- 0
    Run_GLV <- NA
    while(rootSolveCounter < 5 && any(is.na(Run_GLV))) {
      Run_GLV <- tryCatch(
        LawMorton1996_NumIntegration(
          CommunityMat[SpeciesPresent, SpeciesPresent],
          Pool$ReproductionRate[SpeciesPresent],
          CurrentAbundance[SpeciesPresent] + abs(rnorm(n = length(SpeciesPresent))) * sqrt(rootSolveCounter),
          OuterTimeStepSize = IntegratorTimeStep,
          InnerTimeStepSize = if(is.null(InnerTimeStepSize)) {
            IntegratorTimeStep
          } else {InnerTimeStepSize}
        ),
        error = function(e) {
          return(NA)
        })
      rootSolveCounter <- rootSolveCounter + 1
    }
    stopifnot(!is.na(Run_GLV))

    CurrentAbundance_New <- Run_GLV[nrow(Run_GLV), 2:(ncol(Run_GLV))]

    CurrentAbundance_New <- ifelse(
      CurrentAbundance_New > EliminationThreshold * CurrentAbundance[SpeciesPresent] &
        CurrentAbundance_New > 0,
      CurrentAbundance_New,
      0
    )

    # Record Run_GLV if we need to before we update the SpeciesPresent.
    if (!is.null(InnerTimeStepSize)) {
      #TODO Consider if this matrix might benefit from being sparse or similar.
      TimeAbundances[
        TimeAbundances_row:(TimeAbundances_row + nrow(Run_GLV) - 1),
        SpeciesPresent + 1
      ] <- Run_GLV[, 2:ncol(Run_GLV)]

      TimeAbundances[
        TimeAbundances_row:(TimeAbundances_row + nrow(Run_GLV) - 1), 1
      ] <- Run_GLV[, 1] + TimeAbundances_time

      TimeAbundances_row <- TimeAbundances_row + nrow(Run_GLV)
      TimeAbundances_time <- TimeAbundances_time + Run_GLV[nrow(Run_GLV), 1] + 1
    }

    # We ignore the new addition if the new addition died out.
    if (CurrentAbundance_New[length(CurrentAbundance_New)] == 0 &&
        length(CurrentAbundance_New) > 1) {
      atSteadyState <- all(
        round(CurrentAbundance_New[-length(CurrentAbundance_New)] /
                CurrentAbundance[SpeciesPresent][-length(CurrentAbundance_New)],
              -log10(EliminationThreshold)) == 1
      )
    } else {
      atSteadyState <- all(
        round(CurrentAbundance_New / CurrentAbundance[SpeciesPresent],
              -log10(EliminationThreshold)) == 1
      )
    }

    CurrentAbundance[SpeciesPresent] <- CurrentAbundance_New

    SpeciesPresent <- which(CurrentAbundance > 0)
    eventNumber <- eventNumber + 1

    if ("Sequence" %in% ReturnValues) {
      if (all(SpeciesPresent_Old %in% SpeciesPresent) &&
          ID %in% SpeciesPresent_Old) {
        SequenceRetVal$Outcome[eventNumber] <- "Present"
      } else if (all(SpeciesPresent_Old %in% SpeciesPresent) &&
                 !(ID %in% SpeciesPresent)
      ) {
        SequenceRetVal$Outcome[eventNumber] <- "Type 1 (Failure)"
      } else if (all(SpeciesPresent_Old %in% SpeciesPresent) &&
                 (ID %in% SpeciesPresent)
      ) {
        SequenceRetVal$Outcome[eventNumber] <- "Type 2 (Invade)"
      } else {
        SequenceRetVal$Outcome[eventNumber] <- "Type 3 (Contract)"
      }
      SequenceRetVal$Community[eventNumber] <- I(list(SpeciesPresent))
    }

    # We can add a check when ``sufficiently many'' species have been added.
    # Or, if it is cheap, we can just check immediately if a community is invadable.
    # The check only works if the system is essentially at steady state.
    if (eventNumber > nrow(Pool) &&
        atSteadyState &&
        length(SpeciesPresent) > 1 &&
        LawMorton1996_CheckUninvadable(
          AbundanceRow = c(NA, CurrentAbundance),
          Pool = Pool, CommunityMatrix = CommunityMat
        )) {
      break()
    }
  }

  retval <- list()
  if ("Abundance" %in% ReturnValues) {
    if (!is.null(InnerTimeStepSize)) {
      lastRow <- sum(apply(TimeAbundances, MARGIN = 1, FUN = function(x) {
        !all(is.na(x[-1]))
      }))
      retval$Abundance <- TimeAbundances[1:lastRow, ]
    } else {
      retval$Abundance <- c((IntegratorTimeStep + 1) * ArrivalEvents - 1,
                            CurrentAbundance)
    }
  }
  if ("Sequence" %in% ReturnValues) {
    retval$Sequence <- SequenceRetVal
  }
  if ("Pool" %in% ReturnValues) {
    retval$Pool <- Pool
  }
  if ("Matrix" %in% ReturnValues) {
    retval$Matrix <- CommunityMat
  }

  if (!is.null(seed)) {
    set.seed(oldSeed)
  }

  return(retval)
}

LawMorton1996_PermanenceAssembly <- function(
  Basal = NULL,
  Consumer = NULL,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  LogBodySize = c(-2, -1, -1, 0), # c(-2, -1) for Basal, c(-1, 0) for Consumer
  ArrivalEvents = 10,
  ArrivalSampler = c("rearrange", "iid"),
  ReturnValues = c("Community", "Equilibrium",
                   "Sequence", "States",
                   "Pool", "Matrix"),
  Pool = NULL,
  CommunityMat = NULL,
  seed = NULL
) {
  if (is.null(Pool) & !is.null(CommunityMat)) {
    stop("CommunityMat should not be specified if Pool is not specified.")
  }

  if (is.null(Basal) & is.null(Consumer) & is.null(Pool)) {
    stop("Either number of Basal and Consumer species must be specified or a Pool must be provided.")
  }

  if (is.null(Pool)) {
    # Create species pool.
    Pool <- LawMorton1996_species(Basal, Consumer, Parameters, LogBodySize)
    speciesNum <- Basal + Consumer
  } else {
    speciesNum <- nrow(Pool)
  }

  if (is.null(CommunityMat)) {
    # Create interaction matrix.
    CommunityMat <- LawMorton1996_CommunityMat(Pool, Parameters)
    # Constructed so that i is row, j is column.
  }

  if (!is.null(seed)) {
    oldSeed <- .Random.seed
    set.seed(seed)
  }

  # Generate arrival events.
  arrivalSampler <- match.arg(ArrivalSampler, c("rearrange", "iid"))
  if (arrivalSampler == "iid") {
    ArrivalIDs <- sample.int(nrow(Pool), size = ArrivalEvents, replace = TRUE)
  } else if (arrivalSampler == "rearrange") {
    ArrivalIDs <- replicate(
      n = ceiling(ArrivalEvents / nrow(Pool)),
      sample.int(nrow(Pool), replace = FALSE)
    )[
      1:ArrivalEvents
    ]
  }

  community <- NULL
  if ("Sequence" %in% ReturnValues) {
    records <- data.frame(
      Events = 0:ArrivalEvents,
      Addition = c(NA, ArrivalIDs),
      Outcome = factor(NA, levels = c(
        "Present",
        "Type 1 (Failure)",
        "Type 2 (Permanent)",
        "Type 3 (!Permanent)")),
      Community = NA
    )
  }
  eventNumber <- 2

  #TODO This really should be some form of object
  #     In that case, we would have a checkPermanence routine.
  #     The routine would do what we do on lines 258 - 388 that
  #     we reproduce effectively starting on line 431.
  #     (We won't just recursively calculate permanence because
  #     we may not need it if we are in a permanent state already.)
  statesEncountered <- list(
    IDs = list(toString(integer(0))),
    Equilibria = list(NULL), # If does not exist, use NA.
    Permanent = TRUE
    )
  abundance <- 0

  for (ID in ArrivalIDs) {
    if (!(ID %in% community)) {
      # Check if species can increase when rare. Yes -> Add to community.
      if (Pool$ReproductionRate[ID] +
          sum(CommunityMat[ID, community] * abundance, na.rm = TRUE) > 0
          ) { # (Note: 0x0 mat * 0 is 0.)
        # Add species to the community. Note that we sort for uniqueness.
        sorted <- sort(c(community, ID), index.return = TRUE)
        community <- sorted$x
        if (length(community) == 1) {
          abundance <- 0.1
        } else {
          abundance <- c(abundance, 0.1)[sorted$ix] # 0.1 is initial density.
        }

        # Store properties.
        stateNumber <- which(toString(community) == statesEncountered$IDs)
        if (length(stateNumber) == 0) {
          stateNumber <- length(statesEncountered$IDs) + 1
          statesEncountered$IDs[stateNumber] <- toString(community)

          # Best way to avoid ending up at the other boundaries?
          # In practice, the positive (pos) argument is actually non-negative.
          # What happens when we have a community without a reassembly path?
          # Easiest way is to start from the previous equilibrium and run.
          # Should we divide out the equilibria we do not want?
          # Should we do anything about negative populations if "runsteady"?
          rootSolveCounter <- 0
          rootSolveSteadyResult <- NA
          while(rootSolveCounter < 5 && is.na(rootSolveSteadyResult)) {
            rootSolveSteadyResult <- tryCatch(
              rootSolve::steady(
                y = abundance + abs(rnorm(n = length(abundance))) * sqrt(rootSolveCounter),
                func = GeneralisedLotkaVolterra,
                parms = list(a = CommunityMat[community, community],
                             r = Pool$ReproductionRate[community]),
                #pos = TRUE,
                method = "runsteady"
              )$y,
              error = function(e) {
                return(NA)
              })
            rootSolveCounter <- rootSolveCounter + 1
          }
          stopifnot(!is.na(rootSolveSteadyResult))

          statesEncountered$Equilibria[[stateNumber]] <- rootSolveSteadyResult

          statesEncountered$Permanent[stateNumber] <- NA
        }

        # Calculate permanence if necessary.
        if (is.na(statesEncountered$Permanent[stateNumber])) {
          # To calculate permanence, need to consider the powerset of the community.
          # Then, using each retrieved equilibrium from the powerset, perform:
          # minimize z using h_i > 0
          # subject to Sum_i(f_i (x_j) h_i) + z > 0 for each subcommunity j.
          # where f_i is per capita rate of growth of species i,
          #       x_j is a retrieved equilibrium

          # Note that the subcommunities generated are sorted.
          subcommunities <- unlist(
            lapply(1:length(community), function(m, community, ...) {
              lapply(combn(m, ...), function(x) community[x])
            },
            x = 1:length(community), simplify = FALSE, community = community),
            recursive = FALSE
          )

          # Would use lapply, but we may want side effects.

          equilibria <- list()
          for (i in 1:length(subcommunities)) {
            set <- subcommunities[[i]]
            id <- which(toString(set) == statesEncountered$IDs)
            if (length(id) == 0) {
              # A subcommunity not yet observed.
              stateNumberSet <- length(statesEncountered$IDs) + 1

              statesEncountered$IDs[stateNumberSet] <- toString(set)

              parentEquilibrium <- statesEncountered$Equilibria[[stateNumber]]
              parentEquilibrium <- parentEquilibrium[community %in% set]

              # statesEncountered$Equilibria[[stateNumberSet]] <-
              #   rootSolve::steady(
              #     y = parentEquilibrium + 1, # force away from 0's.
              #     func = GeneralisedLotkaVolterra,
              #     parms = list(a = CommunityMat[set, set],
              #                  r = Pool$ReproductionRate[set]),
              #     #pos = TRUE,
              #     method = "runsteady"
              #   )$y

              rootSolveCounter <- 0
              rootSolveSteadyResult <- NA
              while(rootSolveCounter < 5 && is.na(rootSolveSteadyResult)) {
                rootSolveSteadyResult <- tryCatch(
                  rootSolve::steady(
                    y = parentEquilibrium + 1 + abs(rnorm(n = length(parentEquilibrium))) * sqrt(rootSolveCounter),
                    func = GeneralisedLotkaVolterra,
                    parms = list(a = CommunityMat[set, set],
                                 r = Pool$ReproductionRate[set]),
                    #pos = TRUE,
                    method = "runsteady"
                  )$y,
                  error = function(e) {
                    return(NA)
                  })
                rootSolveCounter <- rootSolveCounter + 1
              }
              stopifnot(!is.na(rootSolveSteadyResult))

              statesEncountered$Equilibria[[stateNumberSet]] <- rootSolveSteadyResult

              statesEncountered$Permanent[stateNumberSet] <- NA
              equilibria[[i]] <- statesEncountered$Equilibria[[stateNumberSet]]
            } else {
              # A subcommunity already observed.
              equilibria[[i]] <- statesEncountered$Equilibria[[id]]
            }
          }

          # Adjust equilibria size.
          equilibria <- lapply(
            seq_along(subcommunities),
            function(i, set, nnz, master) {
              retval <- rep(0, length(master))
              setInd <- 1
              masInd <- 1
              while (setInd <= length(set[[i]])) {
                if (set[[i]][setInd] == master[masInd]) {
                  retval[masInd] <- nnz[[i]][setInd]
                  masInd <- masInd + 1
                  setInd <- setInd + 1
                } else if (set[[i]][setInd] < master[masInd]) {
                  setInd <- setInd + 1
                } else {
                  masInd <- masInd + 1
                }
              }
              return(retval)
            },
            set = subcommunities,
            nnz = equilibria,
            master = community
          )


          # Gather coefficients.
          # Row: equilibria
          # Col: per capita growth rate
          coefficients <- lapply(
            equilibria,
            FUN = function(eq, rs, m) {
              # f_col evaluated at equilibrium_row, col = 1:length(community)
              # f_col = r_col + sum()
              rs + m %*% eq
            },
            rs = Pool$ReproductionRate[community],
            m = CommunityMat[community, community]
          )

          # h_i coefficients
          coefficientsMatrix <- matrix(
            unlist(coefficients), byrow = TRUE,
            nrow = length(equilibria), ncol = length(community)
            )

          # add z
          coefficientsMatrix <- cbind(
            coefficientsMatrix,
            rep(1, length(equilibria))
          )

          # add constraint h_i > 0
          coefficientsMatrix <- rbind(
            coefficientsMatrix,
            cbind(diag(ncol = length(community),
                       nrow = length(community)),
                  rep(0, length(community)))
          )

          # add constraint sum(h_i) < 1 (cannot use <=, but equivalent in lp.)
          # same as -sum(h_i) > -1
          coefficientsMatrix <- rbind(
            coefficientsMatrix,
            c(rep(-1, length(community)), 0)
          )

          #TODO Due to numerical problems, we round. There must be a better way.
          # Cannot round without having to round again later, as this has induced
          # different numerical problems (i.e. a result of 1E-10 instead of 0).
          coefficientsMatrix <- round(coefficientsMatrix, 8)


          # Perform optimisation problem.
          output <- lpSolve::lp(
            "min", # minimize
            c(rep(0, length(community)), 1), # 0 * h_1...h_n + z subject to
            coefficientsMatrix,
            ">",
            c(rep(0, length(equilibria)), #sum(fi(xj) * hi) + z > 0
              rep(1E-4, length(community)), # h_i should not be 0.
              -1) # -sum(h_i) > -1
            )

          output <- round(output$solution, 8)

          statesEncountered$Permanent[stateNumber] <-
            if (output[length(output)] <= 0) TRUE else FALSE
        }

        # Check permanence.
        if (statesEncountered$Permanent[stateNumber]) {
          abundance <- statesEncountered$Equilibria[[stateNumber]]

          if ("Sequence" %in% ReturnValues) {
            records$Outcome[eventNumber] <- "Type 2 (Permanent)"
            records$Community[eventNumber] <- I(list(community))
          }
        } else {
          # Need to find an attracting subcommunity.
          # Attracting if permanent and
          # uninvadable by members of the community not in the subcommunity.
          subcommunities <- unlist(
            lapply(1:length(community), function(m, community, ...) {
              lapply(combn(m, ...), function(x) community[x])
            },
            x = 1:length(community), simplify = FALSE, community = community),
            recursive = FALSE
          )

          subcommunitiesPermanent <- rep(NA, length(subcommunities))
          for (scIndex in 1:length(subcommunities)) {
            set <- subcommunities[[scIndex]]
            id <- which(toString(set) == statesEncountered$IDs)

            # Since we have already written down whether the community
            # is permanent or not, we should have all subcommunities
            # listed in the states we have encountered already, alongside
            # their equilibria.
            stopifnot(length(id) == 1)

            if (is.na(statesEncountered$Permanent[id])) {
              # We have not recorded whether this state is permanent or not!
              # Need to perform the permanence check: retrieve all subcommunities
              # check their equilibria, use it to build a coefficients matrix
              # and run a linear program.

              # Note that the subcommunities generated are sorted.
              subsubcommunities <- unlist(
                lapply(1:length(set), function(m, community, ...) {
                  lapply(combn(m, ...), function(x) community[x])
                },
                x = 1:length(set), simplify = FALSE, community = set),
                recursive = FALSE
              )

              # Would use lapply, but we may want side effects.

              equilibria <- list()
              for (i in 1:length(subsubcommunities)) {
                subset <- subsubcommunities[[i]]
                subid <- which(toString(subset) == statesEncountered$IDs)
                # A subsubcommunity already observed.
                equilibria[[i]] <- statesEncountered$Equilibria[[subid]]
              }

              # Adjust equilibria size.
              equilibria <- lapply(
                seq_along(subsubcommunities),
                addMissingEquilibriaEntries,
                set = subsubcommunities,
                nnz = equilibria,
                master = set
              )


              # Gather coefficients.
              # Row: equilibria
              # Col: per capita growth rate
              coefficients <- lapply(
                equilibria,
                FUN = function(eq, rs, m) {
                  # f_col evaluated at equilibrium_row, col = 1:length(community)
                  # f_col = r_col + sum()
                  rs + m %*% eq
                },
                rs = Pool$ReproductionRate[set],
                m = CommunityMat[set, set]
              )

              # h_i coefficients
              coefficientsMatrix <- matrix(
                unlist(coefficients), byrow = TRUE,
                nrow = length(equilibria), ncol = length(set)
              )

              # add z
              coefficientsMatrix <- cbind(
                coefficientsMatrix,
                rep(1, length(equilibria))
              )

              # add constraint h_i > 0
              coefficientsMatrix <- rbind(
                coefficientsMatrix,
                cbind(diag(ncol = length(set),
                           nrow = length(set)),
                      rep(0, length(set)))
              )

              # add constraint sum(h_i) < 1 (cannot use <=, but equivalent in lp.)
              # same as -sum(h_i) > -1
              coefficientsMatrix <- rbind(
                coefficientsMatrix,
                c(rep(-1, length(set)), 0)
              )

              #TODO Due to numerical problems, we round. There must be a better way.
              coefficientsMatrix <- round(coefficientsMatrix, 8)

              # Perform optimisation problem.
              output <- lpSolve::lp(
                "min", # minimize
                c(rep(0, length(set)), 1), # 0 * h_1...h_n + z subject to
                coefficientsMatrix,
                ">",
                c(rep(0, length(equilibria)), #sum(fi(xj) * hi) + z > 0
                  rep(1E-4, length(set)), # h_i should not be 0.
                  -1) # -sum(h_i) > -1
              )

              output <- round(output$solution, 8)

              statesEncountered$Permanent[id] <-
                if (output[length(output)] <= 0) TRUE else FALSE
            }
            subcommunitiesPermanent[scIndex] <- statesEncountered$Permanent[id]
          }

          # Remove subcommunities that are not permanent.
          for (scIndex in length(subcommunitiesPermanent):1) {
            if (!subcommunitiesPermanent[scIndex]) {
              subcommunities[[scIndex]] <- NULL
            }
          }

          # Check for an uninvadable (w.r.t. community) subcommunity.
          subcommunitiesUninvadable <- rep(NA, length(subcommunities))
          for (scIndex in 1:length(subcommunities)) {
            set <- subcommunities[[scIndex]]
            id <- which(toString(set) == statesEncountered$IDs)
            abundanceRow <- statesEncountered$Equilibria[[id]]
            abundanceRow <- addMissingEquilibriaEntries(
              1, list(set), list(abundanceRow), community
            )

            # Add a time column to the abundanceRow
            abundanceRow <- c(NA, abundanceRow)

            subcommunitiesUninvadable[scIndex] <- LawMorton1996_CheckUninvadable(
              AbundanceRow = abundanceRow,
              Pool = Pool[community, ],
              CommunityMatrix = CommunityMat[community, community]
            )
          }

          stopifnot(sum(subcommunitiesUninvadable) >= 1)

          if (sum(subcommunitiesUninvadable) > 1) {
            # Law and Morton 1996 use numerical integration at this point.
            # We have already done this with our steady-state calculation.
            numericalSoln <- which(
              statesEncountered$Equilibria[[stateNumber]] > abundance * 1E-4
            )
            community <- community[numericalSoln]
            abundance <- statesEncountered$Equilibria[[
              which(toString(community) == statesEncountered$IDs)
            ]]

          } else {
            scIndex <- which(subcommunitiesUninvadable)
            community <- subcommunities[[scIndex]]
            abundance <- statesEncountered$Equilibria[[
              which(toString(community) == statesEncountered$IDs)
            ]]
          }

          if ("Sequence" %in% ReturnValues) {
            records$Outcome[eventNumber] <- "Type 3 (!Permanent)"
            records$Community[eventNumber] <- I(list(community))
          }
        }

      } else {
        if ("Sequence" %in% ReturnValues) {
          records$Outcome[eventNumber] <- "Type 1 (Failure)"
          records$Community[eventNumber] <- I(list(community))
        }
      }
    } else {
      if ("Sequence" %in% ReturnValues) {
        records$Outcome[eventNumber] <- "Present"
        records$Community[eventNumber] <- I(list(community))
      }
    }

    # We can add a check when ``sufficiently many'' species have been added.
    # Or, if it is cheap, we can just check immediately if a community is invadable.
    eventNumber <- eventNumber + 1
    uninvadable <- LawMorton1996_CheckUninvadable(
      AbundanceRow = c(NA, addMissingEquilibriaEntries(
        1, list(community), list(abundance), 1:nrow(Pool)
      )),
      Pool = Pool, CommunityMatrix = CommunityMat
    )
    if (eventNumber - 1 > nrow(Pool) &&
        length(community) > 1 &&
        uninvadable
        ) {
      break()
    }
  }

  if (!is.null(seed)) {
    set.seed(oldSeed)
  }

  retval <- list()

  if ("Community" %in% ReturnValues) {
    retval$Community <- community
  }
  if ("Equilibrium" %in% ReturnValues) {
    abundance <- addMissingEquilibriaEntries(
      1, list(community), list(abundance), 1:nrow(Pool)
    )
    retval$Equilibrium <- abundance
  }
  if ("Sequence" %in% ReturnValues) {
    retval$Sequence <- records
  }
  if ("States" %in% ReturnValues) {
    retval$States <- statesEncountered
  }
  if ("Pool" %in% ReturnValues) {
    retval$Pool <- Pool
  }
  if ("Matrix" %in% ReturnValues) {
    retval$Matrix <- CommunityMat
  }

  return(retval)
}
