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

  # Resolve each arrival.
  for (ID in ArrivalIDs) {
    CurrentAbundance[ID] <- CurrentAbundance[ID] + ArrivalDensity

    if (!(ID %in% SpeciesPresent)) {
      SpeciesPresent[length(SpeciesPresent) + 1] <- ID
    }

    Run_GLV <- LawMorton1996_NumIntegration(
      CommunityMat[SpeciesPresent, SpeciesPresent],
      Pool$ReproductionRate[SpeciesPresent],
      CurrentAbundance[SpeciesPresent],
      OuterTimeStepSize = IntegratorTimeStep,
      InnerTimeStepSize = if(is.null(InnerTimeStepSize)) {
        IntegratorTimeStep
      } else {InnerTimeStepSize}
    )

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

    CurrentAbundance[SpeciesPresent] <- CurrentAbundance_New

    SpeciesPresent <- which(CurrentAbundance > 0)
  }

  retval <- list()
  if ("Abundance" %in% ReturnValues) {
    if (!is.null(InnerTimeStepSize)) {
      retval$Abundance <- TimeAbundances
    } else {
      retval$Abundance <- c((IntegratorTimeStep + 1) * ArrivalEvents - 1,
                            CurrentAbundance)
    }
  }
  if ("Sequence" %in% ReturnValues) {
    retval$Sequence <- data.frame(
      Time = (IntegratorTimeStep + 1) * 0:(ArrivalEvents - 1),
      IDs = ArrivalIDs
    )
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
  ReturnValues = c("Community", "Sequence", "Pool", "Matrix"),
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
      Events = 1:ArrivalEvents,
      Addition = ArrivalIDs,
      Outcome = factor(NA, levels = c(
        "Type 1 (Failure)",
        "Type 2 (Permanent)",
        "Type 3 (!Permanent)"))
    )
    eventNumber <- 1
  }

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
          sum(CommunityMat[community, community] * abundance, na.rm = TRUE) > 0
          ) { # (Note: || will skip abundance if community is null.)
        # Add species to the community. Note that we sort for uniqueness.
        community <- sort(c(community, ID))

        # Store properties.
        stateNumber <- which(toString(community) == statesEncountered$IDs)
        if (length(stateNumber) == 0) {
          stateNumber <- length(statesEncountered$IDs) + 1
          statesEncountered$IDs[stateNumber] <- toString(community)

          statesEncountered$Equilibria[stateNumber] <-
            rootSolve::steady(
              y = rep(1000, length(community)),
              func = GeneralisedLotkaVolterra,
              parms = list(a = CommunityMat[community, community],
                           r = Pool$ReproductionRate[community]),
              positive = TRUE
            )$y

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

              statesEncountered$Equilibria[stateNumberSet] <-
                rootSolve::steady(
                  y = rep(1, length(set)),
                  func = GeneralisedLotkaVolterra,
                  parms = list(a = CommunityMat[set, set],
                               r = Pool$ReproductionRate[set]),
                  positive = TRUE
                )

              statesEncountered$Permanent[stateNumberSet] <- NA
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
          coefficientsMatrix <- round(coefficientsMatrix, 6)


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

          statesEncountered$Permanent[stateNumber] <-
            if (output$solution[length(output$solution)] <= 0) TRUE else FALSE
        }

        # Check permanence.
        if (statesEncountered$Permanent[stateNumber]) {
          if ("Sequence" %in% ReturnValues) {
            records$Outcome[eventNumber] <- "Type 2 (Permanent)"
            eventNumber <- eventNumber + 1
          }

          abundance <- statesEncountered$Equilibria[[stateNumber]]
        } else {
          if ("Sequence" %in% ReturnValues) {
            records$Outcome[eventNumber] <- "Type 3 (!Permanent)"
            eventNumber <- eventNumber + 1
          }

          # Need to find an attracting subcommunity.
          # Attracting if permanent and
          # uninvadable by members of the community not in the subcommunity.
        }

      } else {
        if ("Sequence" %in% ReturnValues) {
          records$Outcome[eventNumber] <- "Type 1 (Failure)"
          eventNumber <- eventNumber + 1
        }
      }
    }


  }

  if (!is.null(seed)) {
    set.seed(oldSeed)
  }

  retval <- records

  return(retval)
}
