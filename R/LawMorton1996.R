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
      Outcome = factor(NULL, levels = c(
        "Type 1 (Failure)",
        "Type 2 (Permanent)",
        "Type 3 (!Permanent)"))
    )
    eventNumber <- 1
  }

  statesEncountered <- list(
    IDs = list(integer(0)),
    Equilibria = list(NULL), # If does not exist, use NA.
    Permanent = TRUE
    )

  for (ID in ArrivalIDs) {
    if (!(ID %in% community)) {
      # Check if species can increase when rare. Yes -> Add to community.
      if (is.null(community) ||
          Pool$ReproductionRate[ID] +
          sum(CommunityMat[community, community] * abundance, na.rm = TRUE) > 0
          ) { # (Note: || will skip abundance if community is null.)
        # Add species to the community.
        community <- sort(c(community, ID))

        # Check Permanence.

        # Evaluate abundance for future iterations.

      } else {
        if ("Sequence" %in% ReturnValues) {
          records$Outcome[eventNumber] <- "Type 1 (Failure)"
        }
      }
    }


  }

  if (!is.null(seed)) {
    set.seed(oldSeed)
  }

  return(retval)
}
