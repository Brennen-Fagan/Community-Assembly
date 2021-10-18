Microbiome_ExampleReproductionRate <- function(df, ...) {
  df$ReproductionRate <- runif(n = nrow(df))
  df
}

Microbiome_Species <- function(
  Species,
  #Numeric.
  Size = 1,
  # Either a function or numeric.
  # Function receives dataframe of species IDs.
  # Should return an augmented dataframe with a Size column.
  ReproductionRate = Microbiome_ExampleReproductionRate,
  # Either a function or numeric.
  # Function receives dataframe of species IDs and sizes.
  # Should return an augmented dataframe with ReproductionRate column.
  seed = NULL,
  ... # Other arguments to Size or ReproductionRate
  ) {
  # Parameter checks.
  stopifnot(Species > 0)

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  retval <- data.frame(
    ID = 1:Species
  )

  if (is.function(Size)) {
    retval <- Size(retval, ...)
  } else if (is.numeric(Size)) {
    retval$Size <- Size
  } else {
    stop("Size not numeric or function")
  }

  if (is.function(ReproductionRate)) {
    retval <- ReproductionRate(retval, ...)
  } else if (is.numeric(ReproductionRate)) {
    retval$ReproductionRate <- ReproductionRate
  } else {
    stop("ReproductionRate not numeric or function")
  }

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  stopifnot("ID" %in% names(retval),
            "Size" %in% names(retval),
            "ReproductionRate" %in% names(retval))

  retval
}

Microbiome_ExampleEffectSize <- function(Pool, Matrix, sd = 0.5, ...) {
  #stopifnot("CarryingCapacity" %in% names(list(...)))

  vals <- abs(rnorm(sum(Matrix > 1), 0, sd))
  matTypes <- lapply(1:5, function(i, m) {m == i}, m = Matrix)
  counter <- 1

  Matrix[matTypes[[1]]] <- 0

  Matrix[matTypes[[2]] & lower.tri(Matrix)] <-
    +vals[counter - 1 + 1:(sum(matTypes[[2]])/2)]
  counter <- counter + sum(matTypes[[2]])/2
  Matrix[matTypes[[2]] & upper.tri(Matrix)] <-
    -vals[counter - 1 + 1:(sum(matTypes[[2]])/2)]
  counter <- counter + sum(matTypes[[2]])/2

  Matrix[matTypes[[3]] & lower.tri(Matrix)] <-
    -vals[counter - 1 + 1:(sum(matTypes[[3]])/2)]
  counter <- counter + sum(matTypes[[3]])/2
  Matrix[matTypes[[3]] & upper.tri(Matrix)] <-
    +vals[counter - 1 + 1:(sum(matTypes[[3]])/2)]
  counter <- counter + sum(matTypes[[3]])/2

  Matrix[matTypes[[4]]] <-
    -(vals[counter - 1 + 1:sum(matTypes[[4]])]) #/ list(...)$CarryingCapacity
  counter <- counter + sum(matTypes[[4]])

  Matrix[matTypes[[5]]] <-
    +(vals[counter - 1 + 1:sum(matTypes[[5]])])
  counter <- counter + sum(matTypes[[5]])

  Matrix
}

Microbiome_ExampleSelfRegulation <- function(Pool, Matrix, ...) {
  diag(Matrix) <- rnorm(nrow(Matrix), -1, 0.5)
  Matrix
}

Microbiome_InteractionMat <- function(
  Pool = Microbiome_Species(10),
  FracExploit = 0.5, # Numeric.
  FracCompete = 0.25, # Numeric.
  FracMutual = 0.25, # Numeric.
  Connectance = 1, # Numeric.
  EffectSize = Microbiome_ExampleEffectSize,
  # Either a function or numeric.
  # Function receives Pool and 5-types case matrix.
  # 1: 0/0, 2: +/-, 3: -/+, 4: -/-, 5: +/+, first entry is lower triangular.
  # Returns a Matrix with filled entries.
  SelfRegulation = Microbiome_ExampleSelfRegulation,
  # Either a function or numeric.
  # Function receives Pool and Effects Matrix with a zero diagonal.
  # Returns an effects Matrix with modifications to the diagonal.
  seed = NULL,
  ...
) {
  # Assuming Amensalism, Commensalism do not exist in this system.

  stopifnot("ID" %in% names(Pool),
            "Size" %in% names(Pool),
            "ReproductionRate" %in% names(Pool))
  stopifnot(nrow(Pool) > 1)

  stopifnot(FracExploit >= 0)
  stopifnot(FracCompete >= 0)
  stopifnot(FracMutual >= 0)
  stopifnot(Connectance >= 0, Connectance <= 1)

  # Normalise parameters. This allows ratio specification as well.
  if ((FracSum <- FracExploit + FracCompete + FracMutual) > 1) {
    FracExploit <- FracExploit / FracSum
    FracCompete <- FracCompete / FracSum
    FracMutual <- FracMutual / FracSum
  }

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # Determine types in the lower triangle.
  retval <- matrix(0,
                   nrow = nrow(Pool), ncol = nrow(Pool))

  rowID <- lapply(2:nrow(Pool), seq, to = nrow(Pool))
  colID <- rep(1:(nrow(Pool) - 1), times = unlist(lapply(rowID, length)))
  rowID <- unlist(rowID)
  retval[cbind(rowID, colID)] <- apply(rmultinom(
    length(rowID),
    size = 1, # 1 : NULL, 2:4 : Fracs
    prob = c(1 - Connectance,
             Connectance * c(
               0.5 * FracExploit, # + / -
               0.5 * FracExploit, # - / +
               FracCompete, # - / -
               FracMutual   # + / +
             ))
  ) == 1, MARGIN = 2, which)

  retval <- retval + t(retval)

  retval <- EffectSize(Pool, retval, ...)

  retval <- SelfRegulation(Pool, retval, ...)

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  retval
}

Microbiome_ExampleFunctionalResponseLinear <- function(
  Abundance, InteractionMatrix, ...
) {
  InteractionMatrix %*% Abundance
}

Microbiome_ExampleFunctionalResponseType2 <- function(
  Abundance, InteractionMatrix,
  TypeIndices = NULL,
  # List of boolean indices by type: Exploit+, Exploit-, Compete, Mutual.
  # Note this is column indexed.
  # We do not consider intra- and inter-specific competition to be the same.
  HalfResponse = 1, ...
) {
  # Interpret Competition and Mutualism as the same, but opposite sign.
  # Mutualism: The Benefit i derives from j is
  #   aij * xj / (h + sum(x_sametype))
  # for abundance x, interaction matrix a, half response h.

  # Break into positive/negative exploitation, competition, mutualism.
  if (is.null(TypeIndices)) {
    TypeIndices <- InteractionMatrix
    TypeIndices[TypeIndices == 0] <- NA
    TypeIndices <- sign(TypeIndices)
    TypeIndices <- 2 * TypeIndices + t(TypeIndices)
    # Exploitation+ = 1, Exploitation- = -1, Competition = -3, Mutualism = 3.
    diag(TypeIndices) <- 0 # Exclude diagonal from Competition.
    TypeIndices <- lapply(c(1, -1, # In practice, do not need E-.
                            -3, 3), function(i, m) {
      !is.na(m) & m == i
    }, m = TypeIndices)
    # Note this is column indexed.
  }

  # Type2 is a matrix of entries that scale the Interaction Matrix.
  # For each type, for each row, compute denominator and place in position.
  Type2 <- lapply(
    TypeIndices[c(1, 3, 4)], function(inds, m, orientation) {
    retmat <- apply(matrix(inds, nrow = length(m), ncol = length(m)),
          MARGIN = 1, FUN = function(ind, m, h) {
            val <- sum(m[ind]) + h
            ind[ind] <- 1 / val
            ind
          },
          m = m,
          h = HalfResponse)
  }, m = Abundance)

  # Condense Type2 into a single matrix.
  Type2 <- t(Type2[[1]]) + Type2[[1]] + t(Type2[[2]]) + t(Type2[[3]])
  Type2[Type2 == 0] <- 1 # Do not rescale things that do not have rescalings.

  (InteractionMatrix * Type2) %*% Abundance
}

Microbiome_DynamicsBasic <- function(
  # The intention of this function is that we will call it inside of a solver,
  # such as the deSolve family of functions.
  times = NULL,
  # Not used, maintained for deSolve compatibility.
  Abundance,
  # (Column) Vector of the abundance densities.
  parameters = NULL,
  # Not used, maintained for deSolve compatibility.
  Pool,
  # Used for the ReproductionRate argument.
  # This exposes more than wanted security-wise, but make sense generically.
  InteractionMatrix,
  # length(Abundance) x length(Abundance), diag = self-regulation
  FunctionalResponse = Microbiome_ExampleFunctionalResponseLinear,
  # Example is a generic type 2.
  # Function receives abundance, interaction matrix, and optional args.
  # The function then returns the effect of each population on the others in
  # a matrix.
  ...
) {
  with(as.list(parameters), {list(
    Abundance * (
      Pool$ReproductionRate + FunctionalResponse(
        Abundance, InteractionMatrix, ...
      ) # %*% Abundance # Incorporated into the Functional Response
    )
  )})
}

Microbiome_DynamicsNormalisedCapacity <- function(
  # The intention of this function is that we will call it inside of a solver,
  # such as the deSolve family of functions.
  times = NULL,
  # Not used, maintained for deSolve compatibility.
  Abundance,
  # (Column) Vector of the abundance densities.
  parameters = NULL,
  # Not used, maintained for deSolve compatibility.
  Pool,
  # Used for the ReproductionRate and Size arguments.
  # This exposes more than wanted security-wise, but make sense generically.
  InteractionMatrix,
  # length(Abundance) x length(Abundance), diag = self-regulation
  CarryingCapacity,
  # Density above which the size-biased abundance cannot increase.
  # At that point, we normalise growth by decay.
  FunctionalResponse = Microbiome_ExampleFunctionalResponseLinear,
  # Function receives abundance, interaction matrix, and carrying capacity.
  # The function then returns the effect of each population on the others in
  # a matrix.
  ...
) {
  # Parameter checks.
  stopifnot(CarryingCapacity > 0)

  # Gather calculations.
  Reproduction <- Pool$ReproductionRate
  InteractionStrengths <- FunctionalResponse(
    Abundance = Abundance,
    InteractionMatrix = InteractionMatrix,
    CarryingCapacity = CarryingCapacity,
    ...
    )

  # NEW IDEA: Normalise a la Coyte et al. 2021.
  Change <- Abundance * (Reproduction + InteractionStrengths)
  if (CarryingCapacity - sum(Pool$Size * Abundance) <= 0) {
    Change[Abundance != 0] <- Change[Abundance != 0] - mean(Change[Abundance != 0])
  }

  # Note that we should use an event triggered when the carrying capacity is
  # reached so that the system does not exceed it too much.

  return(list(Change))
}

Microbiome_CheckInvadability <- function(
  Pool, Dynamics, InteractionMatrix, CurrentAbundance, ...
) {
  speciesNum <- nrow(Pool)
  SpeciesPresent <- which(CurrentAbundance > 0)

  # Invadable if any per capita growth rate larger than 0.
  any((Dynamics(Pool = Pool, InteractionMatrix = InteractionMatrix,
                 Abundance = CurrentAbundance + .Machine$double.eps,
                 ...)[[1]]/(
                   CurrentAbundance + .Machine$double.eps
                 ))[!(1:speciesNum %in% SpeciesPresent)] > 0)
}

Microbiome_NumericalAssembly <- function(
  Pool = NULL,
  InteractionMatrix = NULL,
  Dynamics = Microbiome_DynamicsBasic,
  Verbose = FALSE,
  seed = NULL,

  InitialAbundance = NULL,
  EliminationThreshold = 10^-4, # Eliminate if drops by.
  InnerEliminationThreshold = 10^-12, # Eliminate if drops below.
  IntegratorTimeStep = 1000, #
  ArrivalDensity = 0.1,
  # Note, this can really matter for determining what steady state we go to.
  # Further, setting it too low can make it species persist longer, as they
  # do not decline as precipitously (i.e. succeeding EliminationThreshold)
  # even if they do eventually die (i.e. fail InnerEliminationThreshold).
  ArrivalEvents = 60,
  ArrivalSampler = c("rearrange", "iid"),
  InnerTimeStepSize = 100,
  ReturnValues = c("Community", "Equilibrium",
                   "Abundance", "Sequence",
                   "Pool", "Matrix",
                   "Steadystate", "Uninvadable"),

  CheckInvadability = FALSE, # Setting this to TRUE forbids species rescue.
  # Note: There is a serious difference between
  # Microbiome_NumericalAssembly(
  #   seed = 1, Verbose = TRUE, ArrivalEvents = 240, IntegratorTimeStep = 10,
  #   InnerTimeStepSize = 1, CheckInvadability = TRUE)
  # and
  # Microbiome_NumericalAssembly(
  #   seed = 1, Verbose = TRUE, ArrivalEvents = 240, IntegratorTimeStep = 10,
  #   InnerTimeStepSize = 1, CheckInvadability = FALSE)
  # Despite only checking cases where uninvadability should be the result!
  # The "problem" appears to be a result of the time taken for an invader to
  # actually reach an elimination threshold, during which rescue can occur.
  ...
) {

  # Parameters: ################################################################
  if (!is.null(seed)) {
    if (exists(".Random.seed"))
      oldSeed <- .Random.seed
    set.seed(seed)
  }

  # Create a pool.
  if (is.null(Pool)) {
    if (Verbose) print("Creating Pool.")
    Pool <- Microbiome_Species(30)
  }

  speciesNum <- nrow(Pool)

  # Create interaction matrix.
  if (is.null(InteractionMatrix)) {
    if (Verbose) print("Creating Interaction Matrix.")
    InteractionMatrix <- Microbiome_InteractionMat(Pool)
  } else {
    stopifnot(nrow(InteractionMatrix) == speciesNum,
              ncol(InteractionMatrix) == speciesNum)
  }

  ##### Setup: #################################################################
  if (is.null(InitialAbundance)) {
    CurrentAbundance <- rep(0, speciesNum)
  } else {
    stopifnot(length(InitialAbundance) == speciesNum)
    CurrentAbundance <- InitialAbundance
  }

  SpeciesPresent <- which(CurrentAbundance > 0)
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
    if (Verbose) print("Randomising Order.")
    ArrivalIDs <- sample.int(speciesNum, size = ArrivalEvents, replace = TRUE)
  } else if (arrivalSampler == "rearrange") {
    if (Verbose) print("Shuffling Order.")
    ArrivalIDs <- replicate(
      n = ceiling(ArrivalEvents / speciesNum),
      sample.int(speciesNum, replace = FALSE)
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
        "Type 3 (Contract)",
        "Error (Benefaction)")),
      Community = NA
    )
  }

  eventNumber <- 1
  rejections <- 0

  ##### Perform Algorithm: #####################################################
  # Resolve each arrival.
  for (ID in ArrivalIDs) {
    if (Verbose) print(paste("Event:", eventNumber, ID, "arriving."))

    SpeciesPresent_Old <- SpeciesPresent

    if (!(ID %in% SpeciesPresent)) {
      SpeciesPresent[length(SpeciesPresent) + 1] <- ID
    }

    #TODO Make sure the Microbiome_DynamicsNormalisedCapacity case is handled.
    if (CheckInvadability) {
      # Check Invadability
      # Note, we check to prevent cases where the new species dies off.
      # This is a problem when using Microbiome_DynamicsNormalisedCapacity,
      # since the dying off can yield a "reverse" benefaction numerically.
      if (
        !Microbiome_CheckInvadability(
          Pool = Pool[SpeciesPresent, ],
          InteractionMatrix = InteractionMatrix[SpeciesPresent, SpeciesPresent],
          Dynamics = Dynamics,
          CurrentAbundance = CurrentAbundance[SpeciesPresent], ...)
      ) {
        if (Verbose) print("Uninvadable")
        # Undo the previous operation.
        SpeciesPresent <- SpeciesPresent[-length(SpeciesPresent)]
      } else {
        if (Verbose) print("Invadable.")
        # Insert the invader.
        CurrentAbundance[ID] <- CurrentAbundance[ID] + ArrivalDensity
      }
    } else {
      CurrentAbundance[ID] <- CurrentAbundance[ID] + ArrivalDensity
    }

    # Use TryCatch if we start seeing errors out the Wazoo.
    # Note we only pass through the relevant information.
    Run_GLV <- deSolve::ode(
      y = CurrentAbundance[SpeciesPresent],
      times = seq(from = 0, to = IntegratorTimeStep, by = InnerTimeStepSize),
      func = Dynamics,
      parms = list(epsilon = InnerEliminationThreshold),
      events = list(func = function(t, y, parms, ...) {
        y[y < parms$epsilon] <- 0
        y
      }, time = seq(from = 0, to = IntegratorTimeStep, by = InnerTimeStepSize)),
      Pool = Pool[SpeciesPresent, ],
      InteractionMatrix = InteractionMatrix[SpeciesPresent, SpeciesPresent],
      ...
    )

    if (any(is.nan(Run_GLV) | is.infinite(Run_GLV) |
            is.na(Run_GLV) | Run_GLV > 1E10)) {
      warnmsg <- paste("NAN, Inf, NA or Abundance > 10^10 detected.",
                       "Mutual Benefaction likely.",
                       paste0("Rejecting step ", eventNumber, "."),
                       "Consider reducing InnerTimeStepSize.")
      warning(warnmsg)
      if (Verbose) print(warnmsg)

      if (SpeciesPresent[length(SpeciesPresent)] == ID)
        SpeciesPresent <- SpeciesPresent[-length(SpeciesPresent)]
      CurrentAbundance[ID] <- CurrentAbundance[ID] - ArrivalDensity

      eventNumber <- eventNumber + 1

      if ("Sequence" %in% ReturnValues) {
        SequenceRetVal$Outcome[eventNumber] <- "Error (Benefaction)"
        SequenceRetVal$Community[eventNumber] <- I(list(SpeciesPresent))
      }

      atSteadyState <- FALSE

      rejections <- rejections + 1

      next()
    }

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
    if (Verbose) print(paste("Checking steadystate."))
    if (CurrentAbundance_New[length(CurrentAbundance_New)] == 0 &&
        length(CurrentAbundance_New) > 1 && SpeciesPresent[length(SpeciesPresent)] == ID) {
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
    if (Verbose) print(paste(atSteadyState))

    CurrentAbundance[SpeciesPresent] <- CurrentAbundance_New

    SpeciesPresent <- which(CurrentAbundance > 0)
    eventNumber <- eventNumber + 1

    if (Verbose) print(paste("Species present:", paste(SpeciesPresent,
                             collapse = ", ")))

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
        # RMTRCode2::LawMorton1996_CheckUninvadable(
        #   AbundanceRow = c(NA, CurrentAbundance),
        #   Pool = Pool, CommunityMatrix = InteractionMatrix
        # )
        !Microbiome_CheckInvadability(
          Pool = Pool, InteractionMatrix = InteractionMatrix,
          Dynamics = Dynamics, CurrentAbundance = CurrentAbundance, ...)
        # !any((Dynamics(Pool = Pool, InteractionMatrix = InteractionMatrix,
        #               Abundance = CurrentAbundance + .Machine$double.eps,
        #               ...)[[1]]/(
        #                 CurrentAbundance + .Machine$double.eps
        #               ))[!(1:speciesNum %in% SpeciesPresent)] > 0)
    ) {
      break()
    }
  }

  ##### Return: ################################################################
  retval <- list()
  if ("Community" %in% ReturnValues) {
    retval$Community <- SpeciesPresent
  }
  if ("Equilibrium" %in% ReturnValues) {
    retval$Equilibrium <- CurrentAbundance
    retval$Equilibrium[is.na(retval$Equilibrium)] <- 0
  }
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
    retval$Matrix <- InteractionMatrix
  }
  if ("Steadystate" %in% ReturnValues) {
    retval$Steadystate <- atSteadyState
  }
  if ("Uninvadable" %in% ReturnValues) {
    retval$Uninvadable <-
      !Microbiome_CheckInvadability(
        Pool = Pool, InteractionMatrix = InteractionMatrix,
        Dynamics = Dynamics, CurrentAbundance = CurrentAbundance, ...)
  }

  retval$Rejections <- rejections

  if (!is.null(seed)) {
    if (exists("oldSeed"))
      set.seed(oldSeed)
  }

  return(retval)
}

Microbiome_PermanenceAssembly <- function(

) {
  # Create a pool.
  # Create interaction matrix.
  # While community non-empty, community invadable, and community not steady...
  # Draw species from pool and add to community.
  # Evaluate community permanence.
  # If not permanent, compute what it collapses to.
  # Evaluate community size, invadability, and steady-state.
}

Microbiome_SolveForStability <- function(

) {
  # Create a pool.
  # Create interaction matrix.
  # Add all of pool to community.
  # Solve growth rates of community.
  # Check boundedness of community; if unbounded, we discard.
  # Check repulsion of boundary; if not repelling, we discard.
}

Microbiome_CollapseCommunity <- function(

) {
  # Create a pool.
  # Create interaction matrix.
  # Add all of pool to community.
  # Run community.
  # (Do we need to check for extinction? Arguable.)
  # Result is whatever remains after running.
}
