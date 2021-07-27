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

Microbiome_ExampleEffectSize <- function(Pool, Matrix, ...) {
  #stopifnot("CarryingCapacity" %in% names(list(...)))

  vals <- rnorm(sum(Matrix > 1), 0, 0.5)
  counter <- 1
  Matrix[Matrix == 1] <- 0

  Matrix[Matrix == 2 | Matrix == 3] <-
    vals[counter - 1 + 1:sum(Matrix == 2 | Matrix == 3)]
  counter <- counter + sum(Matrix == 2 | Matrix == 3)

  Matrix[Matrix == 4] <-
    -abs(vals[counter - 1 + 1:sum(Matrix == 4)]) #/ list(...)$CarryingCapacity
  counter <- counter + sum(Matrix == 4)

  Matrix[Matrix == 5] <-
    abs(vals[counter - 1 + 1:sum(Matrix == 5)])
  counter <- counter + sum(Matrix == 5)

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

  stopifnot(FracExploit > 0)
  stopifnot(FracCompete > 0)
  stopifnot(FracMutual > 0)
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
  Abundance, InteractionMatrix, CarryingCapacity, ...
) {

}

Microbiome_ExampleFunctionalResponseType2 <- function(
  Abundance, InteractionMatrix, CarryingCapacity, ...
) {

}

Microbiome_DynamicsBasic <- function(
  # The intention of this function is that we will call it inside of a solver,
  # such as the deSolve family of functions.
  times,
  # Not used, maintained for deSolve compatibility.
  Abundance,
  # (Column) Vector of the abundance densities.
  parameters,
  # Not used, maintained for deSolve compatibility.
  Pool,
  # Used for the ReproductionRate argument.
  # This exposes more than wanted security-wise, but make sense generically.
  InteractionMatrix,
  # length(Abundance) x length(Abundance), diag = self-regulation
  FunctionalResponse = Microbiome_ExampleFunctionalResponse,
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
  times,
  # Not used, maintained for deSolve compatibility.
  Abundance,
  # (Column) Vector of the abundance densities.
  parameters,
  # Not used, maintained for deSolve compatibility.
  Pool,
  # Used for the ReproductionRate and Size arguments.
  # This exposes more than wanted security-wise, but make sense generically.
  InteractionMatrix,
  # length(Abundance) x length(Abundance), diag = self-regulation
  CarryingCapacity,
  # Density above which the size-biased abundance cannot increase.
  # At that point, we normalise growth by decay.
  FunctionalResponse = Microbiome_ExampleFunctionalResponse,
  # Function receives abundance, interaction matrix, and carrying capacity.
  # The function then returns the effect of each population on the others in
  # a matrix.
  ...
) {
  # Parameter checks.
  stopifnot(CarryingCapacity > 0)

  # Gather calculations.
  Reproduction <- Abundance * Pool$ReproductionRate
  InteractionStrengths <- FunctionalResponse(
    Abundance, InteractionMatrix, CarryingCapacity
    )

  PosR <- ifelse(Reproduction > 0, Reproduction, 0)
  PosI <- ifelse(InteractionStrengths > 0, InteractionStrengths, 0)
  NegR <- ifelse(Reproduction < 0, Reproduction, 0)
  NegI <- ifelse(InteractionStrengths < 0, InteractionStrengths, 0)

  # Determine Density Losses.
  Loss <- NegR + Abundance %*% NegI

  # Determine Density Gains.
  Gain <- PosR + Abundance %*% PosI

  #TODO The Gain and Loss are derivatives, but Space is a
  # density, so the units do not make sense here.

  # Determine amount of space left in simulation.
  Space <- (
    CarryingCapacity
    - sum(Pool$Size * Abundance)
    + sum(Loss$Size * Abundance)
  )

  # If space gains (density * size) > space left, normalise.
  if (sum(Pool$Size * Gain) > Space) {
    Gain * sum(Gain)/Space + Loss
  } else {
    Gain + Loss
  }
}

Microbiome_NumericalAssembly <- function(
  Pool = Microbiome_Species(10),
  Interactions = Microbiome_InteractionMat,
  Dynamics = Microbiome_Dynamics,
  Verbose = FALSE
) {

  # Create a pool.
  # Create interaction matrix.
  # While community non-empty, community invadable, and community not steady...
  # Draw species from pool and add to community.
  # Run community.
  # Note the carrying capacity addresses the mutual benefaction problem.
  # To address the carrying capacity, one normalises the total abundance gain
  # by the total abundance loss so that the sum of the two is zero.
  # Note we are implicitly assuming that everything is the same size.
  # One cannot replace the space occupied by one microbe with two new microbes.
  # Evaluate community size, invadability, and steady-state.

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
