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
  FunctionalResponse = Microbiome_ExampleFunctionalResponseLinear,
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
    Abundance = Abundance,
    InteractionMatrix = InteractionMatrix,
    CarryingCapacity = CarryingCapacity,
    ...
    )

  # NEW IDEA: Normalise a la Coyte et al. 2021.
  Change <- Abundance * (Reproduction + InteractionStrengths)
  if (CarryingCapacity - sum(Pool$Size * Abundance) <= 0) {
    Change <- Change - mean(Change)
  }

  # Note that we should use an event triggered when the carrying capacity is
  # reached so that the system does not exceed it too much.

  return(list(Change))
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
