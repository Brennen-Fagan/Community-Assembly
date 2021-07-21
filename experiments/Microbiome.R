Microbiome_Species <- function(
  Species,
  #Numeric.
  Size = 1,
  # Either a function or numeric.
  # Function receives dataframe of species IDs.
  # Should return an augmented dataframe with a Size column.
  ReproductionRate = function(df, ...) {
    df$ReproductionRate <- runif(n = nrow(df))
    df
  },
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

Microbiome_InteractionMat <- function(
  Pool = Microbiome_Species(10),
  FracExploit = 0.5, # Numeric.
  FracCompete = 0.25, # Numeric.
  FracMutual = 0.25, # Numeric.
  Connectance = 1, # Numeric.
  EffectSize = function(Pool, Matrix, ...) {
    stopifnot("CarryingCapacity" %in% names(list(...)))

    vals <- rnorm(sum(Matrix > 1), 0, 0.5)
    counter <- 1
    Matrix[Matrix == 1] <- 0

    Matrix[Matrix == 2 | Matrix == 3] <-
      vals[counter - 1 + 1:sum(Matrix == 2 | Matrix == 3)]
    counter <- counter + sum(Matrix == 2 | Matrix == 3)

    Matrix[Matrix == 4] <-
      -abs(vals[counter - 1 + 1:sum(Matrix == 4)]) / list(...)$CarryingCapacity
    counter <- counter + sum(Matrix == 4)

    Matrix[Matrix == 5] <-
      abs(vals[counter - 1 + 1:sum(Matrix == 5)])
    counter <- counter + sum(Matrix == 5)

    Matrix
  },
  # Either a function or numeric.
  # Function receives Pool and 5-types case matrix.
  # 1: 0/0, 2: +/-, 3: -/+, 4: -/-, 5: +/+, first entry is lower triangular.
  # Returns a Matrix with filled entries.
  SelfRegulation = function(Pool, Matrix, ...) {
    diag(Matrix) <- rnorm(nrow(Matrix), -1, 0.5)
    Matrix
  },
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

Microbiome_Dynamics <- function(
  CarryingCapacity,
  FunctionalResponse = function(x, sum, ...) {x/(20 + sum)}, # Generic type 2.
  FunctionalSum = c("Unique", "SameType", "All")
) {
  # Parameter checks.
  stopifnot(CarryingCapacity > 0)
  FunctionalSum <- match.arg(FunctionalSum[1], c("Unique", "SameType", "All"))

}

Microbiome_NumericalAssembly <- function(
  Pool = Microbiome_Species(),
  Interactions = Microbiome_InteractionMat(),
  Dynamics = Microbiome_Dynamics(),
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
