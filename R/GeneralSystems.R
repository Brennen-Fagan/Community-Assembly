# Random Matrices: ############################################################

GeneralSystems_RandomMatrix <- function(n, rfunc, ...) {
  matrix(rfunc(n^2, ...), nrow = n, ncol = n)
}

# Example Usage:
# GeneralSystems_ProportionedMatrix(
#   10,
#   list(
#     function(n){cbind(abs(rnorm(n)), -abs(rnorm(n)))},
#     function(n){cbind(-abs(rexp(n)), abs(rexp(n)))}
#   ),
#   c(0.5, 0.5),
#   list(function(n) {rep(0, n)}), c(1)
# )
GeneralSystems_ProportionedMatrix <- function(
  n, offrfuncs, offps, diagrfuncs, diagps
) {
  cat <- GeneralSystems_CategoricalMatrix(n, offps, diagps)
  return(GeneralSystems_CategoricalToProportioned(cat, offrfuncs, diagrfuncs))
}

GeneralSystems_CategoricalMatrix <- function(
  n, offps, diagps, upper.tri.only = TRUE
) {
  offTypes <- seq_along(offps) %*% rmultinom(
    n = n * (n - 1) / 2,
    size = 1, prob = offps/sum(offps))
  diagTypes <- seq_along(diagps) %*% rmultinom(
    n = n,
    size = 1, prob = diagps/sum(diagps))

  retval <- matrix(NA, nrow = n, ncol = n)
  retval[upper.tri(retval)] <- offTypes
  if (!upper.tri.only) {
    retval[lower.tri(retval)] <- 0
    retval <- retval + t(retval)
  }
  diag(retval) <- diagTypes

  return(retval)
}

GeneralSystems_CategoricalToProportioned <- function(
  categorical, offrfuncs, diagrfuncs
) {
  # Note: ignores the lower.tri.
  stopifnot(is.matrix(categorical))
  stopifnot((n <- nrow(categorical)) == ncol(categorical))
  stopifnot(max(c.u.t. <- categorical[upper.tri(categorical)]) <=
              length(offrfuncs))
  stopifnot(max(c.d. <- diag(categorical)) <= length(diagrfuncs))
  # Transpose of c.u.t.
  c.l.t. <- t(categorical)[lower.tri(categorical)]

  offcalls <- table(categorical[upper.tri(categorical)])
  diagcalls <- table(diag(categorical))

  # 2 columns
  # offvals <- mapply(function(rfunc, n) {rfunc(n)},
  #                   offrfuncs, offcalls, SIMPLIFY = FALSE)
  offvals <- lapply(seq_along(offrfuncs), function(i, rfs, ns) {
    if (i %in% as.numeric(names(ns)))
      rfs[[i]](ns[as.numeric(names(ns)) == i][[1]])
    else NULL
  }, rfs = offrfuncs, ns = offcalls)
  # 1 column
  # diagvals <- mapply(function(rfunc, n) {rfunc(n)},
  #                    diagrfuncs, diagcalls, SIMPLIFY = FALSE)
  diagvals <- lapply(seq_along(diagrfuncs), function(i, rfs, ns) {
    if (i %in% as.numeric(names(ns)))
      rfs[[i]](ns[as.numeric(names(ns)) == i][[1]])
    else NULL
  }, rfs = diagrfuncs, ns = diagcalls)

  retval <- matrix(NA, nrow = n, ncol = n)

  for (type in as.numeric(names(offcalls))) {
    retval[upper.tri(retval)][c.u.t. == type] <- offvals[[type]][, 1]
    retval[lower.tri(retval)][c.l.t. == type] <- offvals[[type]][, 2]
  }
  for (type in as.numeric(names(diagcalls))) {
    diag(retval)[c.d. == type] <- diagvals[[type]]
  }

  return(retval)
}

# Pools with Missing Structures: ##############################################

GeneralSystems_SizePool <- function(
  Species,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  LogBodySize = c(-2, -1, -1, 0), # c(-2, -1) for Basal, c(-1, 0) for Consumer
  seed = NULL
) {
  stopifnot(Species > 0)
  stopifnot(length(Parameters) == 6)

  if (Parameters[4] >= 1) {
    warning("Unrealistic energy consumption efficiency k4: k4 >= 1.")
  }
  stopifnot(Parameters[6] > 0)

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # Assign each species in the pool a body size si.
  # Do so by drawing from a uniform distribution and exponentiating.
  # For all i, j, if si < sj then i may be eaten by j but not vice versa.

  retval <- data.frame(
    ID = 1:(Species),
    # ELSE U.A.R. from over the same range
    Size = 10^(runif(Species, min = min(LogBodySize), max = max(LogBodySize))),
    ReproductionRate = 0
  )

  # For species i,
  # if i is basal, set p to 10^(-1 - 0.25 log10 si),
  # otherwise set p to -0.1.
  # Draw ri from a truncated normal distribution, mean p, std. p k6.
  # (Set the sign of ri to that of p).

  retval$ReproductionRate <- unlist(lapply(
    retval$Size, function(si, k6) {
      p <- 10 ^ (-1 - 0.25 * log10(si))
      return(sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * k6))
    },
    k6 = Parameters[6]
  ))

  # IF (SIZE) i.i.d. from centred and rescaled basal distribution.
  centre <- mean(retval$ReproductionRate)
  maxi <- max(retval$ReproductionRate)
  retval$ReproductionRate <- (retval$ReproductionRate - centre)
  retval$ReproductionRate <- retval$ReproductionRate / max(retval$ReproductionRate)
  retval$ReproductionRate <- retval$ReproductionRate * maxi

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  retval
}

GeneralSystems_BasalConsumerPool <- function(
  Basal,
  Consumer,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  seed = NULL
) {
  stopifnot(Basal > 0)
  stopifnot(Consumer >= 0)
  stopifnot(length(Parameters) == 6)

  if (Parameters[4] >= 1) {
    warning("Unrealistic energy consumption efficiency k4: k4 >= 1.")
  }
  stopifnot(Parameters[6] > 0)

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # Assign each species in the pool a body size si.
  # Do so by drawing from a uniform distribution and exponentiating.
  # For all i, j, if si < sj then i may be eaten by j but not vice versa.

  Species <- data.frame(
    ID = 1:(Basal + Consumer),
    Type = c(rep("Basal", Basal), rep("Consumer", Consumer)),
    ReproductionRate = 0
  )

  # IF (BCS) i.i.d. from either (C) Normal(-0.1, -0.01) or (B) Expon(5.3)
  Species$ReproductionRate[1:Basal] <- rexp(Basal, 5.3)

  if (Consumer > 0) {
    Species$ReproductionRate[(Basal + 1) : (Basal + Consumer)] <- unlist(lapply(
      Species$ID[(Basal + 1) : (Basal + Consumer)], function(i, k6) {
        #NOTE i is a dummy here.
        p <- -0.1
        return(sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * k6))
      },
      k6 = Parameters[6]
    ))
  }

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  Species
}

GeneralSystems_NoStructurePool <- function(
  Species,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  seed = NULL
) {
  stopifnot(Species > 0)
  stopifnot(length(Parameters) == 6)

  if (Parameters[4] >= 1) {
    warning("Unrealistic energy consumption efficiency k4: k4 >= 1.")
  }
  stopifnot(Parameters[6] > 0)

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # ELSE U.A.R from (-0.5, 0.5)
  retval <- data.frame(
    ID = 1:(Species),
    ReproductionRate = runif(Species, -0.5, 0.5)
  )

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  retval
}

# Matrices with Missing Structures: ###########################################

GeneralSystems_BasalConsumerMatrix <- function(
  Pool, seed = NULL
) {
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }

    set.seed(seed)
  }
  #   intraspecific interactions
  # IF (BCS) i.i.d. from either (C) 0 or (B) -reprate * 10^Unif(-2, 0) / 100
  #   interspecific interactions
  # IF (BCS) Basal-Basal Non-interacting, Consumer-* Pred-Prey,
  #   dists abs(Normal(0, 0.5))

  Types <- table(Pool$Type)

  # Offdiag Type 1 == Nothing, Type 2 == Pred-Prey, Type 3 == Prey-Pred
  # Diag Type 1 == Basal Intraspecific, Type 2 == Consumer Intraspecific
  retval <- cbind(
    rbind(
      # B x B
      matrix(1, nrow = Types[1], ncol = Types[1]),
      # C x B
      matrix(2, nrow = Types[2], ncol = Types[1])
    ),
    rbind(
      # B x C
      matrix(2, nrow = Types[1], ncol = Types[2]),
      # C x C
      GeneralSystems_CategoricalMatrix(
        Types[2], offps = c(0.5, 0.5), diagps = 1, FALSE
      ) + 1
    )
  )

  retval <- GeneralSystems_CategoricalToProportioned(
    retval, offrfuncs = list(
      function(n) {cbind(rep(0, n), rep(0, n))},
      # NOTE: effect of column on row. Basals first. col 1 is upper right.
      # => predator on prey first, then prey on predator.
      function(n) {cbind(-abs(rnorm(n, 0, 0.5)), abs(rnorm(n, 0, 0.5)))},
      function(n) {cbind(abs(rnorm(n, 0, 0.5)), -abs(rnorm(n, 0, 0.5)))}
    ),
    diagrfuncs = list(
      function(n) {
        -Pool$ReproductionRate[Pool$Type == "Basal"] * 10^runif(n, -2, 0)
      },
      function(n) {rep(0, n)}
    )
  )

  if (!is.null(seed) & exists("oldSeed")) {
    set.seed(oldSeed)
  }

  return(retval)
}

GeneralSystems_SizeMatrix <- function(
  Pool, Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), seed = NULL
) {
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }

    set.seed(seed)
  }

  #   intraspecific interactions
  # IF (SIZE) i.i.d. from basal distribution (use above reprate!)
  #   interspecific interactions
  # IF (SIZE) Pred-Prey by Size, dists follow LM1996
  retval <- matrix(nrow = nrow(Pool), ncol = nrow(Pool))

  diag(retval) <- - Pool$ReproductionRate * Pool$Size / Parameters[5]

  # Abbreviation from LM 1996
  k <- Parameters

  for (i in 1:(nrow(Pool) - 1)) {
    for(j in (i+1):nrow(Pool)) {
      if (Pool$Size[i] > Pool$Size[j]) {
        prey <- j; pred <- i
      } else {
        prey <- i; pred <- j
      }
      # Effect of Predator on Prey
      retval[prey, pred] <-
        - k[1] * exp(-(log10(k[2] * Pool$Size[prey] / Pool$Size[pred]) / k[3]) ^ 2)
      # Effect of Prey on Predator
      retval[pred, prey] <-
        - retval[prey, pred] * k[4] * Pool$Size[prey] / Pool$Size[pred]
    }
  }

  retval <- matrix(
    sign(retval) * abs(rnorm(nrow(Pool)**2, retval, abs(retval) * k[6])),
    nrow = nrow(Pool), ncol = nrow(Pool)
  )

  if (!is.null(seed) & exists("oldSeed")) {
    set.seed(oldSeed)
  }

  return(retval)
}

GeneralSystems_NoStructureMatrix <- function(
  Pool, Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), seed = NULL
) {
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }

    set.seed(seed)
  }
  # intraspecific interactions
  #    ELSE - reprate / 100 * Unif(0, 1)
  # interspecific interactions
  #    ELSE 1/3 Pred-Prey, 1/3 Prey-Pred, 1/3 Competition,
  #          dists abs(Normal(0, 0.5)).
  retval <- GeneralSystems_CategoricalMatrix(
    n = nrow(Pool), offps = c(1/3, 1/3, 1/3), diagps = 1
  )

  retval <- GeneralSystems_CategoricalToProportioned(
    retval, offrfuncs = list(
      function(n) {cbind(-abs(rnorm(n, 0, 0.5)), -abs(rnorm(n, 0, 0.5)))},
      function(n) {cbind(abs(rnorm(n, 0, 0.5)), -abs(rnorm(n, 0, 0.5)))},
      function(n) {cbind(-abs(rnorm(n, 0, 0.5)), abs(rnorm(n, 0, 0.5)))}
    ),
    diagrfuncs = list(
      function(n) rep(0, n)
    )
  )

  diag(retval) <-
    -Pool$ReproductionRate / Parameters[5] * runif(nrow(Pool), 0, 1)

  if (!is.null(seed) & exists("oldSeed")) {
    set.seed(oldSeed)
  }

  return(retval)
}
