LawMorton1996_species <- function(
  Basal,
  Consumer,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  LogBodySize = c(-2, -1, -1, 0), # c(-2, -1) for Basal, c(-1, 0) for Consumer
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
    Size = 10^(c(runif(Basal, min = LogBodySize[1], max = LogBodySize[2]),
                 runif(Consumer, min = LogBodySize[3], max = LogBodySize[4]))),
    ReproductionRate = 0
  )

  # For species i,
  # if i is basal, set p to 10^(-1 - 0.25 log10 si),
  # otherwise set p to -0.1.
  # Draw ri from a truncated normal distribution, mean p, std. p k6.
  # (Set the sign of ri to that of p).

  Species$ReproductionRate[1:Basal] <- unlist(lapply(
    Species$Size[1:Basal], function(si, k6) {
      p <- 10 ^ (-1 - 0.25 * log10(si))
      return(sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * k6))
    },
    k6 = Parameters[6]
  ))

  if (Consumer > 0) {
    Species$ReproductionRate[(Basal + 1) : (Basal + Consumer)] <- unlist(lapply(
      Species$Size[(Basal + 1) : (Basal + Consumer)], function(si, k6) {
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

LawMorton1996_aij <- function(
  Species_i, Species_j, k = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  CompetitionBasal = 0, Connectance = 1, DiagParam = 0 # Suggested by Jon.
) {
  # returns 'p's for the effect of j on i, aij

  # i == j
  if (Species_i$ID == Species_j$ID) {
    if (Species_i$Type == "Basal") {
      return(- Species_i$ReproductionRate * Species_i$Size / k[5] + DiagParam)
    } else if(Species_i$Type == "Consumer") {
      return(0 + DiagParam)
    } else {
      return(NA)
    }
  }

  # i != j
  # This is not quite fair; this implements amensalism and commensalism
  # rather than proper connectance reduction in the undirected graph.
  if (runif(1) > Connectance) {
    return(0)
  }

  if (Species_i$Type == "Basal") {
    if (Species_j$Type == "Basal") {

      if (CompetitionBasal > 0) {
        # Without a good idea of what they otherwise would look like, I am
        # just going to borrow the predation effects.
        # This is not without precedent; Coyte et al. 2015 used the same dist.
        # but with different sign pairs for different interaction types.
        # We will use the CompetitionBasal parameter to scale this effect though.
        if (Species_i$Size < Species_j$Size) {
          # Predator on Prey
          # Note negativity.
          return(
            - k[1] * exp(-(log10(k[2] * Species_i$Size / Species_j$Size) / k[3]) ^ 2) * CompetitionBasal
          )

        } else if (Species_i$Size > Species_j$Size) {
          # Prey on Predator
          # Note double negativity in the original, we remove one for competition.
          aji <-
            - k[1] * exp(-(log10(k[2] * Species_j$Size / Species_i$Size) / k[3]) ^ 2)
          return(
            aji * k[4] * Species_j$Size / Species_i$Size * CompetitionBasal
          )
        } else {
          return(0)
        }

      } else {
        # "Basal species are assumed to be independent of one another."
        return(0)
      }

    } else if (Species_j$Type == "Consumer") {

      # average effect of consumer j on victim i.
      # aij = -k1 exp( -[ log10 (k2 si / sj) / k3]^2 ) if si < sj,
      # 0 otherwise
      if (Species_i$Size < Species_j$Size) {
        return(
          - k[1] * exp(-(log10(k[2] * Species_i$Size / Species_j$Size) / k[3]) ^ 2)
        )
      } else {
        return(0)
      }

    } else {
      return(NA)
    }
  } else if (Species_i$Type == "Consumer") {
    if (Species_j$Type == "Basal") {

      # Then average effect of victim j on consumer i
      # aij = - aji k4 sj / si if sj < si, 0 otherwise.
      if (Species_j$Size < Species_i$Size) {
        aji <-
          - k[1] * exp(-(log10(k[2] * Species_j$Size / Species_i$Size) / k[3]) ^ 2)
        return(
          - aji * k[4] * Species_j$Size / Species_i$Size
        )
      } else {
        return(0)
      }

    } else if (Species_j$Type == "Consumer") {
      # Then it comes down to their sizes.

      if (Species_i$Size < Species_j$Size) {
        # Predator on Prey
        return(
          - k[1] * exp(-(log10(k[2] * Species_i$Size / Species_j$Size) / k[3]) ^ 2)
        )

      } else if (Species_i$Size > Species_j$Size) {
        # Prey on Predator
        aji <-
          - k[1] * exp(-(log10(k[2] * Species_j$Size / Species_i$Size) / k[3]) ^ 2)
        return(
          - aji * k[4] * Species_j$Size / Species_i$Size
        )
      } else {
        return(0)
      }

    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}

LawMorton1996_CommunityMat <- function(
  Pool, Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), seed = NULL,
  Competition = 0, Mutualism = 0, # Suggestion from Coyte 2015, in [0, 1]
  CompetitionBasal = 0,
  Connectance = 1, # in [0, 1]
  DiagParam = 0) {

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }

    newSeeds <- matrix(runif(
      nrow(Pool) * nrow(Pool)
    ) * 1E8, nrow = nrow(Pool))
  } else {
    newSeeds <- matrix(nrow = nrow(Pool), ncol = nrow(Pool))
  }

  theMatrix <- foreach::foreach(
    i = iterators::iter(Pool, by = 'row'), .combine = 'rbind',
    nsrow = iterators::iter(newSeeds, by = 'row')
  ) %:% foreach::foreach(
    j = iterators::iter(Pool, by = 'row'), .combine = 'c',
    nsval = iterators::iter(nsrow)
  ) %dopar% {
    if (!is.na(nsval)) {
      if (exists(".Random.seed")) {
        oldSeed <- .Random.seed
      }
      set.seed(nsval)
    }
    p <- LawMorton1996_aij(i, j, Parameters, CompetitionBasal,
                           Connectance, DiagParam)
    retval <- ifelse(p == 0,
                     0,
                     sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * Parameters[6]))
    if (!is.na(nsval) & exists("oldSeed")) {
      set.seed(oldSeed)
    }
    retval
  }

  # Need to handle mutualism, competition outside of the entries, since it
  # requires knowledge of the other entry's sign and a shared probability.
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # There are nrow(Pool) * (nrow(Pool) - 1) / 2 Interaction Pairs.
  # Currently, they are all (+, -) exploitation pairs.
  # Following Coyte 2015, we recycle the dists but change the signs.
  probs <- runif(nrow(Pool) * (nrow(Pool) - 1) / 2)
  pindex <- 1
  for (i in 1:(nrow(Pool) - 1)) {
    for (j in (i + 1):nrow(Pool)) {
      # Note i != j.
      if (probs[pindex] < Competition) {
        # Competition
        theMatrix[i, j] <- -abs(theMatrix[i, j])
        theMatrix[j, i] <- -abs(theMatrix[j, i])
      } else if (probs[pindex] < Competition + Mutualism) {
        # Mutualism
        theMatrix[i, j] <- abs(theMatrix[i, j])
        theMatrix[j, i] <- abs(theMatrix[j, i])
      }
      pindex <- pindex + 1
    }
  }

  if (!is.null(seed) & exists("oldSeed")) {
    set.seed(oldSeed)
  }

  return(theMatrix)
}

LawMorton1996_NumIntegration <- function(
  A, R, X,
  OuterTimeStepSize = 1000,
  InnerTimeStepSize = 1,
  Tolerance = 0
) {
  # A = matrix of aij's (community matrix)
  # R = vector of ri's (basic reproduction rate)
  # X = vector of initial abundances.
  # OuterTimeStepSize is the length of the solution,
  # InnerTimeStepSize is the time in between records of the solution.

  ts <- seq(from = 0,
            to = OuterTimeStepSize,
            by = InnerTimeStepSize)

  deSolve::ode(
    X,
    times = ts,
    func = GeneralisedLotkaVolterra,
    parms = list(a = A, r = R, epsilon = Tolerance),
    events = list(func = function(t, y, parms) {
      y[y < parms$epsilon] <- 0
      y
    }, time = ts)
  )
}

LawMorton1996_PlotAbundance <- function(
  Abundance,
  Sequence = NULL,
  guides = TRUE
) {

  colnames(Abundance) <- c("Time", format(1:(ncol(Abundance) - 1)))

  long <- tidyr::pivot_longer(
    as.data.frame(Abundance),
    2:ncol(Abundance),
    names_to = "Species",
    values_to = "Abundance"
  )

  long <- long[!is.na(long$Abundance), ]

  thePlot <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = Time,
      y = Abundance,
      color = Species
    )
  ) + ggplot2::geom_line(
  )

  if (!is.null(Sequence)) {
    thePlot <- thePlot + ggplot2::geom_vline(
      data = Sequence,
      mapping = ggplot2::aes(xintercept = Events),
      linetype = "dashed",
      color = "black"
    )

    timeDiff <- mean(diff(Sequence$Events), na.rm = TRUE)

    thePlot <- thePlot + ggplot2::geom_label(
      data = Sequence,
      mapping = ggplot2::aes(
        x = Events + timeDiff/3,
        label = Addition
      ),
      y = max(long$Abundance, na.rm = TRUE) * 0.85,
      color = "black"
    )
  }

  if (guides == FALSE) {
    thePlot <- thePlot + ggplot2::guides(color = FALSE)
  }

  return(invisible(thePlot))
}
