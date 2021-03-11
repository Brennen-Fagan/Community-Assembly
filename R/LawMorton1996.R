# Law and Morton, 1996 #########################################################
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
    oldSeed <- .Random.seed
    set.seed(seed)
  }

  # Assign each species in the pool a body size si.
  # Do so by drawing from a uniform distribution and exponentiating.
  # For all i, j, if si < sj then i may be eaten by j but not vice versa.

  Species <- data.frame(
    ID = 1:(Basal + Consumer),
    Type = c(rep("Basal", Basal), rep("Consumer", Consumer)),
    Size = exp(c(runif(Basal, min = LogBodySize[1], max = LogBodySize[2]),
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

  Species$ReproductionRate[(Basal + 1) : (Basal + Consumer)] <- unlist(lapply(
    Species$Size[(Basal + 1) : (Basal + Consumer)], function(si, k6) {
      p <- -0.1
      return(sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * k6))
    },
    k6 = Parameters[6]
  ))

  if (!is.null(seed)) {
    set.seed(oldSeed)
  }

  Species
}

LawMorton1996_aij <- function(
  Species_i, Species_j, k = c(0.01, 10, 0.5, 0.2, 100, 0.1) # Table 2 values.
) {
  # returns 'p's for the effect of j on i, aij

  # i == j
  if (Species_i$ID == Species_j$ID) {
    if (Species_i$Type == "Basal") {
      return(- Species_i$ReproductionRate * Species_i$Size / k[5])
    } else if(Species_i$Type == "Consumer") {
      return(0)
    } else {
      return(NA)
    }
  }

  # i != j
  if (Species_i$Type == "Basal") {
    if (Species_j$Type == "Basal") {

      # "Basal species are assumed to be independent of one another."
      return(0)

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

LawMorton1996_CommunityMat <- function(Pool, Parameters, seed = NULL) {
  foreach::foreach(
    i = iterators::iter(Pool, by = 'row'), .combine = 'rbind'
  ) %:% foreach::foreach(
    j = iterators::iter(Pool, by = 'row'), .combine = 'c'
  ) %dopar% {
    if (!is.null(seed)) {
      if (exists(".Random.seed")) {
        oldSeed <- .Random.seed
      }
      set.seed(seed + i$ID + j$ID)
    }
    p <- LawMorton1996_aij(i, j, Parameters)
    retval <- ifelse(p == 0,
           0,
           sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * Parameters[6]))
    if (!is.null(seed) & exists("oldSeed")) {
      set.seed(oldSeed)
    }
    retval
  }
}

GeneralisedLotkaVolterra <- function(
  t, y, parms
) {
  # t is the current time point in the integration,
  # y is the current estimate of the variables in the ODE system.
  # parms is a vector or list of parameters
  # The return value of func should be a list, whose first element is a vector
  # containing the derivatives of y with respect to time, and whose next
  # elements are global values that are required at each point in times.
  # The derivatives must be specified in the same order as... y.
  # -- DESolve::rk documentation
  #
  # Our parms contains interaction matrix a and reproduction rates r.

  with(as.list(parms), {
    list(y * (r + a %*% y))
  })
}

LawMorton1996_NumIntegration <- function(
  A, R, X,
  OuterTimeStepSize = 1000,
  InnerTimeStepSize = 1
) {
  # A = matrix of aij's (community matrix)
  # R = vector of ri's (basic reproduction rate)
  # X = vector of initial abundances.
  # OuterTimeStepSize is the length of the solution,
  # InnerTimeStepSize is the time in between records of the solution.

  deSolve::rk(
    X,
    times = seq(from = 0,
                to = OuterTimeStepSize,
                by = InnerTimeStepSize),
    func = GeneralisedLotkaVolterra,
    parms = list(a = A, r = R)
  )
}

LawMorton1996 <- function(
  Basal = NULL,
  Consumer = NULL,
  Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  LogBodySize = c(-2, -1, -1, 0), # c(-2, -1) for Basal, c(-1, 0) for Consumer
  Integrator = "Numerical",
  EliminationThreshold = 10^-4,
  IntegratorTimeStep = 1000,
  ArrivalDensity = 1,
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
      CurrentAbundance_New > EliminationThreshold * CurrentAbundance[SpeciesPresent],
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

LawMorton1996_PlotAbundance <- function(
  Abundance,
  Sequence = NULL,
  guides = FALSE
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
      mapping = ggplot2::aes(xintercept = Time),
      linetype = "dashed",
      color = "black"
      )

    timeDiff <- mean(diff(Sequence$Time), na.rm = TRUE)

    thePlot <- thePlot + ggplot2::geom_label(
      data = Sequence,
      mapping = ggplot2::aes(
        x = Time + timeDiff/3,
        label = IDs
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

LawMorton1996_CheckUninvadable <- function(
  AbundanceRow,
  Pool,
  CommunityMatrix,
  Threshold = 1E-4,
  EliminateSmallPopulations = FALSE
  ) {
  # It is uninvadable by any other species from the pool because the per
  # capita rate of increase of each species absent from the community is
  # negative at the equilibrium point.
  # The per capita rate of increase is the function f_i = r_i + Sum_j(aij xj).
  Abundance <- AbundanceRow[-1]

  if (EliminateSmallPopulations) {
    Abundance[is.na(Abundance) | Abundance <= Threshold] <- 0
  } else {
    Abundance[is.na(Abundance)] <- 0
  }

  notPresent <- !(!is.na(Abundance) & Abundance > Threshold)
  all(Pool$ReproductionRate[notPresent] +
        CommunityMatrix[notPresent, ] %*% Abundance < 0)
}

# LawMorton1996_CheckPermanence <- function(
#   Pool,
#   CommunityMatrix
# ) {
#
# }
