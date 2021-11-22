ArrivalFUN_Example <- function(Events, Rate) {
  cumsum(rexp(Events, Rate))
}
ExtinctFUN_Example <- ArrivalFUN_Example

PerCapitaDynamics_Type1 <- function(
  ReproductionRate, InteractionMatrix, NumEnvironments
) {
  # Parms contains r = ReproductionRate, a = InteractionMatrix
  function(t, y, parms = NULL) {
    rep(ReproductionRate, NumEnvironments) + InteractionMatrix %*% y
  }
}

CreateEnvironmentInteractions <- function(
  Pool, # Required from outside function.
  NumEnvironments, # Number of environments
  ComputeInteractionMatrix, # Required outside function.
  EnvironmentSeeds = NULL, # If one seed, used to generate seeds for the system.
  ModifyPool = base::identity, # A function that modifies the pool before
                               # computing the interaction matrix.
  # Otherwise, we can use seeds equal to the number of environments
  ...
) {
  ### Generate EnvironmentSeeds if either NULL nor of correct length. ##########
  if (is.null(EnvironmentSeeds)) {
    EnvironmentSeeds <- 1E8 * runif(NumEnvironments)
  }

  if (length(EnvironmentSeeds) != NumEnvironments) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(EnvironmentSeeds)
    EnvironmentSeeds <- 1E8 * runif(NumEnvironments)
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  ### Generate InteractionMatrices for each Environment. #######################
  InteractionMatrices <- lapply(
    1:NumEnvironments,
    function(i, seed, pool, ...) {
      if (!is.null(seed[i])) {
        if (exists(".Random.seed"))
          oldSeed <- .Random.seed
        set.seed(seed[i])
      }

      retval <- ComputeInteractionMatrix(pool, ...)

      if (exists("oldSeed")) {
        set.seed(oldSeed)
      }

      return(retval)
    },
    seed = EnvironmentSeeds,
    pool = Pool,
    ... = ...
  )

  return(list(Mats = InteractionMatrices, Seeds = EnvironmentSeeds))
}

CreateAssemblySequence <- function(
  Species, NumEnvironments,
  # See https://en.wikipedia.org/wiki/Coupon_collector%27s_problem#Extensions_and_generalizations
  ArrivalEvents = NULL, # Recommend, n = nrow(Pool), n (ln n + 5)
  ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
  ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
  ExtinctEvents = NULL, # Note it is very possible to have misses.
  ExtinctRate = NULL, # If NULL, set to arrival rate.
  ExtinctFUN = NULL, # Takes Events and Rate, Returns Event Times.
  HistorySeed = NULL # Use only one seed, this controls all "external" dynamics.
) {
  ### Set up the History, i.e. arrival and extinction events. ##################
  if (is.null(HistorySeed)) {
    HistorySeed <- 1E8 * runif(1)
  }

  if (!is.null(HistorySeed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }

    set.seed(HistorySeed)

    if (is.null(ArrivalEvents)) {
      ArrivalEvents <- 10
    }
    if (is.null(ExtinctEvents)) {
      ExtinctEvents <- ArrivalEvents
    }

    if (is.null(ArrivalRate)) { # Roughly
      ArrivalRate <- 1 # / min(abs(Re(eigen(mean(InteractionMatrices)))))
    }
    if (is.null(ExtinctRate)) {
      ExtinctRate <- ArrivalRate
    }

    # Retrieve Times.
    if (is.null(ArrivalFUN)) {
      ArrivalFUN <- ArrivalFUN_Example
    }
    if (is.null(ExtinctFUN)) {
      ExtinctFUN <- ExtinctFUN_Example
    }
    Arrivals <- ArrivalFUN(ArrivalEvents, ArrivalRate)
    Extincts <- ExtinctFUN(ExtinctEvents, ExtinctRate)

    # Retrieve Species. Retrieve Locations.
    Events <- data.frame(
      Times = c(Arrivals, Extincts),
      Species = sample.int(Species,
                           size = ArrivalEvents + ExtinctEvents,
                           replace = TRUE),
      Environment = sample.int(NumEnvironments,
                               size = ArrivalEvents + ExtinctEvents,
                               replace = TRUE),
      Type = c(rep("Arrival", ArrivalEvents),
               rep("Extinct", ExtinctEvents)),
      Success = NA
      # Arrival AND Invadable <- TRUE
      # Arrival AND Uninvadable <- FALSE
      # Extinct AND Present     (immediately before extinction) <- TRUE
      # Extinct AND Not Present (immediately before extinction) <- FALSE
    )

    # Place in temporal order.
    Events <- with(Events, Events[order(Times), ])

    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  return(list(Events = Events, Seed = HistorySeed))
}

# Note: This is a function that returns a function!
EliminationAndNeutralEvents <- function(
  EventsAndSeed, Species, #InteractionMatrices, #Pool,
  PerCapitaDynamics, EliminationThreshold,
  ArrivalDensity, Verbose = FALSE
) {
  EventDF <- EventsAndSeed$Events # Other list entry is seed.
  PerCapitaDynams <- PerCapitaDynamics
  ArrivalDens <- ArrivalDensity
  Verb <- Verbose

  function(t, y, parms,
           ReturnEvents = FALSE) {
    if (ReturnEvents) {return(EventDF)}
    if (Verb > 0) {
      if (Verb == 1) {
        print(paste("t:", t))
      } else if (Verb > 1) {
        print(paste(Sys.time(), "t:", t))
      }
    }

    y <- ifelse(y <= EliminationThreshold, 0, y)
    event <- which(t == EventDF$Times)
    if (length(event)) {
      abundanceIndex <- (EventDF$Environment[event] - 1) * Species +
        EventDF$Species[event]
      if (EventDF$Type[event] == "Extinct") {
        # Check if already present for records purposes.
        EventDF$Success[event] <<- y[abundanceIndex] > 0
        y[abundanceIndex] <- 0
      } else if (EventDF$Type[event] == "Arrival") {
        # Check if Uninvadable
        # <=> per capita growth rate > 0
        # <=> lim (epsilon -> 0) Dynamics(y + epsilon) > 0
        # <=> PerCapitaDynamics > 0.
        if (
          y[abundanceIndex] > 0 ||
          PerCapitaDynams(0, y, NULL)[abundanceIndex] > 0
        ) {
          y[abundanceIndex] <- y[abundanceIndex] + ArrivalDens
          EventDF$Success[event] <<- TRUE
        } else {
          EventDF$Success[event] <<- FALSE
        }
      }
    }
    return(y)
  }
}

CreateDispersalMatrix <- function(
  EnvironmentDistances, # Distances, which we will invert. Square Matrix >= 0.
  # Entries are zero if two environments are not connected. No self-loops.
  SpeciesSpeeds # Distances/Times.
  # Result is a matrix of frequencies (1/Times) which will be %*% Abundance.
) {
  NumEnvironments <- nrow(EnvironmentDistances)
  stopifnot(NumEnvironments == ncol(EnvironmentDistances))
  stopifnot(
    tryCatch(diag(EnvironmentDistances), error = function(e) {
      Matrix::diag(EnvironmentDistances)
    }) == 0
  )

  dispersalDiags <- NULL
  # Take i to be the (super/sub) diagonal index, aka "band".
  for (i in (NumEnvironments - 1):-(NumEnvironments - 1)) {
    if (i == 0) next
    # Start in top right band, move to bottom left.
    dispersalDiags <- c(
      dispersalDiags,
      list(c(
        unlist(lapply(
          1:NumEnvironments,
          function(index, offset, mat, vec) {
            if (offset != 0 &&
                nrow(mat) >= index + offset &&
                index + offset >= 1
            )
              if (mat[index, index + offset] != 0) { # To prevent Diffusive Infs.
                1/mat[index, index + offset] * vec
              } else {
                0 * vec
              }
          },
          offset = i, mat = EnvironmentDistances, vec = SpeciesSpeeds
        ))
      ))
    )
  }

  # Amount of travel from j to i is d[i,j]
  # Amount of gain to i from j is d[i,j]
  # Amount of travel from i is d[i,i]
  # This way, we can write the change in y from travel as d %*% y
  # Characterising this as proportions of a population, but that is assuming a
  # normalisation that I do not think is strictly necessary.
  # The matrix is sparse and has colsum = 0, diag < 0, offdiag >= 0.
  ns <- length(SpeciesSpeeds) * NumEnvironments
  ks <- length(SpeciesSpeeds) * c((NumEnvironments - 1):1,
                                  -(1:(NumEnvironments - 1)))
  ks <- ks[!unlist(lapply(dispersalDiags, is.null))]
  dispersalDiags <- dispersalDiags[!unlist(lapply(dispersalDiags, is.null))]

  dispersalMatrix <- Matrix::bandSparse(
    n = ns,
    k = ks,
    diagonals = dispersalDiags# c(
    # list(c(rep(0.0001, length(redCom)),   # Island 2 -> Island 1
    #        rep(0.0001, length(redCom)))), # Island 3 -> Island 2
    # list(rep(0, length(redCom))),         # Island 3 -> Island 1
    # list(c(rep(0.0001, length(redCom)),   # Island 1 -> Island 2
    #        rep(0.0001, length(redCom)))), # Island 2 -> Island 3
    # list(rep(0, length(redCom)))          # Island 1 -> Island 3
    # )
  )
  stopifnot(all(dispersalMatrix >= 0))
  diag(dispersalMatrix) <- -Matrix::colSums(dispersalMatrix)
  stopifnot(isTRUE(all.equal(Matrix::colSums(dispersalMatrix),
                             rep(0, ncol(dispersalMatrix)))))

  return(dispersalMatrix)
}

CalculateTrophicStructure <- function(
  Pool,
  NumEnvironments,
  InteractionMatrices,
  EliminationThreshold
  ) {
  # Borrowing from LM1996-NumPoolCom-FoodWebs-2021-07.Rmd
  nrowPool <- nrow(Pool)
  `%>%` <- magrittr::`%>%`

  # This function should be appliable row-wise to the results.
  # One does need to remove the time column, as usual.
  function(y) {
    # Clean up anything not present.
    y <- ifelse(y <= EliminationThreshold, 0, y)

    # Break y up into its environments.
    EnvsY <- lapply(1:NumEnvironments, function(i, ys) {
      ys[(i - 1) * nrowPool + 1:nrowPool]
    }, ys = y)

    EnvsEdgeVertexLists <- lapply(
      seq_along(EnvsY),
      function(i, envY, mats) {
      # "Red"uced "Com"munity; who's present.
      redCom <- which(envY[[i]] > 0)

      if (length(redCom) == 0) {
        return(list(
          Edges = NA,
          Vertices = NA
        ))
      }

      redPool <- Pool[redCom, ]
      redMat <- matrix(mats[[i]][redCom, redCom],
                       nrow = length(redCom),
                       ncol = length(redCom))

      colnames(redMat) <- paste0('s',as.character(redCom))
      rownames(redMat) <- colnames(redMat)

      names(redPool)[1] <- "node"
      redPool$node <- colnames(redMat)

      Graph <- igraph::graph_from_adjacency_matrix(
        redMat, weighted = TRUE
      )

      Graph <- igraph::set.vertex.attribute(
        Graph, "name", value = colnames(redMat)
      )

      redPool$N <- envY[[i]][redCom]

      # For later analysis, take the matrix diagonal.

      redPool$Intraspecific <- diag(redMat)

      GraphAsDataFrame <- igraph::as_data_frame(Graph)

      # Add in abundances for calculating abundance * (gain or loss)
      GraphAsDataFrame <- dplyr::left_join(
        GraphAsDataFrame,
        dplyr::select(redPool, node, N),
        by = c("to" = "node")
      )

      if (("weight" %in% colnames(GraphAsDataFrame))) {
        # We're in a case where there are edges.
        # In the opposite case, we cannot do this part of the calculation.
        # Split data frame.
        ResCon <- GraphAsDataFrame[GraphAsDataFrame$weight > 0,]
        ConRes <- GraphAsDataFrame[GraphAsDataFrame$weight < 0,]

        # Reorder and rename variables.
        ResCon <- dplyr::select(ResCon,
                                to, from, # resource = to, consumer = from,
                                effectPerUnit = weight, resourceAbund = N)
        ConRes <- dplyr::select(ConRes,
                                to, from, # resource = from, consumer = to,
                                effectPerUnit = weight, consumerAbund = N)
        ResCon <- dplyr::mutate(dplyr::group_by(ResCon, from),
                                effectActual = effectPerUnit * resourceAbund,
                                Type = "Exploit+")
        ConRes <- dplyr::mutate(dplyr::group_by(ConRes, from),
                                effectActual = effectPerUnit * consumerAbund,
                                Type = ifelse(from == to,
                                              "SelfReg-",
                                              "Exploit-"))
      }

      IntriG <- with(redPool, data.frame(
        from = node, #resource = node,
        to = node, #consumer = node,
        effectPerUnit = ifelse(ReproductionRate > 0,
                               ReproductionRate, 0),
        effectActual = ifelse(ReproductionRate > 0,
                              N * ReproductionRate, 0),
        Type = "Intrisc+",
        stringsAsFactors = FALSE))
      IntriL <- with(redPool, data.frame(
        from = node, #resource = node,
        to = node, #consumer = node,
        effectPerUnit = ifelse(ReproductionRate < 0,
                               ReproductionRate, 0),
        effectActual = ifelse(ReproductionRate < 0,
                              N * ReproductionRate, 0),
        Type = "Intrisc-",
        stringsAsFactors = FALSE))

      if (exists("ResCon")) {
        ResCon <- dplyr::select(ResCon, -resourceAbund)
      } else {
        ResCon <- data.frame()
      }
      if (exists("ConRes")) {
        ConRes <- dplyr::select(ConRes, -consumerAbund)
      } else {
        ConRes <- data.frame()
      }

      EdgeDataFrame <- dplyr::bind_rows(
        ResCon, ConRes,
        IntriG, IntriL
      )

      EdgeDataFrame <- EdgeDataFrame %>% dplyr::rename(
        # Empirically speaking, to and from appear reversed.
        # A consumer (from) should have a negative effect on resource (to),
        # but the organisation so far marks it as positive. We fix this.
        tempname = to,
        to = from
      ) %>% dplyr::rename(
        from = tempname
      ) %>% dplyr::filter(
        # Remove placeholder entries
        effectPerUnit != 0
      ) %>% dplyr::mutate(
        # Useful to keep effects separate
        effectSign = sign(effectPerUnit)
      ) %>% dplyr::group_by(
        to, effectSign
      ) %>% dplyr::mutate(
        # Perform the post mortem of the most influential from's
        effectEfficiency = effectPerUnit / sum(effectPerUnit),
        effectNormalised = effectActual / sum(effectActual)
      ) %>% dplyr::arrange(to)

      list(
        Edges = EdgeDataFrame,
        Vertices = redPool
      )

    },
    envY = EnvsY,
    mats = InteractionMatrices$Mats)

    EnvsCheddar <- lapply(EnvsEdgeVertexLists, RMTRCode2::toCheddar)

    EnvsTrophic <- lapply(EnvsCheddar,
                          function(x, weight.by) {
                            if (all(!is.na(x)))
                              cheddar::TrophicLevels(x, weight.by = weight.by)
                            else NA
                            },
                          weight.by = "effectNormalised")

    # In principle, I think these are the two return values.
    # In practice, it seems more useful to return the EdgeVertexLists and
    # the Trophic Levels, given the importance of intraspecific interactions.
    # These are what Cheddar does not capture.

    return(list(
      EdgeVertexLists = EnvsEdgeVertexLists,
      TrophicLevels = EnvsTrophic
    ))
  }

}

MultipleNumericalAssembly_Dispersal <- function(
  Pool, # Required from outside function.
  NumEnvironments, # Number of environments
  CharacteristicRate = NULL, # Eigenvalue of the system about stable fixed point
                             # Provide InteractionMatrices if not known.
  InteractionMatrices = NULL, # List of Environment specific matrices
  Events = NULL, # Dataframe of Times, Species, Environment, Type, and Success.

  PerCapitaDynamics, # The dynamical system governing the ecosystem interactions.
  DispersalMatrix, # Matrix, abundance-conserving species-environ.-travel rates.
  # Note that this can be made using CreateDispersalMatrix.
  EliminationThreshold = 10^-4, # Below which species are removed from internals
  ArrivalDensity = EliminationThreshold * 4 * 10 ^ 3, # Traill et al. 2007

  PopulationInitial = NULL, # if NULL, rep(rep(0, nrow(Pool)), NumEnvironments).
  # Otherwise, should be of length nrow(Pool) * NumEnvironments, ordered as
  # each species by environment: (Env 1 Spe 1) (Env 1 Spe 2)... (Env 2 Spe 1)...

  MaximumTimeStep = 1, # Maximum time solver can proceed without elimination.
  BetweenEventSteps = 10, # Number of steps to reach next event to smooth.
  # Otherwise, we can use seeds equal to the number of environments
  HistorySeed = NULL, # Use only one seed, this controls all "external" dynamics.

  Verbose = FALSE,

  # ARGUMENTS REQUIRED IF ALTERNATIVES NOT PROVIDED:
  # See https://en.wikipedia.org/wiki/Coupon_collector%27s_problem#Extensions_and_generalizations
  ArrivalEvents = NULL, # Recommend, n = nrow(Pool), n (ln n + 5)
  ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
  ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
  ExtinctEvents = NULL, # Note it is very possible to have misses.
  ExtinctRate = NULL, # If NULL, set to arrival rate.
  ExtinctFUN = NULL, # Takes Events and Rate, Returns Event Times.

  ... # Arguments to pass through to ComputeInteractionMatrix or for Spatials.

  # Total run length is determined by the last arrival/extinction event +
  # 3 characteristic time lengths. We use three since exp(-3) < 0.05.
  # > 1/exp(3)
  # [1] 0.04978707
) {

  # System Set Up: #############################################################
  # Time scale
  if (is.null(CharacteristicRate) &&
      is.null(InteractionMatrices)
      ) {
    stop(paste("Need to supply Characteristic Rate or Interaction Matrices",
               "in order to establish time scales and time steps."))
  } else if (!is.null(CharacteristicRate)) {
    ReactionTime <- 1 / CharacteristicRate
  } else {
    # We'll take a guess as to how the eigenvalues of Reactions relate to the
    # characteristic time of the system.
    ReactionTime <- 1 / max(
      unlist(lapply(InteractionMatrices$Mats, function(mat) {
        abs(eigen(mat, only.values = TRUE)$values)
      }))
    )
  }

  # Dataframe of Times, Species, Environment, Type, and Success.
  if (is.null(Events))
    Events <- CreateAssemblySequence(
      Species = nrow(Pool),
      NumEnvironments = NumEnvironments,
      ArrivalEvents = ArrivalEvents,
      ArrivalRate = ArrivalRate,
      ArrivalFUN = ArrivalFUN,
      ExtinctEvents = ExtinctEvents,
      ExtinctRate = ExtinctRate,
      ExtinctFUN = ExtinctFUN,
      HistorySeed = HistorySeed
    )

  ### Event and Root Set Up: ###################################################
  # Note that this is contingent on space.
  # For dispersal, we have nodes on a graph that leak into each other.

  deEvents <- list(
    func = EliminationAndNeutralEvents(
      EventsAndSeed = Events, Species = nrow(Pool),
      PerCapitaDynamics = PerCapitaDynamics,
      EliminationThreshold = EliminationThreshold,
      ArrivalDensity = ArrivalDensity,
      Verbose = Verbose
    )
  )

  # Dynamics: ##################################################################
  Dynamics <- function(t, y, parms) {
    list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
      as.numeric(
        y * PerCapitaDynamics(t, y, parms)
        # Transport: Dispersal means movement of abundance between nodes.
        + DispersalMatrix %*% y
      )
    )
  }

  # Timings: ###################################################################
  # Basic Targets: Start time, Event Times, and Final Settling Time.
  Timings <- c(
    0,
    Events$Events$Times,
    tail(Events$Events$Times, n = 1) + ReactionTime * 3
  )

  # # Identify and divide the inter distances, then recombine.
  # Timings <- cumsum(c(
  #   0, rep(diff(Timings) / BetweenEventSteps, each = BetweenEventSteps)
  # ))
  Timings <- c(0, unlist(lapply(
    2:length(Timings),
    function(i, ts) {
      seq(from = ts[i - 1], to = ts[i],
          by = min((ts[i] - ts[i - 1]) / BetweenEventSteps, MaximumTimeStep)
      )
    },
    ts = Timings
  )))

  deEvents$time <- Timings

  # Evaluation: ################################################################
  if (is.null(PopulationInitial)) {
    PopulationInitial <- rep(0, nrow(Pool) * NumEnvironments)
  }

  if (Verbose > 0) {print("Beginning Evaluation.")}
  abundance <- deSolve::lsoda(
    y = PopulationInitial,
    times = Timings,
    func = Dynamics,
    events = deEvents
  )

  # Return values ##############################################################
  retval <- list()
  retval$Events <- deEvents$func(ReturnEvents = TRUE)
  retval$Abundance <- abundance
  retval$NumEnvironments <- NumEnvironments
  retval$ReactionTime <- ReactionTime
  retval$HistorySeed <- Events$Seed
  retval$Parameters <- list(
    EliminationThreshold = EliminationThreshold,
    ArrivalDensity = ArrivalDensity,
    MaximumTimeStep = MaximumTimeStep,
    BetweenEventSteps = BetweenEventSteps
  )
  retval$Ellipsis <- list(...)

  if (!is.null(ArrivalEvents))
    retval$Parameters$ArrivalEvents <- ArrivalEvents
  if (!is.null(ArrivalRate))
    retval$Parameters$ArrivalRate <- ArrivalRate
  if (!is.null(ArrivalFUN))
    retval$Parameters$ArrivalFUN <- ArrivalFUN
  if (!is.null(ExtinctEvents))
    retval$Parameters$ExtinctEvents <- ExtinctEvents
  if (!is.null(ExtinctRate))
    retval$Parameters$ExtinctRate <- ExtinctRate
  if (!is.null(ExtinctFUN))
    retval$Parameters$ExtinctFUN <- ExtinctFUN
  if (!is.null(InteractionMatrices)) {
    retval$Parameters$EnvironmentSeeds <- InteractionMatrices$Seeds
  }

  return(retval)
}
