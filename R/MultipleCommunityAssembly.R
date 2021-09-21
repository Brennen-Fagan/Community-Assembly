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
  ArrivalEvents = 10, # Recommend, n = nrow(Pool), n (ln n + 5)
  ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
  ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
  ExtinctEvents = ArrivalEvents, # Note it is very possible to have misses.
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
  ArrivalDensity
) {
  EventDF <- EventsAndSeed$Events # Other list entry is seed.
  PerCapitaDynams <- PerCapitaDynamics
  ArrivalDens <- ArrivalDensity
  function(t, y, parms,
           ReturnEvents = FALSE) {
    if (ReturnEvents) {return(EventDF)}

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
  stopifnot(diag(EnvironmentDistances) == 0)

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
                index + offset >= 1)
              1/mat[index, index + offset] * vec
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
  dispersalMatrix <- Matrix::bandSparse(
    n = length(SpeciesSpeeds) * NumEnvironments,
    k = length(SpeciesSpeeds) * c((NumEnvironments - 1):1,
                           -(1:(NumEnvironments - 1))),
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

MultipleNumericalAssembly_Dispersal <- function(
  Pool, # Required from outside function.
  NumEnvironments, # Number of environments
  ComputeInteractionMatrix, # Required outside function.

  PerCapitaDynamics, # The dynamical system governing the ecosystem interactions.
  DispersalMatrix, # Matrix, abundance-conserving species-environ.-travel rates.
  # Note that this can be made using CreateDispersalMatrix.
  EliminationThreshold = 10^-4, # Below which species are removed from internals
  ArrivalDensity = EliminationThreshold * 4 * 10 ^ 3, # Traill et al. 2007

  # See https://en.wikipedia.org/wiki/Coupon_collector%27s_problem#Extensions_and_generalizations
  ArrivalEvents = 10, # Recommend, n = nrow(Pool), n (ln n + 5)
  ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
  ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
  ExtinctEvents = ArrivalEvents, # Note it is very possible to have misses.
  ExtinctRate = NULL, # If NULL, set to arrival rate.
  ExtinctFUN = NULL, # Takes Events and Rate, Returns Event Times.

  BetweenEventSteps = 2, # Number of steps to inject between events to smooth.
  EnvironmentSeeds = NULL, # If one seed, used to generate seeds for the system.
  # Otherwise, we can use seeds equal to the number of environments
  HistorySeed = NULL, # Use only one seed, this controls all "external" dynamics.
  Verbose = FALSE,
  ... # Arguments to pass through to ComputeInteractionMatrix or for Spatials.

  # Total run length is determined by the last arrival/extinction event +
  # 3 characteristic time lengths. We use three since exp(-3) < 0.05.
  # > 1/exp(3)
  # [1] 0.04978707
) {

  # System Set Up: #############################################################
  # List of per environment interaction matrices.
  InteractionMatrices <- CreateEnvironmentInteractions(
    Pool = Pool, NumEnvironments = NumEnvironments,
    ComputeInteractionMatrix = ComputeInteractionMatrix,
    EnvironmentSeeds = EnvironmentSeeds, ...
  )
  # Dataframe of Times, Species, Environment, Type, and Success.
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

  # EXPERIMENTAL: if the root function hits 0 from any direction, then the event
  # function is triggered. This will be triggered quite a bit if the system
  # is diffusing insufficiently to pass the threshold in "one go".
  # The problem comes in that this might block all non-assembly species movement
  # so we need to check to see what happens.

  # NOTE: The function in deEvents probably should be converted to a closure.
  # Let's getting a working model up first though.
  deRootFun <- function (t, y, pars) {
    return(y - EliminationThreshold)
  }
  deEvents <- list(
    root = TRUE,
    func = EliminationAndNeutralEvents(
      EventsAndSeed = Events, Species = nrow(Pool),
      PerCapitaDynamics = PerCapitaDynamics,
      EliminationThreshold = EliminationThreshold,
      ArrivalDensity = ArrivalDensity
    ),
    times = Events$Events$Times,
    maxroot = 10000
  )

  # Computation Times: #########################################################
  # Arrange the InteractionMatrices into a sparse matrix.
  Reactions <- Matrix::bdiag(InteractionMatrices$Mats)
  # We'll take a guess as to how the eigenvalues of Reactions relate to the
  # characteristic time of the system.
  ReactionTime <- 1/max(abs(eigen(Reactions)$values))

  # Dynamics: ##################################################################
  Dynamics <- function(t, y, parms) {
     list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
       y * PerCapitaDynamics(t, y, parms)
       # Transport: Dispersal means movement of abundance between nodes.
       + DispersalMatrix %*% y
     )
  }

  # Timings: ###################################################################
  # Basic Targets: Start time, Event Times, and Final Settling Time.
  Timings <- c(
    0,
    Events$Events$Times,
    tail(Events$Events$Times, n = 1) + ReactionTime * 3
  )

  # Identify and divide the inter distances, then recombine.
  Timings <- cumsum(c(
    0, rep(diff(Timings) / BetweenEventSteps, each = BetweenEventSteps)
  ))

  # Evaluation: ################################################################
  abundance <- deSolve::lsoda(
    y = rep(0, nrow(Pool) * NumEnvironments),
    times = Timings,
    func = Dynamics,
    rootfunc = deRootFun,
    events = deEvents
  )

  return(list(
    Events = deEvents$func(ReturnEvents = TRUE),
    Abundance = abundance,
    EnvironmentSeeds = InteractionMatrices$Seeds,
    HistorySeed = Events$Seed
    ))
}


# MultipleNumericalAssembly <- function(
#   Pool, # Required from outside function.
#   NumEnvironments, # Number of environments
#   ComputeInteractionMatrix, # Required outside function.
#
#   Dynamics, # The dynamical system governing the ecosystem interactions.
#   EliminationThreshold = 10^-4, # Below which species are removed from internals
#   ArrivalDensity = EliminationThreshold * 4 * 10 ^ 3, # Traill et al. 2007
#
#   # See https://en.wikipedia.org/wiki/Coupon_collector%27s_problem#Extensions_and_generalizations
#   ArrivalEvents = 10, # Recommend, n = nrow(Pool), n (ln n + 5)
#   ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
#   ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
#   ExtinctEvents = ArrivalEvents, # Note it is very possible to have misses.
#   ExtinctRate = NULL, # If NULL, set to arrival rate.
#   ExtinctFUN = NULL, # Takes Events and Rate, Returns Event Times.
#
#   SpatialArrangement = NULL, # c("None", "Assembly", "Dispersal", "Diffusion")
#
#   EnvironmentSeeds = NULL, # If one seed, used to generate seeds for the system.
#     # Otherwise, we can use seeds equal to the number of environments
#   HistorySeed = NULL, # Use only one seed, this controls all "external" dynamics.
#   Verbose = FALSE,
#   ... # Arguments to pass through to ComputeInteractionMatrix or for Spatials.
#
#   # Total run length is determined by the last arrival/extinction event +
#   # 3 characteristic time lengths. We use three since exp(-3) < 0.05.
#   # > 1/exp(3)
#   # [1] 0.04978707
# ) {
#
#   # System Setup: ##############################################################
#
#   Species <- nrow(Pool)
#
#   ### Generate EnvironmentSeeds if either NULL nor of correct length. ##########
#   if (!is.null(EnvironmentSeeds)) {
#     EnvironmentSeeds <- 1E8 * runif(NumEnvironments)
#   }
#
#   if (length(EnvironmentSeeds) != NumEnvironments) {
#     if (exists(".Random.seed")) {
#       oldSeed <- .Random.seed
#     }
#     set.seed(EnvironmentSeeds)
#     EnvironmentSeeds <- 1E8 * runif(NumEnvironments)
#     if (exists("oldSeed")) {
#       set.seed(oldSeed)
#     }
#   }
#
#   ### Generate InteractionMatrices for each Environment. #######################
#   InteractionMatrices <- lapply(
#     1:NumEnvironments,
#     function(i, seed, pool, ...) {
#       if (!is.null(seed[i])) {
#         if (exists(".Random.seed"))
#           oldSeed <- .Random.seed
#         set.seed(seed)
#       }
#
#       retval <- ComputeInteractionMatrix(pool, ...)
#
#       if (exists("oldSeed")) {
#         set.seed(oldSeed)
#       }
#
#       return(retval)
#     },
#     seed = EnvironmentSeeds,
#     pool = Pool,
#     ... = ...
#   )
#
#   ### Setup the History, i.e. arrival and extinction events. ###################
#   if (is.null(HistorySeed)) {
#     HistorySeed <- 1E8 * runif(1)
#   }
#
#   if (!is.null(HistorySeed)) {
#     if (exists(".Random.seed")) {
#       oldSeed <- .Random.seed
#     }
#
#     set.seed(HistorySeed)
#
#     if (is.null(ArrivalRate)) { # Roughly
#       ArrivalRate <- 1 / min(abs(Re(eigen(mean(InteractionMatrices)))))
#     }
#     if (is.null(ExtinctRate)) {
#       ExtinctRate <- ArrivalRate
#     }
#
#     # Retrieve Times.
#     if (is.null(ArrivalFUN)) {
#       ArrivalFUN <- ArrivalFUN_Example
#     }
#     if (is.null(ExtinctFUN)) {
#       ExtinctFUN <- ExtinctFUN_Example
#     }
#     Arrivals <- ArrivalFUN(ArrivalEvents, ArrivalRate)
#     Extincts <- ExtinctFUN(ExtinctEvents, ExtinctRate)
#
#     # Retrieve Species. Retrieve Locations.
#     Events <- data.frame(
#       Times = c(Arrivals, Extincts),
#       Species = sample.int(Species,
#                            size = ArrivalEvents + ExtinctEvents,
#                            replace = TRUE),
#       Environment = sample.int(Environments,
#                                size = ArrivalEvents + ExtinctEvents,
#                                replace = TRUE),
#       Type = c(rep("Arrival", ArrivalEvents),
#                rep("Extinct", ExtinctEvents)),
#       Invadable = NA
#     )
#
#     # Place in temporal order.
#     Events <- with(Events, Events[order(Times), ])
#
#     if (exists("oldSeed")) {
#       set.seed(oldSeed)
#     }
#   }
#
#   ### Space: ###################################################################
#   # Space determines the number of slots you need for each species abundance.
#   # At minimum we need one slot per species-environment pair.
#   # If we are in the diffusion, rather than assembly or dispersal cases,
#   # then we possibly have environments that are multiple units long that connect
#   # with each other.
#
#   SpatialArrangement <- match.arg(
#     tolower(SpatialArrangement),
#     tolower(c("none", "assembly", "dispersal", "diffusion"))
#   )
#
#   ### Event and Root Set Up: ###################################################
#   # Note that this is contingent on space. When in diffusion, need to consider
#   # the grid size as well as remove it from all appropriate grid cells.
#   # We might be able to add a property to the grid that tells us which
#   # environment we are in, which might help if we have interesting spatial
#   # arrangements.
#
#
#   # EXPERIMENTAL: if the root function hits 0 from any direction, then the event
#   # function is triggered. This will be triggered quite a bit if the system
#   # is diffusing insufficiently to pass the threshold in "one go".
#   # The problem comes in that this might block all non-assembly species movement
#   # so we need to check to see what happens.
#
#   # NOTE: The function in deEvents probably should be converted to a closure.
#   # Let's getting a working model up first though.
#   deRootFun <- function (t, y, pars) {
#     return(y - EliminationThreshold)
#   }
#   deEvents <- list(
#     root = TRUE,
#     func = if (SpatialArrangement == "diffusion") {
#
#     } else {function(t, y, parms,
#                      Events = Events,
#                      Environments = InteractionMatrices,
#                      Pool = Pool) {
#       y <- ifelse(y <= EliminationThreshold, 0, y)
#       event <- which(t == Events$Times)
#       if (length(event)) {
#         abundanceIndex <- (Events$Environment[event] - 1) * Species +
#           Events$Species[event]
#         if (Events$Type[event] == "Extinct") {
#           y[abundanceIndex] <- 0
#         } else if (Events$Type[event] == "Arrival") {
#           # Check if Uninvadable
#           # <=> per capita growth rate > 0
#           # <=> lim (epsilon -> 0) Dynamics(y + epsilon) > 0
#           if (
#             # Pool$ReproductionRate[Events$Species[event]] +
#             # Environments[[Events$Environment[event]]] %*% y[
#             #   (Events$Environment[event] - 1) * Species + 1:Species
#             # ] > 0
#             y[abundanceIndex] > 0 ||
#             Dynamics(y + c(rep(0, abundanceIndex - 1),
#                            .Machine$double.eps,
#                            rep(0, length(y) - abundanceIndex)))[abundanceIndex]
#           ) {
#             y[abundanceIndex] <- y[abundanceIndex] + ArrivalDensity
#             Events$Invadable[event] <<- TRUE
#           } else {
#             Events$Invadable[event] <<- FALSE
#           }
#         }
#       }
#       return(y)
#     }
#     },
#     maxroot = 10000
#   )
#
#   abundance <-
# }


