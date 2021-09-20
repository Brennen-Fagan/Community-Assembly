ArrivalFUN_Example <- function(Events, Rate) {
  cumsum(rexp(Events, Rate))
}
ExtinctFUN_Example <- ArrivalFUN_Example

CreateEnvironmentInteractions <- function(
  Pool, # Required from outside function.
  NumEnvironments, # Number of environments
  ComputeInteractionMatrix, # Required outside function.
  EnvironmentSeeds = NULL, # If one seed, used to generate seeds for the system.
  # Otherwise, we can use seeds equal to the number of environments
  ...
) {
  ### Generate EnvironmentSeeds if either NULL nor of correct length. ##########
  if (!is.null(EnvironmentSeeds)) {
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
        set.seed(seed)
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

  return(InteractionMatrices)
}

CreateAssemblySequence <- function(
  # See https://en.wikipedia.org/wiki/Coupon_collector%27s_problem#Extensions_and_generalizations
  ArrivalEvents = 10, # Recommend, n = nrow(Pool), n (ln n + 5)
  ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
  ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
  ExtinctEvents = ArrivalEvents, # Note it is very possible to have misses.
  ExtinctRate = NULL, # If NULL, set to arrival rate.
  ExtinctFUN = NULL, # Takes Events and Rate, Returns Event Times.
  HistorySeed = NULL # Use only one seed, this controls all "external" dynamics.
) {
  ### Setup the History, i.e. arrival and extinction events. ###################
  if (is.null(HistorySeed)) {
    HistorySeed <- 1E8 * runif(1)
  }

  if (!is.null(HistorySeed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }

    set.seed(HistorySeed)

    if (is.null(ArrivalRate)) { # Roughly
      ArrivalRate <- 1 / min(abs(Re(eigen(mean(InteractionMatrices)))))
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
      Environment = sample.int(Environments,
                               size = ArrivalEvents + ExtinctEvents,
                               replace = TRUE),
      Type = c(rep("Arrival", ArrivalEvents),
               rep("Extinct", ExtinctEvents)),
      Invadable = NA
    )

    # Place in temporal order.
    Events <- with(Events, Events[order(Times), ])

    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  return(Events)
}

# Note: This is a function that returns a function!
EliminationAndNeutralEvents <- function(
  Events, InteractionMatrices, Pool, PerCapitaDynamics, EliminationThreshold,
  ArrivalDensity, Species = nrow(Pool)
) {
  function(t, y, parms,
           Events = Events,
           Environments = InteractionMatrices,
           Pool = Pool,
           PerCapitaDynamics = PerCapitaDynamics,
           ArrivalDensity = ArrivalDensity) {
    y <- ifelse(y <= EliminationThreshold, 0, y)
    event <- which(t == Events$Times)
    if (length(event)) {
      abundanceIndex <- (Events$Environment[event] - 1) * Species +
        Events$Species[event]
      if (Events$Type[event] == "Extinct") {
        y[abundanceIndex] <- 0
      } else if (Events$Type[event] == "Arrival") {
        # Check if Uninvadable
        # <=> per capita growth rate > 0
        # <=> lim (epsilon -> 0) Dynamics(y + epsilon) > 0
        # <=> PerCapitaDynamics > 0.
        if (
          y[abundanceIndex] > 0 ||
          PerCapitaDynamics(y)[abundanceIndex] > 0
        ) {
          y[abundanceIndex] <- y[abundanceIndex] + ArrivalDensity
          Events$Invadable[event] <<- TRUE
        } else {
          Events$Invadable[event] <<- FALSE
        }
      }
    }
    return(y)
  }
}

CreateDispersalMatrix <- function(
  EnvironmentDistances, # Distances, which we will invert. Matrix >= 0.
  # Entries are zero if two environments are not connected.
  SpeciesSpeeds # Distances/Times.
  # Result is a matrix of frequencies (1/Times) which will be %*% Abundance.
) {
  # So the question is how to format the diagonals then.
  dispersalDiags <- NULL
  for (i in (length(Communities) - 1):-(length(Communities) - 1)) {
    if (i == 0) next
    # Start in top right band, move to bottom left.
    dispersalDiags <- c(
      dispersalDiags,
      list(c(
        unlist(lapply(
          1:length(Communities),
          function(index, offset, mat, vec) {
            if (offset != 0 &&
                nrow(mat) >= index + offset &&
                index + offset >= 1)
              mat[index, index + offset] * vec
          },
          offset = i, mat = DispersalIsland, vec = redDisPool
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
    n = length(redCom) * length(Communities),
    k = length(redCom) * c((length(Communities) - 1):1,
                           -(1:(length(Communities) - 1))),
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
  stopifnot(all(Matrix::colSums(dispersalMatrix) == 0))
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
    Pool, NumEnvironments, ComputeInteractionMatrix,
    EnvironmentSeeds, ...
  )
  # Dataframe of Times, Species, Environment, Type, and Invadable.
  Events <- CreateAssemblySequence(
    ArrivalEvents, ArrivalRate, ArrivalFUN,
    ExtinctEvents, ExtinctRate, ExtinctFUN,
    HistorySeed
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
      Events, InteractionMatrices, Pool,
      PerCapitaDynamics, EliminationThreshold, ArrivalDensity
    ),
    maxroot = 10000
  )

  # Computation Times: #########################################################
  # Arrange the InteractionMatrices into a sparse matrix.
  Reactions <- Matrix::bdiag(InteractionMatrices)
  # We'll take a guess as to how the eigenvalues of Reactions relate to the
  # characteristic time of the system.
  ReactionTime <- 1/max(abs(eigen(Reactions)$values))

  # Dynamics: ##################################################################
  Dynamics <- function(t, y, parms) {
     list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
       y * PerCapitaDynamics(y)
       # Transport: Dispersal means movement of abundance between nodes.
       + DispersalMatrix %*% y
     )
  }

  # Evaluation: ################################################################
  abundance <- deSolve::lsoda(
    y = rep(0, nrow(Pool) * NumEnvironments),
    times = seq(from = 0, by = ReactionTime,
                to = max(Events$Times) + ReactionTime * 3),
    func = Dynamics,
    rootfunc = deRootFun,
    events = deEvents
  )

  return(list(
    Events = Events,
    Abundance = abundance
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


