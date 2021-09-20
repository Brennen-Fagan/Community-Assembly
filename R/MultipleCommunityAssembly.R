ArrivalFUN_Example <- function(Events, Rate) {
  cumsum(rexp(Events, Rate))
}
ExtinctFUN_Example <- ArrivalFUN_Example

MultipleNumericalAssembly <- function(
  Pool, # Required from outside function.
  NumEnvironments, # Number of environments
  ComputeInteractionMatrix, # Required outside function.

  Dynamics, # The dynamical system governing the ecosystem interactions.
  EliminationThreshold = 10^-4, # Below which species are removed from internals
  ArrivalDensity = EliminationThreshold * 4 * 10 ^ 3, # Traill et al. 2007

  # See https://en.wikipedia.org/wiki/Coupon_collector%27s_problem#Extensions_and_generalizations
  ArrivalEvents = 10, # Recommend, n = nrow(Pool), n (ln n + 5)
  ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
  ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
  ExtinctEvents = ArrivalEvents, # Note it is very possible to have misses.
  ExtinctRate = NULL, # If NULL, set to arrival rate.
  ExtinctFUN = NULL, # Takes Events and Rate, Returns Event Times.

  SpatialArrangement = NULL, # c("None", "Assembly", "Dispersal", "Diffusion")

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

  # System Setup: ##############################################################

  Species <- nrow(Pool)

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

  ### Space: ###################################################################
  # Space determines the number of slots you need for each species abundance.
  # At minimum we need one slot per species-environment pair.
  # If we are in the diffusion, rather than assembly or dispersal cases,
  # then we possibly have environments that are multiple units long that connect
  # with each other.

  SpatialArrangement <- match.arg(
    tolower(SpatialArrangement),
    tolower(c("none", "assembly", "dispersal", "diffusion"))
  )

  ### Event and Root Set Up: ###################################################
  # Note that this is contingent on space. When in diffusion, need to consider
  # the grid size as well as remove it from all appropriate grid cells.
  # We might be able to add a property to the grid that tells us which
  # environment we are in, which might help if we have interesting spatial
  # arrangements.


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
    func = if (SpatialArrangement == "diffusion") {

    } else {function(t, y, parms,
                     Events = Events,
                     Environments = InteractionMatrices,
                     Pool = Pool) {
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
          if (
            # Pool$ReproductionRate[Events$Species[event]] +
            # Environments[[Events$Environment[event]]] %*% y[
            #   (Events$Environment[event] - 1) * Species + 1:Species
            # ] > 0
            y[abundanceIndex] > 0 ||
            Dynamics(y + c(rep(0, abundanceIndex - 1),
                           .Machine$double.eps,
                           rep(0, length(y) - abundanceIndex)))[abundanceIndex]
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
    },
    maxroot = 10000
  )

  abundance <-
}


