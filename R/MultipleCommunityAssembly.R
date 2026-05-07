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
  ExtinctionProportion = 1, # Proportion to remove on an "Extinct" event.

  PopulationInitial = NULL, # if NULL, rep(rep(0, nrow(Pool)), NumEnvironments).
  # Otherwise, should be of length nrow(Pool) * NumEnvironments, ordered as
  # each species by environment: (Env 1 Spe 1) (Env 1 Spe 2)... (Env 2 Spe 1)...

  TimeInitial = NULL,
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

  if(!is.null(TimeInitial)) {
    eventsToKeep <- Events$Events$Times >= TimeInitial
    Events$Events <- Events$Events[eventsToKeep, ]
  }

  ### Event and Root Set Up: ###################################################
  # Note that this is contingent on space.
  # For dispersal, we have nodes on a graph that leak into each other.

  deEvents <- list(
    func = EliminationAndNeutralEvents(
      EventsAndSeed = Events, Species = nrow(Pool),
      PerCapitaDynamics = PerCapitaDynamics,
      EliminationThreshold = EliminationThreshold,
      ArrivalDensity = ArrivalDensity,
      ExtinctionProportion = ExtinctionProportion,
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
    if(is.null(TimeInitial)) 0 else TimeInitial,
    Events$Events$Times,
    tail(Events$Events$Times, n = 1) + ReactionTime * 3
  )

  # # Identify and divide the inter distances, then recombine.
  # Timings <- cumsum(c(
  #   0, rep(diff(Timings) / BetweenEventSteps, each = BetweenEventSteps)
  # ))
  # Note: if approximately exactly divisible, a time can occur twice.
  # I am unsure if this can cause a double arrival, but unique should handle it.
  Timings <- unique(unlist(lapply(
    2:length(Timings),
    function(i, ts) {
      seq(from = ts[i - 1], to = ts[i],
          by = min((ts[i] - ts[i - 1]) / BetweenEventSteps, MaximumTimeStep)
      )
    },
    ts = Timings
  )))

  # Bit of an interesting bug. If two consecutive timings are too close, then
  # the solver cannot run. How close is too close? I found a failed run and
  # plot(log10(diff(unique(Timings)))) showed that most are at least 1E-5 of
  # each other. The outlier in this case was 1E-15 apart.
  # It is hard to know how large a threshold to pick from one datum, but
  # sqrt(.Machine$double.eps) is about 1.5E-8, which seems like a good start.

  if (any(diff(Timings) < sqrt(.Machine$double.eps))) {
    # Note that this removes the left entry, which is desired since the
    # right entry is probably the actual event.
    Timings <- Timings[-which(diff(Timings) < sqrt(.Machine$double.eps))]
  }

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
    parms = NULL,
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
    ExtinctionProportion = ExtinctionProportion,
    ArrivalDensity = ArrivalDensity,
    MaximumTimeStep = MaximumTimeStep,
    BetweenEventSteps = BetweenEventSteps,
    DispersalMatrix = DispersalMatrix
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
