EventsTemp <- Events
EventsTemp$Events <- EventsTemp$Events[1:10, ]
result <- theFun(
  pool, NumEnvironments = Environments,
  InteractionMatrices = InteractionMatrices,
  Events = EventsTemp,
  PerCapitaDynamics = PerCapitaDynamics2,
  DispersalMatrix = DispersalMatrix,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  ExtinctionProportion = ExtinctionProportion,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = TRUE
)
