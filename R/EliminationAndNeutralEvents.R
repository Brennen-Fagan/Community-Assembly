# Note: This is a function that returns a function!
EliminationAndNeutralEvents <- function(
  EventsAndSeed, Species, #InteractionMatrices, #Pool,
  PerCapitaDynamics, EliminationThreshold,
  ArrivalDensity, ExtinctionProportion = 1, Verbose = FALSE
) {
  EventDF <- EventsAndSeed$Events # Other list entry is seed.
  PerCapitaDynams <- PerCapitaDynamics
  ArrivalDens <- ArrivalDensity
  Verb <- Verbose
  if (ExtinctionProportion <= 1)
    RemainPercentage <- 1 - ExtinctionProportion
  else
    stop("ExtinctionProportion (= x) unclear. Please keep 0 <= x <= 1.")

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
        y[abundanceIndex] <- RemainPercentage * y[abundanceIndex]
      } else if (EventDF$Type[event] == "Arrival") {
        # Check if Uninvadable
        # <=> per capita growth rate > 0
        # <=> lim (epsilon -> 0) Dynamics(y + epsilon) > 0
        # <=> PerCapitaDynamics > 0.
        if (
          y[abundanceIndex] > 0 ||
          PerCapitaDynams(t, y, parms)[abundanceIndex] > 0
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
