FindSteadyStateFromEstimate <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
  Tolerance = 1E-1,
  MaxAttempts = 1E2,
  maxRandVal = 1E6,
  Verbose = FALSE
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }
  
  if (is.character(Populations)) {
    init <- CsvRowSplit(Populations)
  } else {
    init <- Populations
  }
  init_min <- min(c(init, 0))
  init_max <- max(c(init, maxRandVal))
  
  anyZeroOrNotSame <- TRUE # Set T to make at least one look.
  epsilon <- Tolerance
  attempt <- 1
  
  # Run and check to see if anyone dies (bad) or changes (not at steady).
  while (anyZeroOrNotSame && MaxAttempts > attempt) {
    if (Verbose) print(paste0(attempt, ":"))
    init_old <- init
    
    init <- #rootSolve::runsteady(
      deSolve::ode(
        times = c(0, 1E4),
        y = init,
        func = Dynamics,
        parms = list(a = InteractionMatrix[com, com],
                     r = Pool$ReproductionRate[com],
                     epsilon = epsilon)#,
        #positive = TRUE
      )#$y
    init <- init[nrow(init), -1]
    if (Verbose) print(init)
    
    if (any(init < epsilon)) {
      # Someone died, reset to random location.
      if (Verbose) print("Died")
      anyZeroOrNotSame <- TRUE
      init[init < epsilon] <- runif(
        n = sum(init < epsilon),
        min = init_min,
        max = init_max
      )
    } else if (any(round(init / init_old, 1) != 1)) {
      # Not done, keep going.
      if (Verbose) print("Changed")
      anyZeroOrNotSame <- TRUE
    } else {
      anyZeroOrNotSame <- FALSE
    }
    
    attempt <- attempt + 1
  }
  
  if (attempt >= MaxAttempts) {
    warning(paste("Failed to converge after", attempt, "attempts."))
  }
  
  return(init)
}