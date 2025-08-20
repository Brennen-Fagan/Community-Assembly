CalculateLocalInvasibles_BareBones <- function(
  Abundance, # Not Time, All Environments
  PerCapitaDynamics, # function(t, y, parms)
  Environments,
  ArrivalDensity,
  DispersalMatrix,
  Tolerance = sqrt(.Machine$double.eps),
  EliminationThreshold = 1E-4
) {
  stopifnot(
    (n <- length(Abundance)) == ncol(DispersalMatrix),
    n == nrow(DispersalMatrix)
  )
  
  # Who are you?
  candidates <- which(
    abs(Abundance - rep(0, n)) < Tolerance
  )
  
  if (length(candidates) == 0) {return(NULL)}
  
  # Will the patch accept you? Escape early if it accepts no-one.
  pcd <- PerCapitaDynamics(0, Abundance, list())
  
  if (!any(pcd[candidates] > 0)) {
    retval <- pcd[candidates] > 0
    names(retval) <- candidates
    return(retval)
  }
  
  # Prepare for mini-assembly process.
  tempDynamics <- compiler::cmpfun(function(t, y, parms) {
    list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
      as.numeric(
        y * PerCapitaDynamics(t, y, parms)
        # Transport: Dispersal means movement of abundance between nodes.
        + DispersalMatrix %*% y
      )
    )
  })
  
  # Time steps: more = longer calculation but smoother.
  tempTimings <- (0:10)/1
  
  invadable <- unlist(lapply(candidates, function(i) {
    if (pcd[i] < 0) return(FALSE)
    
    testAbund <- Abundance
    testAbund[i] <- testAbund[i] + ArrivalDensity
    
    retAbund <- deSolve::lsoda(
      y = testAbund,
      times = tempTimings,
      func = tempDynamics
    )
    
    retAbund[, -1] <- ifelse(retAbund[, -1] <= EliminationThreshold,
                             0, retAbund[, -1])
    
    # Must grow at least once over the timescale and not die out by the end.
    value <- any(diff(retAbund[, i + 1]) >= 0) &
      all(tail(retAbund[, i + 1]) > 0)
    
    names(value) <- NULL
    return(value)
  }))
  
  return(invadable)
}

CalculateLocalInvadables_BareBones <- CalculateLocalInvasibles_BareBones