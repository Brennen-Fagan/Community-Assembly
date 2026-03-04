checkDominantEigenvalue <- function(x, f) {
  max(abs(eigen(rootSolve::jacobian.full(x, f), only.values = TRUE)$values))
}

checkAllDominantEigenvalues <- function(
  Abundance,
  Events = NULL,
  PerCapitaDynamics,
  DispersalMatrix,
  MaxChangeSize = 1e-5
) {
  # For production version, would need to update to handle t, parms.
  Dynamics <- function(t, y, parms) {
    list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
      as.numeric(
        y * PerCapitaDynamics(t, y, parms)
        # Transport: Dispersal means movement of abundance between nodes.
        + DispersalMatrix %*% y
      )
    )
  }

  if (!is.null(Events)) {
    # Look at the times right before a successful event occurs to simplify.
    CandidateEvents <- Events[Events$Success,]
    Candidates <- which(Abundance[, 1] %in% CandidateEvents$Times) - 1 # rows
    if (length(Candidates) != nrow(CandidateEvents)) {
      warning("Some Candidate Times not found in Events.")
    }
  } else {
    Candidates <- 1:nrow(Abundance) # rows -- not times.
  }
  # Check to see if there was much change at the Candidate times.
  MaxChange <- apply(diff(Abundance[,-1])[Candidates, ], 1,
                     function(x) max(abs(x)))
  Candidates <- Candidates[MaxChange < MaxChangeSize]

  if (is.null(Events)) {
    # Take last value in a cluster (shouldn't be clusters if !is.null(events)).
    Candidates <- Candidates[diff(c(Candidates, Inf)) != 1]
  } else {
    # Above accommodates diffs.
    # Now move to last time to accommodate eigenvalue procedure.
    Candidates <- Candidates + 1
  }

  allDominantEigenvalues <- apply(
    Abundance[Candidates, -1], 1,
    function(x) checkDominantEigenvalue(x, Dynamics)
  )

  return(data.frame(
    Time = Abundance[Candidates, 1],
    Row = Candidates,
    Eigenvalue = allDominantEigenvalues
  ))
}
