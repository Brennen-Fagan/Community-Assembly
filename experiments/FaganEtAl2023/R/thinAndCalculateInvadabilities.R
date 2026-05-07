thinAndCalculateInvadabilities <- function(loaded, dyn, dis,
                                           preferred_rows_per_event,
                                           divide_time_by,
                                           burnout) {
  # We can't handle all of the data that we are going to be looking at;
  # a small sample had ~180k rows for ~5.3k events = ~34 rows per event.
  # To reduce it, we will divide time up so that there are about the
  # preferred number of rows per event.
  bythin <- floor((nrow(loaded$Abundance)
                   / nrow(loaded$Events))
                  / preferred_rows_per_event)
  
  loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                           to = nrow(loaded$Abundance),
                                           by = bythin), ]
  
  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- loaded$Abundance[, -1] <
    loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
  loaded$Abundance[, -1][toEliminate] <- 0
  loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divide_time_by
  
  # Identify and grab only the last pre-burn-out time:
  target <- which.max(loaded$Abundance[, 1] > burnout / divide_time_by) - 1
  
  if (target == 0) {
    stop("Beginning of burn-out not found.")
  }
  
  # # Convert to Binary, since we are only interested in richness.
  # loaded$Abundance[, -1] <- loaded$Abundance[, -1] > 0
  
  # loaded$Abundance is now prepared.
  # Now, we calculate richnesses.
  print("prep")
  
  ### Invadability: #####################################################
  invadability <- CalculateLocalInvadables_KnockOn(
    Abundance = loaded$Abundance[target, -1],
    PerCapitaDynamics = dyn,
    Environments = loaded$NumEnvironments,
    ArrivalDensity = loaded$Parameters$ArrivalDensity,
    DispersalMatrix = dis,
    TimeScale = loaded$ReactionTime
  )
  
  binary <- rep(FALSE, length(loaded$Abundance[target, -1]))
  binary[invadability$candidates] <- invadability$invadable
  
  invadabilityMat <- (
    Matrix::drop0(Matrix::Matrix(binary, nrow = 1, ncol = length(binary)))
  )
  
  # candidateRegional = candidates who are not present on any patch
  # invadableRegional = above candidates, but who can invade any patch.
  # Matrix(Invadability$Invadabilities[[1]]$invadability[,], nrow = 100, ncol = 10)
  # Compare Matrix(Invadability$Invadabilities[[1]]$species, nrow = 100, ncol = 10)
  # Then candidates are then the all 1 rows of
  # Matrix(loaded$Abundance[target, -1] < EliminationThreshold, nrow = 100, ncol = 10)
  # and they are invadable if there are any 1's in there corresponding row of
  # Matrix(Invadability$Invadabilities[[1]]$invadability[,], nrow = 100, ncol = 10)
  
  theSpecies <- 1:((ncol(loaded$Abundance) - 1) / loaded$NumEnvironments)
  ### Return Invadabilities: ###############################################
  return(list(
    invadability = invadabilityMat, # Sparse
    time = loaded$Abundance[target, 1], # Not Sparse.
    species = rep(
      theSpecies,
      loaded$NumEnvironments
    ), # i.e. 1 2 3 1 2 3 1 2 3, Not Sparse.
    environment = rep(
      1:loaded$NumEnvironments,
      each = ((ncol(loaded$Abundance) - 1) / loaded$NumEnvironments)
    ), # i.e. 1 1 1 2 2 2 3 3 3, Not Sparse.
    speciesRegional = theSpecies,
    invadabilityRegional = as.numeric(
      apply(Matrix(loaded$Abundance[target, -1] <
                     loaded$Parameters$EliminationThreshold,
                   nrow = 100, ncol = 10), 1, all) *
        apply(Matrix(invadabilityMat[,], nrow = 100, ncol = 10), 1, any)
    ),
    effectAbundance = invadability$effectAbundance,
    effectRichness = invadability$effectRichness,
    effectEstablish = invadability$effectEstablish
  ))
}