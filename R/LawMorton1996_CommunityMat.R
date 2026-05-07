LawMorton1996_CommunityMat <- function(
  Pool, Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), seed = NULL,
  Competition = 0, Mutualism = 0, # Suggestion from Coyte 2015, in [0, 1]
  CompetitionBasal = 0,
  Connectance = 1, # in [0, 1]
  DiagParam = 0,
  ConstrainP = NULL) {
  
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    
    set.seed(seed)
  }
  
  newSeeds <- matrix(runif(
    nrow(Pool) * nrow(Pool)
  ) * 1E8, nrow = nrow(Pool))
  
  theMatrix <- foreach::foreach(
    i = iterators::iter(Pool, by = 'row'), .combine = 'rbind',
    nsrow = iterators::iter(newSeeds, by = 'row')
  ) %:% foreach::foreach(
    j = iterators::iter(Pool, by = 'row'), .combine = 'c',
    nsval = iterators::iter(nsrow)
  ) %dopar% {
    if (!is.na(nsval)) {
      if (exists(".Random.seed")) {
        oSeed <- .Random.seed
      }
      set.seed(nsval)
    }
    p <- LawMorton1996_aij(i, j, Parameters, CompetitionBasal,
                           Connectance, DiagParam)
    
    retval <- ifelse(
      p == 0,
      0,
      if(is.null(ConstrainP)) {
        sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * Parameters[6])
      } else {
        with(
          list(Range = UpdateTruncNormBounds(
            ConstrainP, 0, Inf, abs(p), abs(p) * Parameters[6]
          )),
          sign(p) * rtruncnorm(Range[1], Range[2], abs(p), abs(p) * Parameters[6])
        )
      }
    )
    if (!is.na(nsval) & exists("oSeed")) {
      set.seed(oSeed)
    }
    retval
  }
  
  # Need to handle mutualism, competition outside of the entries, since it
  # requires knowledge of the other entry's sign and a shared probability.
  
  # There are nrow(Pool) * (nrow(Pool) - 1) / 2 Interaction Pairs.
  # Currently, they are all (+, -) exploitation pairs.
  # Following Coyte 2015, we recycle the dists but change the signs.
  probs <- runif(nrow(Pool) * (nrow(Pool) - 1) / 2)
  pindex <- 1
  for (i in 1:(nrow(Pool) - 1)) {
    for (j in (i + 1):nrow(Pool)) {
      # Note i != j.
      if (probs[pindex] < Competition) {
        # Competition
        theMatrix[i, j] <- -abs(theMatrix[i, j])
        theMatrix[j, i] <- -abs(theMatrix[j, i])
      } else if (probs[pindex] < Competition + Mutualism) {
        # Mutualism
        theMatrix[i, j] <- abs(theMatrix[i, j])
        theMatrix[j, i] <- abs(theMatrix[j, i])
      }
      pindex <- pindex + 1
    }
  }
  
  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }
  
  return(theMatrix)
}