

Mutualism_CommunityMat <- function(
  Pool,
  MinimumIntraguild = -2/2.5, # While one can try to play with runif vector args
  MaximumIntraguild = -1/2.5, # I advise against it, since they are
  MinimumInterguild = 2/2.5,  # used for sub-blocks at this time.
  MaximumInterguild = 3/2.5,  # Function might need to be rewritten for easier variation.
  IntraspeciesCompetitionMultiplier = 2, # rep(Mult, times = table(Pool$Type))
  seed = NULL
) {
  
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }
  
  # Intraguild interactions
  retmat <- matrix(0, nrow(Pool), nrow(Pool))
  competitors <- outer(Pool$Type, Pool$Type, function(i, j) i == j)
  retmat[competitors] <- runif(sum(competitors),
                               min = MinimumIntraguild,
                               max = MaximumIntraguild)
  
  # Intraspecies (subset of Intraguild) interactions.
  diag(retmat) <- diag(retmat) * IntraspeciesCompetitionMultiplier
  
  # Takes care of block diagonal. Now off-diagonal (interguild) blocks.
  retmat[!competitors] <- runif(sum(!competitors),
                                min = MinimumInterguild,
                                max = MaximumInterguild)
  
  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }
  
  return(retmat)
}