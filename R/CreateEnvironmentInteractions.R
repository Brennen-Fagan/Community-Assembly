CreateEnvironmentInteractions <- function(
  Pool, # Required from outside function.
  NumEnvironments, # Number of environments
  ComputeInteractionMatrix, # Required outside function.
  EnvironmentSeeds = NULL, # If one seed, used to generate seeds for the system.
  ModifyPool = base::identity, # A function that modifies the pool before
  # computing the interaction matrix.
  # Otherwise, we can use seeds equal to the number of environments
  ...
) {
  ### Generate EnvironmentSeeds if either NULL nor of correct length. ##########
  if (is.null(EnvironmentSeeds)) {
    EnvironmentSeeds <- 1E8 * runif(NumEnvironments)
  }
  
  if (length(EnvironmentSeeds) != NumEnvironments) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(EnvironmentSeeds)
    EnvironmentSeeds <- 1E8 * runif(NumEnvironments)
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }
  
  ### Generate InteractionMatrices for each Environment. #######################
  InteractionMatrices <- lapply(
    1:NumEnvironments,
    function(i, seed, pool, ...) {
      if (!is.null(seed[i])) {
        if (exists(".Random.seed"))
          oldSeed <- .Random.seed
        set.seed(seed[i])
      }
      
      # Development Note: there is some ambiguity here as to whether we should
      # simply modify the pool and resample the interaction matrix conditioned
      # on the new pool, or if, instead, we should rescale the existing
      # interactions, sampling only those that fundamentally change.
      # The former (below) is far easier, but the latter is potentially more
      # correct.
      retval <- ComputeInteractionMatrix(ModifyPool(pool), ...)
      
      if (exists("oldSeed")) {
        .Random.seed <<- oldSeed
      }
      
      return(retval)
    },
    seed = EnvironmentSeeds,
    pool = Pool,
    ... = ...
  )
  
  return(list(Mats = InteractionMatrices, Seeds = EnvironmentSeeds))
}