

Mutualism_saturation <- function(
  SpeciesTypes,
  MinimumSaturation = 0.1, # Note, using runif so vectors propagate as usual.
  MaximumSaturation = 1,
  seed = NULL
) {
  # T + F suggest [0.1, 1], while Z et al. suggest 0.2.
  # (Note: same dims, since T + F use inverse in equation but in assignment.)
  # Set intraguild to Infinity (should be no mutualism!).
  
  stopifnot(length(SpeciesTypes) > 0)
  
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }
  
  retmat <- matrix(
    runif(n = length(SpeciesTypes)^2,
          min = MinimumSaturation, max = MaximumSaturation),
    nrow = length(SpeciesTypes),
    ncol = length(SpeciesTypes)
  )
  
  diag(retmat) <- Inf
  
  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }
  
  return(retmat)
}