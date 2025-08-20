PerCapitaDynamics_Mutualistic1 <- function(
  ReproductionRate, InteractionMatrix, NumEnvironments,
  SpeciesTypes, Saturations
) {
  # This is just a workhorse/organising function (factory).
  # No further parametrisation necessary.
  # ReproductionRate = r,
  # InteractionMatrix = Block matrix, Beta blocks on diagonal, Gamma off diag.
  # NumEnvironments = number of patches, as usual
  # Species Types = Species in each Guild (usually length 2: Plant, Pollinator)
  # Subtlety: Mutualism plateaus per-guild in this implementation.
  # Saturations = Matrix of guild-guild saturation values.
  force(ReproductionRate)
  force(Saturations)
  # See forcing evaluation of function factories.
  # This might be a bug elsewhere (e.g. if arguments aren't evaluated before
  # the produced functions are used).
  
  st <- c(0, cumsum(rep(SpeciesTypes, NumEnvironments)))
  # Block i, j: (st[i] + 1):st[i+1], (st[j] + 1):st[j+1]
  stind <- seq_along(st[-length(st)])
  
  # Break into blocks
  guilds <- outer(
    stind, stind, Vectorize(
      function(i, j) {list(
        InteractionMatrix[(st[i] + 1):st[i+1], (st[j] + 1):st[j+1]]
      )}
    ))
  
  # Linear Terms
  intraguild <- function(y) {Matrix::bdiag(diag(guilds)) %*% y}
  
  # Saturating Terms
  interguild <- function(i, j, y) {
    # i == j => Intraguild, covered above
    if ( i == j )  return(
      matrix(0, nrow = st[i+1] - st[i], ncol = 1)
    )
    
    mutualism <- guilds[i, j][[1]] %*% y[(st[j] + 1):st[j+1]]
    
    # Note: risk of Inf * 0; which we try to catch.
    
    mutualism <- ifelse(
      mutualism == 0,
      0,
      mutualism / (1 + Saturations[
        ((i - 1) %% length(SpeciesTypes)) + 1,
        ((j - 1) %% length(SpeciesTypes)) + 1]
        * mutualism)
    )
    
    return(mutualism)
  }
  
  interactions <- function(y) {
    inter <- outer(stind, stind, Vectorize(
      function(i, j) {
        interguild(i, j, y)
      }
    ))
    # Matrix of lists of column vectors
    # Reorganise by row into matrices, then sum for final (per capita) effect.
    inter <- unlist(apply(inter, 1, function(x) rowSums(do.call(cbind, x))))
    return(inter + intraguild(y))
  }
  
  # Reconstruct as a single function
  function(t, y, parms = NULL) {
    rep(ReproductionRate, NumEnvironments) + interactions(y)
  }
}