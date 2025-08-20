PerCapitaDynamics_Mutualistic2 <- function(
  ReproductionRate, InteractionMatrix, NumEnvironments,
  SpeciesTypes, Saturations
) {
  # Variant of PerCapitaDynamics_Mutualistic1 based on the version in T&F2010.
  
  # This is just a workhorse/organising function (factory).
  # No further parametrisation necessary.
  # ReproductionRate = r,
  # InteractionMatrix = Block matrix, Beta blocks on diagonal, Gamma off diag.
  # NumEnvironments = number of patches, as usual
  # Species Types = Species in each Guild (usually length 2: Plant, Pollinator)
  # Subtlety: Mutualism plateaus per-guild but contributions are per species.
  # Saturations = Matrix of species-species saturation values across patchs.
  #               Saturation is amongst substitutes for a function.
  # Note, the transposal of the indices between plants and animals is equivalent
  # to considering a symmetric matrix. I don't see why this symmetry should be
  # exact (cij == cji) as an assumption and so have not included it here.
  # It can still be enforced with the inputs.
  
  # Note also that we have a strict definition of guild. To be a member of a
  # specific guild, you need to have the same interaction types as all other
  # members of the same guild with other guilds. All pollinators within one
  # guild pollinate the same guild of plants (although not necessarily the
  # exact same set of plants). But members of a guild must compete with each
  # other since members of the same species compete with each other.
  # "Intraguild predation" then either requires a multiplex network (not here)
  # or implementation as two separate guilds with otherwise the same
  # interactions amongst other guilds but with members of the top guild eating
  # members of the other focal/intermediate guild.
  
  force(ReproductionRate)
  force(Saturations)
  # See forcing evaluation of function factories.
  # This might be a bug elsewhere (e.g. if arguments aren't evaluated before
  # the produced functions are used). Note that the guards below should render
  # forcings unnecessary.
  
  stopifnot(sum(SpeciesTypes) == length(ReproductionRate),
            length(ReproductionRate) * NumEnvironments == nrow(Saturations),
            nrow(Saturations) == ncol(Saturations),
            nrow(Saturations) == nrow(InteractionMatrix),
            ncol(Saturations) == ncol(InteractionMatrix)
  )
  
  TypePatches <- rep(SpeciesTypes, NumEnvironments)
  st <- c(0, cumsum(TypePatches))
  # Block i, j: (st[i] + 1):st[i+1], (st[j] + 1):st[j+1]
  stind <- seq_along(st[-length(st)])
  
  # Break into blocks
  guilds <- outer(
    stind, stind, Vectorize(
      function(i, j) {list(
        InteractionMatrix[(st[i] + 1):st[i+1], (st[j] + 1):st[j+1]]
      )}
    ))
  
  # Double-check that the guilds rule above is satisfied.
  # There shouldn't be any guild that has both + and - interactions with another
  # guild (including itself). 0's are permitted.
  stopifnot(all(
    unlist(lapply(guilds, function(guild) {
      !(any(sign(guild) == -1) && any(sign(guild) == 1))
    }))
  ))
  
  # Linear Terms
  intraguild <- function(y) {Matrix::bdiag(diag(guilds)) %*% y}
  
  # Saturating Terms
  interguild <- function(i, j, y) {
    # i == j => Intraguild, covered above
    if ( i == j )  return(
      matrix(0, nrow = st[i+1] - st[i], ncol = 1)
    )
    
    ### Differences from Mutualistic1: ###########
    if(any(guilds[i, j][[1]] > 0)) {
      # Gains a benefit, in which case divide amongst providers of benefit.
      
      numerator <- matrix(# Required to make sure structure is matrix in edge case
        guilds[i, j][[1]] *
          Saturations[(st[i] + 1):st[i + 1], (st[j] + 1):st[j + 1]],
        # Hadamard TF2010's alpha and c
        nrow = length((st[i] + 1):st[i + 1]),
        ncol = length((st[j] + 1):st[j + 1])
      )
      
      divisor <- (
        1 + Saturations[(st[i] + 1):st[i + 1], (st[j] + 1):st[j + 1]] * (
          matrix(# Required to make sure structure is matrix in edge case.
            Saturations[(st[i] + 1):st[i + 1], (st[j] + 1):st[j + 1]] > 0,
            nrow = length((st[i] + 1):st[i + 1]),
            ncol = length((st[j] + 1):st[j + 1])
          ) %*%
            y[(st[j] + 1):st[j + 1]]
        )[1:TypePatches[i]] # 1 + alpha[i,j] * sum_{k, alpha[i,k] > 0} (y_k)
      ) # Treats y as a column vector and performs dot product. Results in column.
      # R default multiplication is down columns. We remove column struct to use.
      
    } else if(any(guilds[i, j][[1]] < 0)) {
      # Receives a penalty, in which case divide among possible recipients.
      
      numerator <- matrix(# Required to make sure structure is matrix in edge case
        guilds[i, j][[1]] *
          Saturations[(st[i] + 1):st[i + 1], (st[j] + 1):st[j + 1]],
        # Hadamard TF2010's alpha and c
        nrow = length((st[i] + 1):st[i + 1]),
        ncol = length((st[j] + 1):st[j + 1])
      )
      
      divisor <- (
        1 + Saturations[(st[i] + 1):st[i + 1], (st[j] + 1):st[j + 1]] * t(replicate(
          TypePatches[i], # Should be length((st[i] + 1):st[i + 1]).
          y[(st[i] + 1):st[i + 1]] %*%
            matrix(# Required to make sure structure is matrix in edge case.
              Saturations[(st[i] + 1):st[i + 1], (st[j] + 1):st[j + 1]] > 0,
              nrow = length((st[i] + 1):st[i + 1]),
              ncol = length((st[j] + 1):st[j + 1])
            )
        )[,,] # 1 + alpha[i,j] * sum_{k, alpha[i,k] > 0} (y_k)
        )) # Treats y as a row vector and performs dot product. Results in row.
      # Repeat row structure to turn into matrix (accessed with [,,]).
      
    } else {
      # zeros matrix
      numerator <- matrix(# Required to make sure structure is matrix in edge case.
        0,
        nrow = length((st[i] + 1):st[i + 1]),
        ncol = length((st[j] + 1):st[j + 1])
      )
      divisor <- matrix(# Required to make sure structure is matrix in edge case.
        1,
        nrow = length((st[i] + 1):st[i + 1]),
        ncol = length((st[j] + 1):st[j + 1])
      )
    }
    
    
    mutualism <- (numerator / divisor) %*% y[(st[j] + 1):st[j+1]]
    
    ##############################################
    
    return(mutualism)
  }
  
  interactions <- function(y) {
    inter <- outer(stind, stind, Vectorize(
      function(i, j) {
        interguild(i, j, y)
      },
      SIMPLIFY = FALSE # Causes bugs when all have same dimensions.
    ))
    # Matrix of lists of column vectors
    # Reorganise by row into matrices, then sum for final (per capita) effect.
    # We again run into a simplification problem here in R4.0.3, but there is
    # no argument that we can use to circumvent it until 4.1.0.
    # Circumvent with recursive unlisting.
    inter <- unlist(apply(inter, 1,
                          function(x) list(Matrix::rowSums(do.call(cbind, x)))))
    return(inter + intraguild(y))
  }
  
  # Reconstruct as a single function
  function(t, y, parms = NULL) {
    rep(ReproductionRate, NumEnvironments) + interactions(y)
  }
}
PerCapitaDynamics_MutualisticTF2010 <- PerCapitaDynamics_Mutualistic2