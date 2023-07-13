# For mutualism, we are assuming that neutral events remain the same.
# It is the dynamics that we are changing:
#   Producers: Positive Intrinsic Growth Rate
#   Pollinators: Negative Intrinsic Growth Rate
#   Both: intra-specific and intra-guild competition, positive mutualism
# Eq: \dot{x}_i = x_i (r - Sum beta_ij x_j + f(x_i, x))
# where f(x_i, x) = (\sum_j gamma_ij x_j)/(1 + h_i \sum_j gamma_ij x_j)
# using Zhang et al. 2022.
# Note the other major (related) approach is Thebault and Fontaine (2010).
# This file is a derivative of Mutualism.R, where initial testing was done.

Mutualism_species <- function(
  SpeciesTypes,
  MinimumRepRates = c(0, -0.2/2.5),
  MaximumRepRates = c(0.01, 0),
  seed = NULL
) {
  #TODO: Might be nice to add some niche specificity in here.
  # Then competition might be distance to niche and pollinators might vary
  # between specific and general.

  # Assuming a bottom-up assembly procedure.
  # As such, we need some method of starting the invasions.
  # One method is to seed the system with a starting network, but this might
  # lead to uninvasible barren patches.
  # Instead, I'll assume small but positive reproductive rates.
  # Otherwise, I'll stick closely to Thebault and Fontaine.

  stopifnot(length(SpeciesTypes) > 0,
            length(SpeciesTypes) == length(MinimumRepRates),
            length(SpeciesTypes) == length(MaximumRepRates))

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # Assuming Basal/Plant/Producers first.
  r <- c(
    runif(n = SpeciesTypes[1],
          min = MinimumRepRates[1],
          max = MaximumRepRates[1]),
    unlist(lapply(seq_along(SpeciesTypes[-1]), function(i, n, mn, mx) {
      runif(n[i], min = mn[i], max = mx[i])
    },
    n = SpeciesTypes[-1],
    mn = MinimumRepRates[-1],
    mx = MaximumRepRates[-1]
    ))
  )

  Species <- data.frame(
    ID = 1:(sum(SpeciesTypes)),
    Type = c(
      rep("Producer", SpeciesTypes[1]),
      unlist(lapply(seq_along(SpeciesTypes)[-1], function(i) {
        rep(i, SpeciesTypes[i])
      }))
    ),
    ReproductionRate = r
  )

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  return(Species)
}

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
      set.seed(oldSeed)
    }
  }

  return(retmat)
}

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
      set.seed(oldSeed)
    }
  }

  return(retmat)
}

Mutualism_CommunityMat_ByBlock <- function(
  Pool,
  MinimumGuildMatrix, # Entry i,j is the min/max effect of species in guild j on
  MaximumGuildMatrix, # species in guild i. Allows mutual., pred., comp., etc.
  IntraspeciesMultiplier = 2, # rep(Mult, times = table(Pool$Type))
  seed = NULL
) {
  # Requirements:
  stopifnot(length(unique(Pool$Type)) == nrow(MinimumGuildMatrix),
            nrow(MinimumGuildMatrix) == ncol(MinimumGuildMatrix),
            nrow(MinimumGuildMatrix) == ncol(MaximumGuildMatrix),
            nrow(MaximumGuildMatrix) == ncol(MaximumGuildMatrix),
            length(IntraspeciesMultiplier) == 1 ||
              length(IntraspeciesMultiplier) == nrow(Pool$Types))

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # Intraguild interactions
  retmat <- matrix(0, nrow(Pool), nrow(Pool))
  dict <- unique(Pool$Type)

  for (j in 1:ncol(MinimumGuildMatrix)) {
    for (i in 1:nrow(MinimumGuildMatrix)) {
      targets <- outer(Pool$Type, Pool$Type,
                       function(x, y) x == dict[i] & y == dict[j])
      retmat[targets] <- runif(sum(targets),
                               min = MinimumGuildMatrix[i, j],
                               max = MaximumGuildMatrix[i, j])
    }}

  # Intraspecies (subset of Intraguild) interactions.
  diag(retmat) <- diag(retmat) * IntraspeciesMultiplier

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  return(retmat)
}

# test:
# pool <- Mutualism_species(c(5, 10))
# stopifnot(
# # Won't be exact same, former goes diagonal first...
  # all(sign(Mutualism_CommunityMat(pool, seed = 1)) ==
  #       sign(Mutualism_CommunityMat_ByBlock(
  #         pool,
  #         MinimumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
  #                                     c(-2/2.5, 2/2.5, 2/2.5, -2/2.5)),
  #         MaximumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
  #                                     c(-1/2.5, 3/2.5, 3/2.5, -1/2.5)),
  #         seed = 1
  #       )))
# )

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

# Version 2 is slow as in > 2 days for 1 run slow.
PerCapitaDynamics_Mutualistic3 <- function(
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
    guild <- guilds[i, j][[1]]
    if ( i == j || !any(guild != 0))  return(
      matrix(0, nrow = st[i+1] - st[i], ncol = 1)
    )

    ### Differences from Mutualistic1: ###########
    if(any(guild > 0)) {
      # Gains a benefit, in which case divide amongst providers of benefit.

      numerator <- matrix(# Required to make sure structure is matrix in edge case
        guild *
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

    } else if(any(guild < 0)) {
      # Receives a penalty, in which case divide among possible recipients.

      numerator <- matrix(# Required to make sure structure is matrix in edge case
        guild *
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
      return(matrix(0, nrow = st[i+1] - st[i], ncol = 1))
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
PerCapitaDynamics_MutualisticTF2010_2 <- PerCapitaDynamics_Mutualistic3
