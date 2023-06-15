# For mutualism, we are assuming that neutral events remain the same.
# It is the dynamics that we are changing:
#   Producers: Positive Intrinsic Growth Rate
#   Pollinators: Negative Intrinsic Growth Rate
#   Both: intra-specific and intra-guild competition, positive mutualism
# Eq: \dot{x}_i = x_i (r - Sum beta_ij x_j + f(x_i, x))
# where f(x_i, x) = (\sum_j gamma_ij x_j)/(1 + h_i \sum_j gamma_ij x_j)
# using Zhang et al. 2022.
# Thebault and Fontaine 2010 suggest beta is mostly 0 except for intraspecific
# interactions and with some variation in the form of the denominator.
# The largest point of deviation appears to be the numerator:
#    ``c_{i,j} [is the] Maximum rate of... benefit for the interaction''
# for Thebault and Fontaine, but Zhang et al. leave it as 1, which forces scales
# for everything else. T + F set it in [2, 3], which suggests dividing r, beta
# by 2.5.

# Note: These initial guesses do appear to work, but I suspect that they could
# be better. I'm not sure if it is a feature or not to have species unable to
# stay about the extinction threshold without their mutualist (should they have
# been able to invade in the first place? should the calculated fixed point be
# an additional criterion?).

# Data_2023-05-02 had almost all species able to live together.
# Afterwards, I increased the competition by a factor of 10 while keeping
# intraspecific:interspecific competition ratio the same.
# Data_2023-05-04 achieved the same result. I'm going to reduce the
# intra:inter ratio from 10:1 to 2:1
# Data_2023-05-05 achieved the same. Going to take a different tactic.
# Increasing basal reproductive rate, but also see Mutualism2.R.
# Data_2023-05-16 didn't inspire much difference (but waiting for 10^0).
# Changing tactic again: increase number of species (more chances for diffs)
# and look at different dispersal ranges. Reverting the 05-05 change.

Mutualism_species <- function(
  SpeciesTypes,
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

  stopifnot(length(SpeciesTypes) > 0)

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # Assuming Basal/Plant/Producers first.
  r <- c(runif(n = SpeciesTypes[1], min = 0, max = 0.01),
         unlist(lapply(SpeciesTypes[-1], function(n) {
           runif(n, min = -0.2/2.5, max = 0)
         })))

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

Mutualism_saturation <- function(SpeciesTypes, seed = NULL) {
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
          min = 0.1, max = 1),
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
  Pool, seed = NULL
) {

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # We'll start with a ratio of 1:10 for inter:intraspecific competition.
  retmat <- matrix(0, nrow(Pool), nrow(Pool))
  competitors <- outer(Pool$Type, Pool$Type, function(i, j) i == j)
  retmat[competitors] <- -runif(sum(competitors),
                               min = 1 / (2.5 ),#* 10),
                               max = 2 / (2.5 ))#* 10))
  # See initial comments for explanation of 2.5.

  diag(retmat) <- diag(retmat) * 2#10

  # Takes care of block diagonal. Now mutualism off-diagonal blocks.

  retmat[!competitors] <- runif(sum(!competitors),
                                min = 2 / 2.5, max = 3 / 2.5)

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      set.seed(oldSeed)
    }
  }

  return(retmat)
}

PerCapitaDynamics_Mutualistic1 <- function(
  ReproductionRate, InteractionMatrix, NumEnvironments,
  SpeciesTypes, Saturations
) {
  # ReproductionRate = r,
  # InteractionMatrix = Block matrix, Beta blocks on diagonal, Gamma off diag.
  # NumEnvironments = number of patches, as usual
  # Species Types = Species in each Guild (usually length 2: Plant, Pollinator)
  # Subtlety: Mutualism plateaus per-guild in this implementation.
  # Saturations = Matrix of guild-guild saturation values.

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


