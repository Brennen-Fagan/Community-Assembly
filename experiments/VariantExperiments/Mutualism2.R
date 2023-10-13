# Following Mutualism.R's failure to separate medium from no dispersal cases,
# I'm interpreting the failure as insufficient control and competition.
# Let's set the abundance of producers to 10 without pollinators or competitors.
# Set abundance of producers to 100 with a single pollinator.
# Set abundance of producers to 1 in presence of a competitor.
# Set abundance of producers to 10 with pollinator, competitor and c's p.

# 1st: To accomplish the first, set reproduction to 10 * intraspecific.
# 3rd: with one other, it's pretty clear interspecific > intraspecific.
#      (By a factor of 9). If we accept that interspecific < intraspecific,
#      by a factor of around 10, then you need about 90 species to trump it,
#      which might be the source of our problems.
#      Instead, set interspecific to be a bit less than intraspecific. 80%?
# 4th: With target abund = 10, competitor = 10, pollinators = 1, intra = 1/10 r,
#      and inter (fudging number of species!) = 9/10 r (r = rep rate), rearrange
#      to get 9 r / (1 - 9 r h) for mutualism. Bottom needs to be > 0, so
#      saturation < 1 / (9 r), say 1 / (10 r) at most. The sum due to mutualism
#      then should be 9 r / (1/10) = 90r.

targetRepRate <- 0.1

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
  r <- c(runif(n = SpeciesTypes[1],
               min = targetRepRate / 2,
               max = targetRepRate * 2),
         unlist(lapply(SpeciesTypes[-1], function(n) {
           -runif(n,
                 min = targetRepRate / 2,
                 max = targetRepRate * 2)
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
          min = min(1 / (2 * targetRepRate) * 1/10, 0.1),
          max = min(2 / targetRepRate * 1/10, 1)) , # See comment #4 above.
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
                                min = targetRepRate / 2,
                                max = targetRepRate * 2) / 10 * 0.8
  # See initial comments (#3) for explanation.

  diag(retmat) <- -Pool$ReproductionRate / 10 # See init. comments (#1).

  # Takes care of block diagonal. Now mutualism off-diagonal blocks.

  retmat[!competitors] <- runif(sum(!competitors),
                                min = 9 * targetRepRate,
                                max = 90 * targetRepRate) # see comment (#4).

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
