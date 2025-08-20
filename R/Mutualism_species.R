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
      .Random.seed <<- oldSeed
    }
  }

  return(Species)
}







