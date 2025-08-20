LawMorton1996_species <- function(
    Basal,
    Consumer,
    Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
    LogBodySize = c(-2, -1, -1, 0), # c(-2, -1) for Basal, c(-1, 0) for Consumer
    seed = NULL
) {
  stopifnot(Basal > 0)
  stopifnot(Consumer >= 0)
  stopifnot(length(Parameters) == 6)

  if (Parameters[4] >= 1) {
    warning("Unrealistic energy consumption efficiency k4: k4 >= 1.")
  }
  stopifnot(Parameters[6] > 0)

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  # Assign each species in the pool a body size si.
  # Do so by drawing from a uniform distribution and exponentiating.
  # For all i, j, if si < sj then i may be eaten by j but not vice versa.

  Species <- data.frame(
    ID = 1:(Basal + Consumer),
    Type = c(rep("Basal", Basal), rep("Consumer", Consumer)),
    Size = 10^(c(runif(Basal, min = LogBodySize[1], max = LogBodySize[2]),
                 runif(Consumer, min = LogBodySize[3], max = LogBodySize[4]))),
    ReproductionRate = 0
  )

  # For species i,
  # if i is basal, set p to 10^(-1 - 0.25 log10 si),
  # otherwise set p to -0.1.
  # Draw ri from a truncated normal distribution, mean p, std. p k6.
  # (Set the sign of ri to that of p).

  Species$ReproductionRate[1:Basal] <- unlist(lapply(
    Species$Size[1:Basal], function(si, k6) {
      p <- 10 ^ (-1 - 0.25 * log10(si))
      return(sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * k6))
    },
    k6 = Parameters[6]
  ))

  if (Consumer > 0) {
    Species$ReproductionRate[(Basal + 1) : (Basal + Consumer)] <- unlist(lapply(
      Species$Size[(Basal + 1) : (Basal + Consumer)], function(si, k6) {
        p <- -0.1
        return(sign(p) * rtruncnorm(0, Inf, abs(p), abs(p) * k6))
      },
      k6 = Parameters[6]
    ))
  }

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }

  Species
}







