# Recreate the communities. ####################################################
set.seed(3680180)

basal2 <- c(5, 10, 15)
consumer2 <- c(20, 40, 60)
logBodySize = c(-2, -1, -1, 0)
parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)

seedsPrep2 <- runif(2 * length(basal2) * length(consumer2)) * 1E8

pool <- with(list(
  b = rep(basal2, times = length(consumer2)),
  c = rep(consumer2, each = length(basal2)),
  s = seedsPrep2[1:(length(basal2) * length(consumer2))],
  ),
  LawMorton1996_species(
    Basal = b[6],
    Consumer = c[6],
    Parameters = parameters,
    LogBodySize = logBodySize,
    seed = s[6]
  )
)

comMat <- with(list(
  s = seedsPrep2[
    (length(basal2) * length(consumer2) + 1):(
      2 * length(basal2) * length(consumer2))
  ]
  ),
  LawMorton1996_CommunityMat(
    Pool = pool,
    Parameters = parameters,
    seed = s
  )
)

# From LawMorton1996-NumericalTables-Parallel.Rmd
communities <- c(
  "8, 10, 12, 14, 15, 16, 39, 43",
  "8, 12, 14, 15, 16, 38, 39"
)

# > runif(1) * 1E8
# [1] 20649322
set.seed(20649322)

populations <- lapply(communities, function(strCom) {
  # Convert to indices.
  com <- as.numeric(
    unlist(strsplit(strCom, split = ", "))
  )

  # Attempt to calculate the (assumed) stable steady state.
  # Make a guess.
  init <- runif(n = length(com), min = 1, max = 1000)

  anyZeroOrNotSame <- TRUE # Set T to make at least one look.

  # Run and check to see if anyone dies (bad) or changes (not at steady).
  while (anyZeroOrNotSame) {
    init_old <- init
    init <- LawMorton1996_NumIntegration(
      A = comMat[com, com],
      R = pool$ReproductionRate[com],
      X = init,
      Tolerance = 1E-10
    )

    if (any(init == 0)) {
      # Someone died, reset to random location.
      anyZeroOrNotSame <- TRUE
      init <- runif(n = length(com), min = 1, max = 1000)
    } else if (any(round(init / init_old) != 1)) {
      # Not done, keep going.
      anyZeroOrNotSame <- TRUE
    } else {
      anyZeroOrNotSame <- FALSE
    }
  }

  return(init)
})

# Reduce to necessary information. #############################################
# Total list of species present.
redCom <- sort(unique(as.numeric(
  unlist(lapply(communities, strsplit, split = ", "))
)))

redPool <- pool[redCom]
redComMat <- comMat[redCom, redCom]

# Convert to new indices
redComs <- lapply(communities, function(strCom) {
  which(redCom %in% as.numeric(
    unlist(strsplit(strCom, split = ", "))
  ))
})

# Sanity Check: Lengths should match.
stopifnot(all(
  unlist(lapply(populations, length)) == unlist(lapply(redComs, length))
))

# Evaluation Parameters. #######################################################
# Amount of travel from island i to island j is [i, j].
# Amount of gain from island i to island j is [j, i].
# Characterising this as proportions of a population, but that is assuming a
# normalisation that I do not think is strictly necessary.
dispersalMatrix <- matrix(c(
  0.9, 0.1,
  0.1, 0.9
), byrow = TRUE)

# Sanity Check: Row sums are 1.
stopifnot(all(rowSums(dispersalMatrix) == 1))

dispersalMatrixConverted <- dispersalMatrix
diag(dispersalMatrixConverted) <- 0
dispersalMatrixConverted <- -dispersalMatrixConverted

dynSys <- function(t, y, parms) {
  with(as.list(parms), {
    # indexes i, j, patches k, l.
    # dy_(i,k)/dt = reproduction     : y_(i,k) * r
    #             + interaction      : Sum_j a_ji y_(j,k) y_(i,k)
    #             + dispersal gain   : Sum_{l != k} d_lk y_i
    #             + dispersal loss   :
    list(y * (r + a %*% y))
  })
}

parmesan <- list(
  r = redPool$ReproductionRate,
  a = redComMat
)

# Repeat to match number of entries we are working with.
parmesan <- lapply(parmesan, rep, times = nrow(dispersalMatrix))

abundance_init <- unlist(lapply(
  seq_along(redComs),
  function(i, com, pop, numPops) {
    # Interlace 0's with population values
    k <- 1
    retval <- rep(0, numPops)
    for (j in 1:numPops) {
      if (j %in% com) {
        retval[j] <- pop[k]
        k <- k + 1
      }
    }
    return(retval)
  },
  com = redComs,
  pop = populations,
  numPops = length(redCom)
))
