# Settings. ####################################################################
islands1 <- c(2, 3) # Community Numbers, No interaction
islands2 <- c(2, 3) # Community Numbers, interaction
islands3 <- c(2, 1, 3) # Community Numbers, buffer/empty island in middle

Tolerance <- 1E-6

OuterTimeStepSize <- 2E4
InnerTimeStepSize <- 5E2

ts <- seq(from = 0,
          to = OuterTimeStepSize,
          by = InnerTimeStepSize)

# Recreate the communities. ####################################################
#library(RMTRCode2)

set.seed(3680180)

basal2 <- c(5, 10, 15)
consumer2 <- c(20, 40, 60)
logBodySize = c(-2, -1, -1, 0)
parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)

seedsPrep2 <- runif(2 * length(basal2) * length(consumer2)) * 1E8

pool <- with(list(
  b = rep(basal2, times = length(consumer2)),
  c = rep(consumer2, each = length(basal2)),
  s = seedsPrep2[1:(length(basal2) * length(consumer2))]
  ),
  RMTRCode2::LawMorton1996_species(
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
  RMTRCode2::LawMorton1996_CommunityMat(
    Pool = pool,
    Parameters = parameters,
    seed = s[6]
  )
)

# From LawMorton1996-NumericalTables-Parallel.Rmd
communities <- c(
  "",
  "8, 10, 12, 14, 15, 16, 39, 43",
  "8, 12, 14, 15, 16, 38, 39"
)

inits <- c(
  "",
  "21, 32, 80, 818, 121, 18, 12, 20",
  "83, 138, 938, 52, 5, 14, 26"
)

stopifnot(length(communities) == length(inits))

# > runif(1) * 1E8
# [1] 20649322
set.seed(20649322)

populations <- lapply(
  seq_along(communities),
  function(i, communities, inits) {
    # Convert to indices.
    com <- as.numeric(
      unlist(strsplit(communities[i], split = ", "))
    )
    init <- as.numeric(
      unlist(strsplit(inits[i], split = ", "))
    )

    anyZeroOrNotSame <- TRUE # Set T to make at least one look.
    epsilon <- 1E-1

    # Run and check to see if anyone dies (bad) or changes (not at steady).
    while (anyZeroOrNotSame) {
      init_old <- init
      # Initial attempt: random search plus dynamics.
      # init <- LawMorton1996_NumIntegration(
      #   A = comMat[com, com],
      #   R = pool$ReproductionRate[com],
      #   X = init,
      #   OuterTimeStepSize = 1000,
      #   InnerTimeStepSize = 1,
      #   Tolerance = 1E-10
      # )
      # init <- init[nrow(init),]
      # init <- init[-1]
      #
      # Applying a more specialised tool.
      init <- rootSolve::steady(
        y = init,
        func = RMTRCode2::GeneralisedLotkaVolterra,
        parms = list(a = comMat[com, com],
                     r = pool$ReproductionRate[com],
                     epsilon = epsilon),
        positive = TRUE
      )$y
      # print(init)

      if (any(init < epsilon)) {
        # Someone died, reset to random location.
        # print("Died")
        anyZeroOrNotSame <- TRUE
        init[init < epsilon] <- runif(n = sum(init < epsilon), min = 2, max = 10000)
      } else if (any(round(init / init_old, 1) != 1)) {
        # Not done, keep going.
        #print(init / init_old)
        #print(round(init / init_old, 1))
        # print("Changed")
        anyZeroOrNotSame <- TRUE
      } else {
        anyZeroOrNotSame <- FALSE
      }
    }

    return(init)
  },
  communities = communities,
  inits = inits
)

productivity <- lapply(
  seq_along(communities),
  function(i, communities, pops) {
    com <- as.numeric(
      unlist(strsplit(communities[i], split = ", "))
    )

    if (length(com) == 0) {
      return(NA)
    }

    pop <- pops[[i]]

    comMatPos <- comMat[com, com]; comMatPos[comMatPos < 0] <- 0
    poolRepPos <- pool$ReproductionRate[com]; poolRepPos[poolRepPos < 0] <- 0

    parmesanPos <- list(
      a = comMatPos,
      r = poolRepPos
    )

    return(
      sum(RMTRCode2::GeneralisedLotkaVolterra(0, pop, parmesanPos)[[1]]
          * pop / sum(pop))
    )

  },
  communities = communities,
  pops = populations
)

# Reduce to necessary information. #############################################
# Total list of species present.
redCom <- sort(unique(as.numeric(
  unlist(lapply(communities, strsplit, split = ", "))
)))

redPool <- pool[redCom,]
redComMat <- comMat[redCom, redCom]

# Convert to new indices
redComs1 <- lapply(communities[islands1], function(strCom) {
  which(redCom %in% as.numeric(
    unlist(strsplit(strCom, split = ", "))
  ))
})
redComs2 <- lapply(communities[islands2], function(strCom) {
  which(redCom %in% as.numeric(
    unlist(strsplit(strCom, split = ", "))
  ))
})
redComs3 <- lapply(communities[islands3], function(strCom) {
  which(redCom %in% as.numeric(
    unlist(strsplit(strCom, split = ", "))
  ))
})

redPops1 <- populations[islands1]
redPops2 <- populations[islands2]
redPops3 <- populations[islands3]

# Sanity Check: Lengths should match.
stopifnot(all(
  unlist(lapply(redPops1, length)) == unlist(lapply(redComs1, length))
))
stopifnot(all(
  unlist(lapply(redPops2, length)) == unlist(lapply(redComs2, length))
))
stopifnot(all(
  unlist(lapply(redPops3, length)) == unlist(lapply(redComs3, length))
))

# Evaluation Parameters. #######################################################
# Amount of travel from j to i is d[i,j]
# Amount of gain to i from j is d[i,j]
# Amount of travel from i is d[i,i]
# This way, we can write the change in y from travel as d %*% y
# Characterising this as proportions of a population, but that is assuming a
# normalisation that I do not think is strictly necessary.
# The matrix is sparse and has colsum = 0, diag < 0, offdiag >= 0.
dispersalMatrix1 <- Matrix::bandSparse(
  n = length(redCom) * length(islands2),
  k = length(redCom) * c(1:(length(islands2) - 1),
                         -(1:(length(islands2) - 1))),
  diagonals = c(
    list(rep(0, length(redCom))), # Island 2 -> Island 1
    list(rep(0, length(redCom)))  # Island 1 -> Island 2
  )
)
dispersalMatrix2 <- Matrix::bandSparse(
  n = length(redCom) * length(islands2),
  k = length(redCom) * c(1:(length(islands2) - 1),
                         -(1:(length(islands2) - 1))),
  diagonals = c(
    list(rep(0.0001, length(redCom))), # Island 2 -> Island 1
    list(rep(0.0001, length(redCom)))  # Island 1 -> Island 2
  )
)
dispersalMatrix3 <- Matrix::bandSparse(
  n = length(redCom) * length(islands3),
  k = length(redCom) * c(1:(length(islands3) - 1),
                         -(1:(length(islands3) - 1))),
  diagonals = c(
    list(c(rep(0.0001, length(redCom)),   # Island 2 -> Island 1
           rep(0.0001, length(redCom)))), # Island 3 -> Island 2
    list(rep(0, length(redCom))),         # Island 3 -> Island 1
    list(c(rep(0.0001, length(redCom)),   # Island 1 -> Island 2
           rep(0.0001, length(redCom)))), # Island 2 -> Island 3
    list(rep(0, length(redCom)))          # Island 1 -> Island 3
  )
)
#   matrix(c(
#   0.9, 0.1,
#   0.01, 0.99
# ), byrow = TRUE, nrow = 2)

stopifnot(all(dispersalMatrix1 >= 0))
stopifnot(all(dispersalMatrix2 >= 0))
stopifnot(all(dispersalMatrix3 >= 0))

diag(dispersalMatrix1) <- -Matrix::colSums(dispersalMatrix1)
diag(dispersalMatrix2) <- -Matrix::colSums(dispersalMatrix2)
diag(dispersalMatrix3) <- -Matrix::colSums(dispersalMatrix3)

# Sanity Check: Row sums are 0.
stopifnot(all(Matrix::colSums(dispersalMatrix1) == 0))
stopifnot(all(Matrix::colSums(dispersalMatrix2) == 0))
stopifnot(all(Matrix::colSums(dispersalMatrix3) == 0))

dynSys <- function(t, y, parms) {
  with(as.list(parms), {
    list(as.numeric(y * (r + a %*% y) + d %*% y))
    # as.numeric since the solver doesn't know what Matrix::Matrices are.
  })
}

parmesan1 <- list(
  r = rep(redPool$ReproductionRate, length(islands1)),
  a = Matrix::bdiag(rep(list(redComMat), length(islands1))),
  d = dispersalMatrix1,
  epsilon = Tolerance
)
parmesan2 <- list(
  r = rep(redPool$ReproductionRate, length(islands2)),
  a = Matrix::bdiag(rep(list(redComMat), length(islands2))),
  d = dispersalMatrix2,
  epsilon = Tolerance
)
parmesan3 <- list(
  r = rep(redPool$ReproductionRate, length(islands3)),
  a = Matrix::bdiag(rep(list(redComMat), length(islands3))),
  d = dispersalMatrix3,
  epsilon = Tolerance
)

# Technically, a bit of extra work being done here since we already copied the
# populations per island. The result is somewhat more readable though.
abundance_init1 <- unlist(lapply(
  seq_along(redComs1),
  function(i, com, pop, numPops) {
    # Interlace 0's with population values
    k <- 1
    retval <- rep(0, numPops)
    for (j in 1:numPops) {
      if (j %in% com[[i]]) {
        retval[j] <- pop[[i]][k]
        k <- k + 1
      }
    }
    return(retval)
  },
  com = redComs1,
  pop = redPops1,
  numPops = length(redCom)
))
abundance_init2 <- unlist(lapply(
  seq_along(redComs2),
  function(i, com, pop, numPops) {
    # Interlace 0's with population values
    k <- 1
    retval <- rep(0, numPops)
    for (j in 1:numPops) {
      if (j %in% com[[i]]) {
        retval[j] <- pop[[i]][k]
        k <- k + 1
      }
    }
    return(retval)
  },
  com = redComs2,
  pop = redPops2,
  numPops = length(redCom)
))
abundance_init3 <- unlist(lapply(
  seq_along(redComs3),
  function(i, com, pop, numPops) {
    # Interlace 0's with population values
    k <- 1
    retval <- rep(0, numPops)
    for (j in 1:numPops) {
      if (j %in% com[[i]]) {
        retval[j] <- pop[[i]][k]
        k <- k + 1
      }
    }
    return(retval)
  },
  com = redComs3,
  pop = redPops3,
  numPops = length(redCom)
))

abundance1 <- deSolve::ode(
  abundance_init1,
  times = ts,
  func = dynSys,
  parms = parmesan1,
  events = list(func = function(t, y, parms) {
    y[y < parms$epsilon] <- 0
    y
  }, time = ts)
)
abundance2 <- deSolve::ode(
  abundance_init2,
  times = ts,
  func = dynSys,
  parms = parmesan2,
  events = list(func = function(t, y, parms) {
    y[y < parms$epsilon] <- 0
    y
  }, time = ts)
)
abundance3 <- deSolve::ode(
  abundance_init3,
  times = ts,
  func = dynSys,
  parms = parmesan3,
  events = list(func = function(t, y, parms) {
    y[y < parms$epsilon] <- 0
    y
  }, time = ts)
)

abundance1_island <- lapply(
  1:length(islands1), function(i, x) {
    x[, c(1, (1:length(redCom)) + 1 + (i - 1) * length(redCom))]
  },
  x = abundance1
)
abundance2_island <- lapply(
  1:length(islands2), function(i, x) {
    x[, c(1, (1:length(redCom)) + 1 + (i - 1) * length(redCom))]
  },
  x = abundance2
)
abundance3_island <- lapply(
  1:length(islands3), function(i, x) {
    x[, c(1, (1:length(redCom)) + 1 + (i - 1) * length(redCom))]
  },
  x = abundance3
)

RMTRCode2::LawMorton1996_PlotAbundance(abundance1_island[[1]]) -> figure1_1
RMTRCode2::LawMorton1996_PlotAbundance(abundance1_island[[2]]) -> figure1_2

RMTRCode2::LawMorton1996_PlotAbundance(abundance2_island[[1]]) -> figure2_1
RMTRCode2::LawMorton1996_PlotAbundance(abundance2_island[[2]]) -> figure2_2

RMTRCode2::LawMorton1996_PlotAbundance(abundance3_island[[1]]) -> figure3_1
RMTRCode2::LawMorton1996_PlotAbundance(abundance3_island[[2]]) -> figure3_2
RMTRCode2::LawMorton1996_PlotAbundance(abundance3_island[[3]]) -> figure3_3


print(figure1_1 + ggplot2::ggtitle(paste("1 Island, Init:", toString(redComs1[[1]]))))
print(figure1_2 + ggplot2::ggtitle(paste("1 Island, Init:", toString(redComs1[[2]]))))
print(figure2_1 + ggplot2::ggtitle(paste("2 Island, Init:", toString(redComs2[[1]]))))
print(figure2_2 + ggplot2::ggtitle(paste("2 Island, Init:", toString(redComs2[[2]]))))
print(figure3_1 + ggplot2::ggtitle(paste("3 Island, Init:", toString(redComs3[[1]]))))
print(figure3_2 + ggplot2::ggtitle(paste("3 Island, Init:", toString(redComs3[[2]]))))
print(figure3_3 + ggplot2::ggtitle(paste("3 Island, Init:", toString(redComs3[[3]]))))

print("1 Island 1 Species:")
print(RMTRCode2::SpeciesPresent(abundance1_island[[1]]))
print("1 Island 2 Species:")
print(RMTRCode2::SpeciesPresent(abundance1_island[[2]]))

print("2 Island 1 Species:")
print(RMTRCode2::SpeciesPresent(abundance2_island[[1]]))
print("2 Island 2 Species:")
print(RMTRCode2::SpeciesPresent(abundance2_island[[2]]))

print("3 Island 1 Species:")
print(RMTRCode2::SpeciesPresent(abundance3_island[[1]]))
print("3 Island 2 Species:")
print(RMTRCode2::SpeciesPresent(abundance3_island[[2]]))
print("3 Island 3 Species:")
print(RMTRCode2::SpeciesPresent(abundance3_island[[3]]))
