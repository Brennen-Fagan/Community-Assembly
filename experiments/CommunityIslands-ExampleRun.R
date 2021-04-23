# Settings. ####################################################################
islands <- c(1, 2) # Community Numbers

Tolerance <- 1E-6

OuterTimeStepSize <- 5E3
InnerTimeStepSize <- 1E2

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
    seed = s[6]
  )
)

# From LawMorton1996-NumericalTables-Parallel.Rmd
communities <- c(
  "8, 10, 12, 14, 15, 16, 39, 43",
  "8, 12, 14, 15, 16, 38, 39"
)

inits <- c(
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
        func = GeneralisedLotkaVolterra,
        parms = list(a = comMat[com, com],
                     r = pool$ReproductionRate[com],
                     epsilon = epsilon),
        positive = TRUE
      )$y
      print(init)

      if (any(init < epsilon)) {
        # Someone died, reset to random location.
        print("Died")
        anyZeroOrNotSame <- TRUE
        init[init < epsilon] <- runif(n = sum(init < epsilon), min = 2, max = 10000)
      } else if (any(round(init / init_old, 1) != 1)) {
        # Not done, keep going.
        #print(init / init_old)
        #print(round(init / init_old, 1))
        print("Changed")
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

# Reduce to necessary information. #############################################
# Total list of species present.
redCom <- sort(unique(as.numeric(
  unlist(lapply(communities, strsplit, split = ", "))
)))

redPool <- pool[redCom,]
redComMat <- comMat[redCom, redCom]

# Convert to new indices
redComs <- lapply(communities[islands], function(strCom) {
  which(redCom %in% as.numeric(
    unlist(strsplit(strCom, split = ", "))
  ))
})

redPops <- populations[islands]

# Sanity Check: Lengths should match.
stopifnot(all(
  unlist(lapply(redPops, length)) == unlist(lapply(redComs, length))
))

# Evaluation Parameters. #######################################################
# Amount of travel from j to i is d[i,j]
# Amount of gain to i from j is d[i,j]
# Amount of travel from i is d[i,i]
# This way, we can write the change in y from travel as d %*% y
# Characterising this as proportions of a population, but that is assuming a
# normalisation that I do not think is strictly necessary.
# The matrix is sparse and has colsum = 0, diag < 0, offdiag >= 0.
dispersalMatrix <- Matrix::bandSparse(
  n = length(redCom) * length(islands),
  k = length(redCom) * c(1:(length(islands) - 1),
                         -(1:(length(islands) - 1))),
  diagonals = c(
    list(rep(0.0001, length(redCom))), # Island 2 -> Island 1
    list(rep(0.0001, length(redCom))) # Island 1 -> Island 2
  )
)
#   matrix(c(
#   0.9, 0.1,
#   0.01, 0.99
# ), byrow = TRUE, nrow = 2)

stopifnot(all(dispersalMatrix >= 0))

diag(dispersalMatrix) <- -Matrix::colSums(dispersalMatrix)

# Sanity Check: Row sums are 0.
stopifnot(all(Matrix::colSums(dispersalMatrix) == 0))

dynSys <- function(t, y, parms) {
  with(as.list(parms), {
    list(as.numeric(y * (r + a %*% y) + d %*% y))
    # as.numeric since the solver doesn't know what Matrix::Matrices are.
  })
}

parmesan <- list(
  r = rep(redPool$ReproductionRate, length(islands)),
  a = Matrix::bdiag(rep(list(redComMat), length(islands))),
  d = dispersalMatrix,
  epsilon = Tolerance
)

# Technically, a bit of extra work being done here since we already copied the
# populations per island. The result is somewhat more readable though.
abundance_init <- unlist(lapply(
  seq_along(redComs),
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
  com = redComs,
  pop = redPops,
  numPops = length(redCom)
))

abundance <- deSolve::ode(
  abundance_init,
  times = ts,
  func = dynSys,
  parms = parmesan,
  events = list(func = function(t, y, parms) {
    y[y < parms$epsilon] <- 0
    y
  }, time = ts)
)

abundance_island <- lapply(
  1:length(islands), function(i, x) {
    x[, c(1, (1:length(redCom)) + 1 + (i - 1) * length(redCom))]
  },
  x = abundance
)

LawMorton1996_PlotAbundance(abundance_island[[1]]) -> figure1
LawMorton1996_PlotAbundance(abundance_island[[2]]) -> figure2


print(figure1 + ggplot2::ggtitle(paste("Init:", toString(redComs[[islands[1]]]))))
print(figure2 + ggplot2::ggtitle(paste("Init:", toString(redComs[[islands[2]]]))))

print("Island 1 Species:")
print(SpeciesPresent(abundance_island[[1]]))
print("Island 2 Species:")
print(SpeciesPresent(abundance_island[[2]]))
