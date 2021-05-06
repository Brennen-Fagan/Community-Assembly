# Reproduce CommunityIslands-ExampleRun2.R with
# code in CommunityIslands.R.

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
  Basal = b[2],
  Consumer = c[2],
  Parameters = parameters,
  LogBodySize = logBodySize,
  seed = s[2]
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
  seed = s[2]
)
)

# From LawMorton1996-NumericalTables-Parallel.Rmd
communities <- c(
  "",
  "2, 4, 6, 12, 29",
  "2, 4, 6, 13, 29"
)

inits <- c(
  "",
  "10, 10, 10, 10, 10",
  "10, 10, 10, 10, 10"
)

stopifnot(length(communities) == length(inits))

# > runif(1) * 1E8
# [1] 20663521
set.seed(20663521)

populationsTest <- lapply(
  seq_along(communities),
  function(i, Pool, InteractionMatrix, Communities, Populations) {
    RMTRCode2::FindSteadyStateFromEstimate(
      Pool, InteractionMatrix, Communities[i], Populations[i]
    )
  },
  Pool = pool,
  InteractionMatrix = comMat,
  Communities = communities,
  Populations = inits
)

# stopifnot(all(unlist(lapply(1:3, function(i) {all(populations[[i]] == populationsTest[[i]])}))))

productivityTest <- lapply(
  seq_along(communities),
  function(i, Pool, InteractionMatrix, Communities, Populations) {
    Productivity(Pool, InteractionMatrix, Communities[i], Populations[[i]])
  },
  Pool = pool,
  InteractionMatrix = comMat,
  Communities = communities,
  Populations = populationsTest
)

# Interestingly, there is an implementation difference here.

abundance1Test <- RMTRCode2::IslandDynamics(
  Pool = pool,
  InteractionMatrix = comMat,
  Communities = list(
    communities[2],
    communities[3]
  ),
  Populations = list(
    populationsTest[[2]],
    populationsTest[[3]]
  ),
  DispersalPool = 0.0001, # rep(0.0001, nrow(pool)),
  DispersalIsland = matrix(0, nrow = 2, ncol = 2),
  Tolerance = 1E-6
)
# stopifnot(all(abundance1 == abundance1Test))

abundance2Test <- RMTRCode2::IslandDynamics(
  Pool = pool,
  InteractionMatrix = comMat,
  Communities = list(
    communities[2],
    communities[3]
  ),
  Populations = list(
    populationsTest[[2]],
    populationsTest[[3]]
  ),
  DispersalPool = 0.0001, # rep(0.0001, nrow(pool)),
  DispersalIsland = matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2),
  Tolerance = 1E-6
)
# stopifnot(all(abundance2 == abundance2Test))

abundance3Test <- RMTRCode2::IslandDynamics(
  Pool = pool,
  InteractionMatrix = comMat,
  Communities = list(
    communities[2],
    communities[1],
    communities[3]
  ),
  Populations = list(
    populationsTest[[2]],
    populationsTest[[1]],
    populationsTest[[3]]
  ),
  DispersalPool = 0.0001, # rep(0.0001, nrow(pool)),
  DispersalIsland = matrix(c(
    0, 1, 0, # Island 2 -> 1
    1, 0, 1, # Island 1 -> 2, Island 3 -> 2
    0, 1, 0  # Island 2 -> 3
  ), nrow = 3, ncol = 3, byrow = TRUE),
  Tolerance = 1E-6
)
# stopifnot(all(abundance3 == abundance3Test))
