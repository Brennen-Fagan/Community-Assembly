# Relies on Chunk 19 at most of LM1996-NumPoolCom-Questions-2021-05.Rmd
# We will take as an example run
# Combination 5 Dataset 2,
# Communities 7 and 12
# which appear to be mutually invadable.
exampleDatasetID <- 2
exampleCombnNum <- 5
exampleCommunities <- communitiesAll %>% dplyr::filter(
  TotalID == paste(exampleCombnNum, exampleDatasetID)
)
exampleCommunityNums <- c(7, 12)
examplePool <- pools[[exampleDatasetID]][[exampleCombnNum]]
exampleMat <- mats[[exampleDatasetID]][[exampleCombnNum]]

stopifnot(
  exampleCommunities$Communities[exampleCommunityNums[1]] !=
    exampleCommunities$Communities[exampleCommunityNums[2]]
    )

exampleDisPool <- 0.0001 # Uniform over all species
exampleDisIslandList <- list(
  # Remember that they are column ordered.
  # Should also be row-to col-from.
  # Directed 1 -> 2
  matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2),
  # Directed 1 <- 2
  matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2),
  # Undirected
  matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
)
Tolerance <- 1E-6

examplepreprocessed1 <- RMTRCode2::IslandPreprocess(
  Pool = examplePool,
  InteractionMatrix = exampleMat,
  Communities = exampleCommunities$Communities[exampleCommunityNums],
  Populations = exampleCommunities$CommunityAbund[exampleCommunityNums],
  DispersalPool = exampleDisPool,
  DispersalIsland = exampleDisIslandList[[1]],
  Tolerance = Tolerance
)

examplepreprocessed2 <- RMTRCode2::IslandPreprocess(
  Pool = examplePool,
  InteractionMatrix = exampleMat,
  Communities = exampleCommunities$Communities[exampleCommunityNums],
  Populations = exampleCommunities$CommunityAbund[exampleCommunityNums],
  DispersalPool = exampleDisPool,
  DispersalIsland = exampleDisIslandList[[2]],
  Tolerance = Tolerance
)

examplepreprocessed3 <- RMTRCode2::IslandPreprocess(
  Pool = examplePool,
  InteractionMatrix = exampleMat,
  Communities = exampleCommunities$Communities[exampleCommunityNums],
  Populations = exampleCommunities$CommunityAbund[exampleCommunityNums],
  DispersalPool = exampleDisPool,
  DispersalIsland = exampleDisIslandList[[3]],
  Tolerance = Tolerance
)

ArrivalEvents <- 100

# Need to extract the entire community
### This is in $redCom or $redComs.
# Need to create one copy per island link.
### Done in next step.
# Need to shuffle the copies.
invaders1 <- lapply(
  exampleDisIslandList[[1]],
  function(i, len, events) {
    if (i > 0) {
      replicate(
        n = ceiling(events / len),
        sample.int(len, replace = FALSE)
      )[
        1:events
      ]
    } else {
      NULL
    }
  },
  len = length(examplepreprocessed1$redCom),
  events = ArrivalEvents
)
invaders2 <- lapply(
  exampleDisIslandList[[2]],
  function(i, len, events) {
    if (i > 0) {
      replicate(
        n = ceiling(events / len),
        sample.int(len, replace = FALSE)
      )[
        1:events
      ]
    } else {
      NULL
    }
  },
  len = length(examplepreprocessed2$redCom),
  events = ArrivalEvents
)
invaders3 <- lapply(
  exampleDisIslandList[[3]],
  function(i, len, events) {
    if (i > 0) {
      replicate(
        n = ceiling(events / len),
        sample.int(len, replace = FALSE)
      )[
        1:events
      ]
    } else {
      NULL
    }
  },
  len = length(examplepreprocessed3$redCom),
  events = ArrivalEvents
)

# Iterate across the copies in parallel.
for (event in 1:ArrivalEvents) {
  # The elements in the copies tell us who is
  # attempting to invade.
  for (island in 1:nrow(islands)) {

  }
}
# We then check to see if invasion begins.
# Since copies are per link, we do not need to
# make sure that the link exists.
# We instead need to check to make sure that the
# species can leave the host island (sufficient
# abundance) and whether the victim island is
# invadable by this type.
# If either condition fails, the invasion does
# not take place.
# All systems (with or w/o invasions!) are then
# progressed in parallel using GLV dynamics.
# Once the invasions are done, check whether the
# system is invadable by any species not in the
# system and, if uninvadable, check the
# populations against the previous populations
# to see if the system is in steady-state.
