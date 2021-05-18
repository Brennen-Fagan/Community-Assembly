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

examplepreprocessed1 <- RMTRCode2::IslandPreprocess(
  Pool = examplePool,
  InteractionMatrix = exampleMat,
  Communities = exampleCommunities$Communities[exampleCommunityNums],
  Populations = exampleCommunities$CommunityAbund[exampleCommunityNums],
  DispersalPool = exampleDisPool,
  DispersalIsland = exampleDisIslandList[[1]]
)

examplepreprocessed2 <- RMTRCode2::IslandPreprocess(
  Pool = examplePool,
  InteractionMatrix = exampleMat,
  Communities = exampleCommunities$Communities[exampleCommunityNums],
  Populations = exampleCommunities$CommunityAbund[exampleCommunityNums],
  DispersalPool = exampleDisPool,
  DispersalIsland = exampleDisIslandList[[2]]
)

examplepreprocessed3 <- RMTRCode2::IslandPreprocess(
  Pool = examplePool,
  InteractionMatrix = exampleMat,
  Communities = exampleCommunities$Communities[exampleCommunityNums],
  Populations = exampleCommunities$CommunityAbund[exampleCommunityNums],
  DispersalPool = exampleDisPool,
  DispersalIsland = exampleDisIslandList[[3]]
)
