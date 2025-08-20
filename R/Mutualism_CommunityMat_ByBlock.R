

Mutualism_CommunityMat_ByBlock <- function(
  Pool,
  MinimumGuildMatrix, # Entry i,j is the min/max effect of species in guild j on
  MaximumGuildMatrix, # species in guild i. Allows mutual., pred., comp., etc.
  IntraspeciesMultiplier = 2, # rep(Mult, times = table(Pool$Type))
  seed = NULL
) {
  # Requirements:
  stopifnot(length(unique(Pool$Type)) == nrow(MinimumGuildMatrix),
            nrow(MinimumGuildMatrix) == ncol(MinimumGuildMatrix),
            nrow(MinimumGuildMatrix) == ncol(MaximumGuildMatrix),
            nrow(MaximumGuildMatrix) == ncol(MaximumGuildMatrix),
            length(IntraspeciesMultiplier) == 1 ||
              length(IntraspeciesMultiplier) == nrow(Pool$Types))
  
  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }
  
  # Intraguild interactions
  retmat <- matrix(0, nrow(Pool), nrow(Pool))
  dict <- unique(Pool$Type)
  
  for (j in 1:ncol(MinimumGuildMatrix)) {
    for (i in 1:nrow(MinimumGuildMatrix)) {
      targets <- outer(Pool$Type, Pool$Type,
                       function(x, y) x == dict[i] & y == dict[j])
      retmat[targets] <- runif(sum(targets),
                               min = MinimumGuildMatrix[i, j],
                               max = MaximumGuildMatrix[i, j])
    }}
  
  # Intraspecies (subset of Intraguild) interactions.
  diag(retmat) <- diag(retmat) * IntraspeciesMultiplier
  
  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }
  
  return(retmat)
}

# test:
# pool <- Mutualism_species(c(5, 10))
# stopifnot(
# # Won't be exact same, former goes diagonal first...
# all(sign(Mutualism_CommunityMat(pool, seed = 1)) ==
#       sign(Mutualism_CommunityMat_ByBlock(
#         pool,
#         MinimumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
#                                     c(-2/2.5, 2/2.5, 2/2.5, -2/2.5)),
#         MaximumGuildMatrix = matrix(byrow = TRUE, nrow = 2, ncol = 2,
#                                     c(-1/2.5, 3/2.5, 3/2.5, -1/2.5)),
#         seed = 1
#       )))
# )
