LawMorton1996_CheckUninvadable <- function(
  AbundanceRow,
  Pool,
  CommunityMatrix,
  Threshold = 1E-4,
  EliminateSmallPopulations = FALSE
) {
  # It is uninvadable by any other species from the pool because the per
  # capita rate of increase of each species absent from the community is
  # negative at the equilibrium point.
  # The per capita rate of increase is the function f_i = r_i + Sum_j(aij xj).
  Abundance <- AbundanceRow[-1]

  if (EliminateSmallPopulations) {
    Abundance[is.na(Abundance) | Abundance <= Threshold] <- 0
  } else {
    Abundance[is.na(Abundance)] <- 0
  }

  notPresent <- !(!is.na(Abundance) & Abundance > Threshold)
  all(Pool$ReproductionRate[notPresent] +
        CommunityMatrix[notPresent, ] %*% Abundance < 0)
}

LawMorton1996_CheckPermanence <- function(
  Pool,
  CommunityMatrix
) {
  # We assume that the pool and community matrix are only for the current community.
  # We extract all subcommunities that form the boundary (2^nrow(Pool) - 1).
  # We then evaluate each one's equilibrium point.
  # Using each equilibrium point, we then perform the linear program.
  # Minimize z
  # subject to, for j = 1 : 2^nrow(Pool), sum(f_i (eqpoint_j) h_i) + z > 0
  # and (using h_i subject to) 0 < h_i <= 1 for i = 1 : nrow(Pool).

  if (is.null(nrow(Pool)) || nrow(Pool) == 0) {
    return(TRUE)
  }

  # Calculate powerset. Skip empty set.
  sets <- unlist(
    lapply(1:(nrow(Pool)), combn,
           x = 1:nrow(Pool), simplify = FALSE),
    recursive = FALSE
  )

  # Find equilibria of each set and evaluate the per capita growth rates.
  equilibria <- lapply(sets, function(set, pool, mat) {
    p <- pool[set,]
    m <- mat[set, set]

    rootSolve::steady(
      y = rep(1, length(set)),
      func = GeneralisedLotkaVolterra,
      parms = list(a = m, r = p$ReproductionRate),
      positive = TRUE
    )
  },
  pool = Pool, mat = CommunityMatrix)

}
