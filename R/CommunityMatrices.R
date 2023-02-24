# Foreach operators
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`
`%:%` <- foreach::`%:%`

SpeciesPresent <- function(Abundance, Threshold = 1E-4) {
  # resultsNotNAs <- apply(
  #   Abundance[, -1],
  #   MARGIN = 1,
  #   FUN = function(row, epsilon) {
  #     which(!is.na(row) & row > epsilon)
  #   }, epsilon = Threshold)
  resultsNotNAs <- foreach::foreach(
    row = iterators::iter(Abundance[, -1], by = 'row')
  ) %do% {which(!is.na(row) & row > Threshold)}

  tempSet <- resultsNotNAs[[1]]
  sets <- data.frame(
    TimeStart = Abundance[1, 1],
    TimeStop = NA,
    Present = I(list(tempSet))
  )
  setsIndex <- 1
  for (i in 2:length(resultsNotNAs)) {
    if (!setequal(tempSet, resultsNotNAs[[i]])) {
      sets$TimeStop[setsIndex] <- Abundance[i - 1, 1]
      setsIndex <- setsIndex + 1

      tempSet <- resultsNotNAs[[i]]

      sets <- rbind(sets, data.frame(
        TimeStart = Abundance[i, 1],
        TimeStop = NA,
        Present = I(list(tempSet))
      ))
    }
  }
  sets <- sets[1:setsIndex, ]

  sets
}

GeneralisedLotkaVolterra <- function(
  t, y, parms
) {
  # t is the current time point in the integration,
  # y is the current estimate of the variables in the ODE system.
  # parms is a vector or list of parameters
  # The return value of func should be a list, whose first element is a vector
  # containing the derivatives of y with respect to time, and whose next
  # elements are global values that are required at each point in times.
  # The derivatives must be specified in the same order as... y.
  # -- DESolve::rk documentation
  #
  # Our parms contains interaction matrix a and reproduction rates r.

  with(as.list(parms), {
    list(y * (r + a %*% y))
  })
}

addMissingEquilibriaEntries <- function(i, set, nnz, master) {
  # meant to be used in lapply
  # i is the index of the set we are using
  # set is a list of subcommunities
  # nnz is the nonzero values that we are going to add zeros to
  # master is the full community that we need to add zeros for.
  stopifnot(is.list(set))
  stopifnot(is.list(nnz))
  retval <- rep(0, length(master))
  setInd <- 1
  masInd <- 1
  while (setInd <= length(set[[i]])) {
    if (set[[i]][setInd] == master[masInd]) {
      retval[masInd] <- nnz[[i]][setInd]
      masInd <- masInd + 1
      setInd <- setInd + 1
    } else if (set[[i]][setInd] < master[masInd]) {
      setInd <- setInd + 1
    } else {
      masInd <- masInd + 1
    }
  }
  return(retval)
}
