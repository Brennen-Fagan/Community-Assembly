# Foreach operators
`%dopar%` <- foreach::`%dopar%`
`%:%` <- foreach::`%:%`

CommunityMatrix <- function(method, ...) {
  if (method %in% c(
    "SpeciesPool", "Species Pool", "Law", "Morton", "LawMorton1996"
  )) {
    return(LawMorton1996(...))
  } else if (method %in% c(

  )) {

  } else {
    stop("Method not recognised.")
  }
}

SpeciesPresent <- function(Abundance, Threshold = 1E-4) {
  resultsNotNAs <- apply(
    Abundance[, -1],
    MARGIN = 1,
    FUN = function(row, epsilon) {
      which(!is.na(row) & row > epsilon)
    }, epsilon = Threshold)

  tempSet <- resultsNotNAs[[1]]
  sets <- data.frame(
    TimeStart = results$Abundance[1, 1],
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
        TimeStart = results$Abundance[i, 1],
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
