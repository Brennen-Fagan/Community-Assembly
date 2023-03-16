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
