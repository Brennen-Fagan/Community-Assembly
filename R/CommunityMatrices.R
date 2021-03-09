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
