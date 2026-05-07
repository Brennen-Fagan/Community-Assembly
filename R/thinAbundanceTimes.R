thinAbundanceTimes <- function(abundance, threshold, times) {
  time <- abundance[, 1]

  if (any(times < min(time))) {
    warning("times before start of time requested. Removing.")
    times[times < min(time)] <- NULL
  }
  if (any(times > max(time))) {
    warning("times after end of time requested. Removing.")
    times[times > max(time)] <- NULL
  }

  rows <- sapply(times, function(x, y) {which.max(y >= x)}, y = time)

  abundance <- abundance[rows, ]
  abundance[, 1] <- times
  # Technically an approximation, but we should be high resolution
  # enough for it not to be a problem.

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- abundance[, -1] < threshold # & abundance[, -1] > 0
  abundance[, -1][toEliminate] <- 0

  return(abundance)
}
