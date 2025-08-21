thinAbundanceEqualTimeSteps <- function(
  abundance, threshold, preferredTimeStep, preferredStart = NULL,
  includeMinTime = FALSE, minTimePlus = NULL
) {
  time <- abundance[, 1]
  consistentDistance <- max(diff(time), preferredTimeStep)
  if (consistentDistance != preferredTimeStep) {
    warning(paste0("preferredTimeStep is too small, minimum:",
                   consistentDistance))
  }
  if (is.null(preferredStart)) preferredStart <- min(time)
  targets <- seq(from = preferredStart,
                 to = max(time), by = consistentDistance)
  if (includeMinTime && preferredStart != min(time))
    targets <- c(min(time), targets)
  if (!is.null(minTimePlus)) {
    targets <- sort(c(targets, min(time)+minTimePlus))
  }

  rows <- unique( # Just In Case?
    sapply(targets, function(x, y) {which.max(y >= x)}, y = time)
  )

  abundance <- abundance[rows, ]
  abundance[, 1] <- targets
  # Technically an approximation, but we should be high resolution
  # enough for it not to be a problem.

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- abundance[, -1] < threshold # & abundance[, -1] > 0
  abundance[, -1][toEliminate] <- 0

  return(abundance)
}
