standardizeGapSize <- function(abundance, gap, timespan = NULL, epsilon = NULL) {
  # format: abundance <- cbind(Time column, species densities)
  # format: gap <- desired step size
  # format: epsilon = NULL
  if(!is.null(epsilon)) {
    abundance[, -1][abundance[, -1] < epsilon] <- 0
  }

  if (is.null(timespan)) {
    timespan <- range(abundance[, 1])
  }
  timesout <- seq(timespan[1], timespan[2], by = gap)
  splines <- lapply(2:ncol(abundance), function(col) {
    spline(abundance[, 1], abundance[, col], xout = timesout)$y
  })

  retval <- cbind(timesout, do.call(cbind, splines))
  if(!is.null(epsilon)) {
    retval[, -1][retval[, -1] < epsilon] <- 0
  }
  return(retval)
}
