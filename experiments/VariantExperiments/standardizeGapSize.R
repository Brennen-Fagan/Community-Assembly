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

datfolder <-
  "TSTS_Simulations_18-1_9-9_2024-08-06"
results <-
  lapply(dir(datfolder, "Simulation", full.names = T),
         function(x) {nm <- load(x); return(get(nm))})
preferredGap <-
  quantile(unlist(lapply(results,
                         function(r) diff(r$Abundance[, 1]))), p = 0.03)
preferredTimes <-
  c(500, min(unlist(lapply(results,
                    function(r) r$Abundance[nrow(r$Abundance), 1]/2))))

temp <- lapply(
  seq_along(results),
  function(i) standardizeGapSize(
    results[[i]]$Abundance,
    gap = preferredGap,
    timespan = preferredTimes,
    epsilon = results[[i]]$Parameters$EliminationThreshold
  )
)

min(unlist(lapply(temp, nrow))) -> minnrow
nonzeros2 <- lapply(temp, function(x) which(!(apply(x, MARGIN = 2, var) == 0)))
nonzeros2 <- Reduce(intersect, nonzeros2)
nonzeros2
estimate_neff <-
  stableGR::n.eff(lapply(temp, function(x) x[1:minnrow, nonzeros2[-1]]))

# Error in solve.default(Tee, W) :
#   system is computationally singular: reciprocal condition number = 6.49124e-19
#???

estimate_GR <-
  stableGR::stable.GR(lapply(temp, function(x) x[1:minnrow, nonzeros2[-1]]))
# Works, but takes a while!
# Returns are a bit strange: mpsrf < 1.
