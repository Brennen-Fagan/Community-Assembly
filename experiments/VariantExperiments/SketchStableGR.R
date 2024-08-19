datfolder <-
  "TSTS_Simulations_18-1_9-9_2024-08-06"
results <-
  lapply(dir(datfolder, "Simulation", full.names = T),
         function(x) {nm <- load(x); return(get(nm))})
interventions <-
  lapply(dir(datfolder, "Intervention", full.names = T),
         function(x) {nm <- load(x); return(get(nm))})
preferredGap <-
  quantile(unlist(lapply(results,
                         function(r) diff(r$Abundance[, 1]))), p = 0.05)
preferredTimes <-
  c(500, min(unlist(lapply(interventions, function(r) r$Abundance[1, 1]))))

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

# Full system:
# Error in solve.default(Tee, W) :
#   system is computationally singular: reciprocal condition number = 6.49124e-19
# Burnt in and Halved:
# Error in solve.default(Tee, W) :
#   system is computationally singular: reciprocal condition number = 8.56199e-18
# Not sure how to proceed.

estimate_GR <-
  stableGR::stable.GR(lapply(temp, function(x) x[1:minnrow, nonzeros2[-1]]))
# Works, but takes a while!
# Returns are a bit strange: mpsrf < 1.

# Need to reimplement beta diversity but for contrasting between two results.
