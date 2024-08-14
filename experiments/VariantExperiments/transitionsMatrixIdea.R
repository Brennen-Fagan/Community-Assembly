library(dplyr)
source("abundanceToOccupancy64.R")
library(ctmcd)
# Problem: ctmcd likely assumes evenly spaced observation times.
# So the below doesn't work; we'd need to interpolate (theoretically trivial
# for occupancy, probably linearly for abundance, but not computationally 
# trivial in either case).

transitions <- dplyr::bind_rows(
  lapply(seq_along(results), function(i) {
    data.frame(
      Time = results[[i]]$Abundance[, 1],
      State = abundanceToOccupancy64(
        results[[i]]$Abundance, 
        epsilon = results[[i]]$Parameters$EliminationThreshold
      )
    ) %>% dplyr::mutate(
      # State = factor(State), 
      NextState = lead(State),
      # NumState = as.numeric(State),
      # NumNext = as.numeric(NextState)
      Simulation = i
    ) 
  }
  )
) %>% dplyr::mutate(
  FacState = factor(State), # Keep track of Character -> Factor -> Numeric.
  NumState = as.numeric(State), # Make more human readable
  NumNext = as.numeric(factor(NextState, levels = levels(FacState)))
)

transitionsMatrix <- with(transitions, table(NumState, NumNext))
