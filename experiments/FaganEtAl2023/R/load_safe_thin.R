load_safe_thin <- function(fname, bythin, divtime, burn) {
  loaded <- tryCatch({load(fname)},
                     error = function(e) {
                       print(fname)
                       print(e)
                       return(NA)
                     })
  if (is.na(loaded)) {
    return(NA)
  } else {
    loaded <- get(loaded)
    loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                             to = nrow(loaded$Abundance),
                                             by = bythin), ]
    
    toEliminate <-
      loaded$Abundance[, -1] <
      loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
    loaded$Abundance[, -1][toEliminate] <- 0
    
    loaded$Abundance <- loaded$Abundance[
      loaded$Abundance[, 1] > burn,
      ]
    
    loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divtime
    return(loaded)
  }
}