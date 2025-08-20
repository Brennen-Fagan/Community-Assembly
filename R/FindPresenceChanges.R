FindPresenceChanges <- function(loaded) {
  abund <- loaded$Abundance
  abund[, -1] <- abund[, -1] > loaded$Parameters$EliminationThreshold
  numSpecies <- (ncol(abund) - 1) / loaded$NumEnvironments
  
  reports <- lapply(1:numSpecies, function(i, a) {
    # Looking for a data.frame(
    #  Time, Species, Environment, Type = c(Immigration, Extirpation),
    #  Regional = c(True, False), Neutral = c(True, False)
    # )
    
    time <- a[, 1]
    x <- a[, i + numSpecies * (1:loaded$NumEnvironments - 1) + 1]
    xdiff <- apply(x, 2, function(y) y - dplyr::lag(y))
    
    # First time after a change occurs. Note this double lists.
    changelocs <- apply(xdiff, 2, function(y) list(which(y != 0)))
    if(length(unlist(changelocs)) == 0) {
      return(NULL)
    }
    
    reportspatch <- lapply(1:loaded$NumEnvironments, function(j, y, locs) {
      if (length(locs[[j]][[1]]) == 0) return(NULL)
      stopifnot(((as.numeric(colnames(y)[j]) - 1) %/% numSpecies) + 1 == j)
      
      data.frame(
        Time = time[locs[[j]][[1]] - 1],
        Species = ((i - 1) %% numSpecies) + 1,
        Environment = as.character(j),
        Type = dplyr::case_when(
          y[locs[[j]][[1]], j] ==  1 ~ "Immigration",
          y[locs[[j]][[1]], j] == -1 ~ "Extirpation",
          TRUE ~ "OOPS"
        )
      )
    }, y = xdiff, locs = changelocs) %>% dplyr::bind_rows()
    
    changelocs <- sort(unique(unlist(changelocs)))
    
    # Regional Report:
    changesregional <- lapply(changelocs, function(loc, y) {
      # For each change loc, look at the timestep before and of.
      target <- y[-1:0 + loc, ]
      # If all 0 before and any 1 of -> "Immigration",
      # If any 1 before and all 0 of -> "Extirpation"
      # Otherwise, discard.
      dplyr::case_when(
        all(target[1, ] == 0) && any(target[2, ] == 1) ~ "Immigration",
        any(target[1, ] == 1) && all(target[2, ] == 0) ~ "Extirpation",
        TRUE ~ "Discard"
      )
    }, y = x)
    
    # print(changelocs)
    # print(time[changelocs])
    # print(changesregional)
    
    reportsregion <- data.frame(
      Time = time[changelocs - 1],
      Species = ((i - 1) %% numSpecies) + 1,
      Type = unlist(changesregional),
      Regional = TRUE
    ) %>% dplyr::filter(Type != "Discard")
    
    reports <- dplyr::full_join(
      reportspatch, reportsregion,
      by = c("Time", "Species", "Type")
    ) %>% dplyr::mutate(
      Regional = ifelse(is.na(Regional), FALSE, Regional)
    )
    
    return(reports)
    
  }, a = abund)
  
  dplyr::bind_rows(reports)
}