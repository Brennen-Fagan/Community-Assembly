# See MNA-Image-ExampleOutcome-Presence.R and
# MNA-Image-ExampleOutcome-TimeJaccard.R for examples in context.

Calculate_TimeJaccard <- function(loaded, nspecies, minTime = NULL) {
  loaded$Abundance[, -1] <-
    loaded$Abundance[, -1] > loaded$Parameters$EliminationThreshold

  if (#"ReactionTime" %in% names(loaded) &&
      is.null(minTime)) {
    minTime = max(loaded$Abundance[, 1])/101 # Slice into 100 comparisons
  }

  ### Calculate Diversity: ##################################################
  diversity <- lapply(
    1:loaded$NumEnvironments,
    function(i, abund, numSpecies) {
      print(i)
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      env_basal <- env[, 1:nspecies[1]]
      env_consumer <- env[, nspecies[1] + 1:nspecies[2]]

      consistentDistance <- max(diff(time), minTime)
      targets <- seq(from = min(time), to = max(time), by = consistentDistance)
      rows <- unique( # Just In Case?
        sapply(targets, function(x, y) {which.max(y >= x)}, y = time)
      )

      helper <- function(i, r, m)
        as.numeric(suppressWarnings(vegan::vegdist(
          x = rbind(m[r[i - 1],], m[r[i],]), method = "jaccard"
        )))

      # Each Row is an "Environment"
      # We need to have consistent "distances" between each "environment".
      # We only need one "Environment" and its Next "Environment".
      Jaccard <- sapply(2:length(rows), helper, r = rows, m = env)
      Jaccard_basal <- sapply(2:length(rows), helper, r = rows, m = env_basal)
      Jaccard_cons <- sapply(2:length(rows), helper, r = rows, m = env_consumer)

      data.frame(Time = time[rows[-length(rows)]],
                 JaccardTemporal = Jaccard,
                 JaccardTemporal_Basal = Jaccard_basal,
                 JaccardTemporal_Consumer = Jaccard_cons,
                 Environment = as.character(i),
                 stringsAsFactors = FALSE)
    },
    abund = loaded$Abundance,
    numSpecies = sum(nspecies)
  )

  diversity <- dplyr::bind_rows(diversity) %>% tidyr::pivot_longer(
    cols = tidyr::starts_with("Jaccard"),
    names_to = "Measurement", values_to = "Value"
  )
  diversity_avg <- diversity %>% dplyr::group_by(
    Time, Measurement
  ) %>% dplyr::summarise(
    Value = mean(Value),
    .groups = "drop"
  ) %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    Environment = "Mean"
  )

  ### Return Diversities: ###############################################
  return(dplyr::bind_rows(
    diversity,
    diversity_avg
  ))
}
