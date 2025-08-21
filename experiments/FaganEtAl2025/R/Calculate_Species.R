Calculate_Species <- function(result, bintimes = FALSE) {
  SpeciesPerEnvironment <- lapply(
    1:result$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      # Need to retrieve Position and Value
      species <- apply(
        cbind(time, env), MARGIN = 1,
        FUN = function(x) {
          time <- x[1]
          dat <- x[-1]
          if (any(dat > 0)) {
            positions <- (which(dat > 0))
            values <- dat[positions]
            data.frame(
              Time = time,
              Species = positions,
              Abundance = values,
              row.names = NULL
            )
          } else {NULL}
          # Returns as list
        }
      )
      return(
        dplyr::bind_rows(species) |> dplyr::mutate(
          Environment = i
        )
      )
    },
    abund = result$Abundance,
    numSpecies = (ncol(result$Abundance) - 1) / result$NumEnvironments
  )
  
  if (bintimes) {
    # Should equalise time steps.
    SpeciesPerEnvironment <- lapply(
      SpeciesPerEnvironment, function(SPE) {
        SPE |> dplyr::mutate(
          TimeFloor = floor(Time*10)/10
        ) |> dplyr::group_by(
          TimeFloor, Species, Environment
        ) |> dplyr::summarise(
          Abundance = median(Abundance, na.rm = TRUE)
        )
      })
  }
  
  return(dplyr::bind_rows(SpeciesPerEnvironment))
}