Productivity <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }
  
  if (is.character(Populations)) {
    pop <- CsvRowSplit(Populations)
  } else {
    pop <- Populations
  }
  
  if (length(com) == 0) {
    return(NA)
  }
  
  comMatPos <- InteractionMatrix[com, com]; comMatPos[comMatPos < 0] <- 0
  poolRepPos <- Pool$ReproductionRate[com]; poolRepPos[poolRepPos < 0] <- 0
  
  parameters <- list(
    a = comMatPos,
    r = poolRepPos
  )
  
  return(
    sum(Dynamics(0, pop, parameters)[[1]] * pop / sum(pop))
  )
}