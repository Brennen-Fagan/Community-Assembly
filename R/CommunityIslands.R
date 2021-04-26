# Islands and Interactions
# So for our initial interaction, we need to have a graph of patches.
# Each patch contains an ecosystem.
# Patches are effectively homogeneous at this point.
# Heterogeneity is hard to decide how to introduce, given that the obvious
# answer would be to modify the effective species carrying capacities, but
# those are implicit in the model, appearing due to reproduction rates and
# intra-species competition.
# Patches also need their dynamics under which the community is run.
# Then there are dispersal rates which are some function of species and distance
# from the patch perspective.
# From this perspective, islands are closely related sets of patches in that
# their adjacent distances are 0 (which does not mean instantaneous travel!).
#       Patch                        Patch
#   |-----------|                |-----------|
#   | Community |    Dispersal   | Community |
#   | Dynamics  | -------------> | Dynamics  |
#   |           | <------------- |           |
#   |           |                |           |
#   |-----------|                |-----------|
#
#
# Inputs: Species Pool, Species Interaction Matrix, which imply dynamics
# Species present on each patch, Patch distances

CsvRowSplit <- function(csv) {
  return(as.numeric(unlist(strsplit(csv, split = ", "))))
}

FindSteadyStateFromEstimate <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
  Tolerance = 1E-1
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }

  if (is.character(Populations)) {
    init <- CsvRowSplit(Populations)
  } else {
    init <- Populations
  }
  init_min <- min(init)
  init_max <- max(init)

  anyZeroOrNotSame <- TRUE # Set T to make at least one look.
  epsilon <- Tolerance

  # Run and check to see if anyone dies (bad) or changes (not at steady).
  while (anyZeroOrNotSame) {
    init_old <- init

    init <- rootSolve::steady(
      y = init,
      func = Dynamics,
      parms = list(a = InteractionMatrix[com, com],
                   r = Pool$ReproductionRate[com],
                   epsilon = epsilon),
      positive = TRUE
    )$y
    print(init)

    if (any(init < epsilon)) {
      # Someone died, reset to random location.
      print("Died")
      anyZeroOrNotSame <- TRUE
      init[init < epsilon] <- runif(
        n = sum(init < epsilon),
        min = init_min,
        max = init_max
      )
    } else if (any(round(init / init_old, 1) != 1)) {
      # Not done, keep going.
      print("Changed")
      anyZeroOrNotSame <- TRUE
    } else {
      anyZeroOrNotSame <- FALSE
    }
  }

  return(init)
}

Productivity <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
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

IslandDynamics <- function(
  Pool,
  InteractionMatrix,

) {
  # Calculate how much of what is on each island.
  #
}
