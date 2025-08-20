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

IslandDynamics <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  DispersalPool, # Species related dispersal rates
                 # Should have length == nrow(Pool). Multiplied by entries of
  DispersalIsland, # Island related dispersal rates. Is a matrix, row = to.
  Dynamics = IslandLotkaVolterra,
  Tolerance = 1E-1,
  Times = seq(from = 0,
              to = 2E4,
              by = 5E2),
  #DispersalRates, # For now, a matrix, to-from notation, for each island.
  #(NOT species on each island, since that requires the user knowing things in advance...)
  #DispersalRates, # List of matrices: column = species, row = (TO other) island, entry = travel rate
  Method = NULL, # Passed through to deSolve.
  TimesEvents = NULL,
  Verbose = FALSE
) {
  preprocessed <- IslandPreprocess(
    Pool = Pool,
    InteractionMatrix = InteractionMatrix,
    Communities = Communities,
    Populations = Populations,
    DispersalPool = DispersalPool,
    DispersalIsland = DispersalIsland,
    Tolerance = Tolerance,
    Verbose = Verbose
  )

  abundance <- with(
    preprocessed,
    deSolve::ode(
      abundance_init,
      times = Times,
      func = Dynamics,
      parms = parameters,
      events = list(func = function(t, y, parms) {
        with(parms, {
          if (exists("Verbose")) {
            if (Verbose) print(paste(t, "Zeroing:",
                                     toString(which(y < parms$epsilon))))
          }
        })
        y[y < parms$epsilon] <- 0
        y
      }, time = if(is.null(TimesEvents)) Times else TimesEvents),
      method = Method
    )
  )

  return(abundance)
}


