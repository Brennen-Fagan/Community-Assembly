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
  IslandDispersalMatrix, # Amount of population on island i goes to j as [i, j].
  IslandPopulationList, # Of same length as nrow(IslandDistanceMatrix).
  Dynamics = GeneralisedLotkaVolterra
) {
  # Calculate how much of what is on each island.
  #
}
