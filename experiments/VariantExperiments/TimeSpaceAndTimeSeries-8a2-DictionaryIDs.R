# Okay, so our last set of runs with the new apparatus all seem to give smaller
# than expected richness values. But if I rerun an old pool, patch, and event
# system, I get the same values (more-or-less).
# Furthermore, attaching affinities that we expect to do nothing, does do
# nothing, and the statistical properties of the pool's sizes, reproduction 
# rates, consumer-consumed interaction strengths, event types, event-species
# frequencies, and waiting times between events all seem not to have varied
# when I'm careful as I was in the last set of runs.
# So what's going on???

# What's left: Number of species, Number of patches, characteristic time scales.
# None of these *should* affect the richness values -- not like we're seeing.
# But they're the only things that look viable to check.

library(dplyr)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Dictionaries.R"))

# Non-spatial Experiments, as agreed email and meeting 2024-09-26.
# I've added in the low intensity, which we otherwise agreed to remove,
# because we are specifically looking sometimes at what happens when the pool
# and patches have a complete mismatch in adaptation values.
experiments <- list(
  ppDO =
    poolpatchDictionaryOrigin %>% dplyr::filter(
      BasalConsumerRatio == 0.5, # 2 Consumers for Each Basal
      NSpecies == 100,
      # Standard, and only implemented, LM1996 Function and Parameters,
      PoolDispersalSpeed == 1,
      NumberEnvironments == 10
    ),
  
  dynDO =
    dynamicsDictionaryOrigin[1, ], # No other options implemented.
  
  eDO =
    eventsDictionaryOrigin %>% dplyr::filter(
      EventsNumberMultiplier == 2,
      ImmigrationMultiplier == 1, 
      ExtirpationMultiplier == 1,
      ExtirpationProportion == 1 # Extirpation == Extinction in a 1 patch system.
    ),
  
  icDO = initialConditionsDictionaryOrigin %>% dplyr::filter(
    (Species == "None") 
  ),
  
  dispDO =
    dispersalDictionaryOrigin %>% dplyr::filter(
      Configuration == "None" # Doesn't make sense for 1 patch systems.
    ),
  
  aDO =
    affinityDictionaryOrigin %>% dplyr::filter(
      SpeciesAffinities %in% c("rep_0"),
      PatchAffinities %in% c("rep_0.5")
    ),
  
  distDO =
    distanceDictionaryOrigin %>% dplyr::filter(
      rhofunction %in% c("rho.10.1.2.euclidean")
    ),
  
  iPDO =
    # Note: Every combination of iPDO with ppDO.
    interventionPatchDictionaryOrigin %>% dplyr::filter(
      PatchAffinities %in% c("rep_0", "rep_0.25", "rep_0.5", "rep_0.75", "rep_1"),
      InterventionLocation == 1,
      InterventionPercentage == 1
    ),
  
  iTDO =
    interventionTimeDictionaryOrigin %>% dplyr::filter(
      Method == "mean"
    ),
  
  # interventionDispersalDictionaryChoice
  iDispChoice = "p",
  # interventionDistanceDictionaryChoice
  iDistChoice = "p"
)