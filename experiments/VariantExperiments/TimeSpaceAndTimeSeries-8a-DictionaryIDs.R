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
      NSpecies == 200,
      # Standard, and only implemented, LM1996 Function and Parameters,
      PoolDispersalSpeed == 1,
      NumberEnvironments == 1
    ),

  dynDO =
    dynamicsDictionaryOrigin[1, ], # No other options implemented.

  eDO =
    eventsDictionaryOrigin %>% dplyr::filter(
      EventsNumberMultiplier == 2, # Longer simulation. No changes to balance.
      ExtirpationProportion == 1 # Extirpation == Extinction in a 1 patch system.
    ),

  dispDO =
    dispersalDictionaryOrigin %>% dplyr::filter(
      Configuration == "None" # Doesn't make sense for 1 patch systems.
    ),

  aDO =
    affinityDictionaryOrigin %>% dplyr::filter(
      SpeciesAffinities %in% c("evensplit_01", "rep_0", "rep_1", "runif"),
      PatchAffinities %in% c("rep_0", "rep_0.25", "rep_0.5", "rep_0.75", "rep_1")
    ),

  distDO =
    distanceDictionaryOrigin %>% dplyr::filter(
      rhofunction %in% c("rho.2.1.2.euclidean",
                         "rho.5.1.2.euclidean",
                         "rho.10.1.2.euclidean")
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
