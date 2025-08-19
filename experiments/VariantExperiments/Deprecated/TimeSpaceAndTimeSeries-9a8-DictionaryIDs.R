library(dplyr)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9-Dictionaries.R"))

# Initial results from the 9a round suggest increasing the consumer body range
# upper bound and decrease the consumer consumption range k3.
modifiedcase <- basecase
modifiedcase$PoolConsumerLogBodySize <- "c(-1, 0.5)"
modifiedcase$PoolK3ConsumerPredationRange <- 10^-1
modifiedcase$InteractionK3ConsumerPredationRange <- 10^-1
modifiedcase$InteractionK5BasalBiomass <- 30
modifiedcase$InteractionEliminationThreshold <- 10^-5
modifiedcase$ColonizationPropaguleSize <-
  modifiedcase$InteractionEliminationThreshold * 4 * 10^3

experiments <- list(list(
  ppDO =
    poolpatchDictionaryOrigin %>% dplyr::filter(
      BasalConsumerRatio == 0.5, # 2 Consumers for Each Basal
      NSpecies == 200,
      # Standard, and only implemented, LM1996 Function and Parameters,
      PoolDispersalSpeed == 1,
      NumberEnvironments == 1,

      apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
        #             (despite documentation saying otherwise)
        .cols = dplyr::any_of(names(modifiedcase)),
        .fns = function(colval) {
          colname <- dplyr::cur_column()
          modifiedcase[[colname]] == colval
        }
      ), 1, all)
    ),

  dynDO =
    dynamicsDictionaryOrigin %>% dplyr::filter(
      apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
        #             (despite documentation saying otherwise)
        .cols = dplyr::any_of(names(modifiedcase)),
        .fns = function(colval) {
          colname <- dplyr::cur_column()
          modifiedcase[[colname]] == colval
        }
      ), 1, all)
    ),

  eDO =
    eventsDictionaryOrigin %>% dplyr::filter(
      EventsNumberMultiplier == 20, # Longer simulation. Keep number same while
      ImmigrationMultiplier == 0.1, # Decreasing rate of occurrence.
      ExtirpationMultiplier == 0.1,
      ExtirpationProportion == 1, # Extirpation == Extinction in a 1 patch system.
      apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
        #             (despite documentation saying otherwise)
        .cols = dplyr::any_of(names(modifiedcase)),
        .fns = function(colval) {
          colname <- dplyr::cur_column()
          modifiedcase[[colname]] == colval
        }
      ), 1, all)
    ),

  icDO = initialConditionsDictionaryOrigin %>% dplyr::filter(
    Species == "None"
  ),

  dispDO =
    dispersalDictionaryOrigin %>% dplyr::filter(
      Configuration == "None" # Doesn't make sense for 1 patch systems.
    ),

  distDO =
    distanceDictionaryOrigin %>% dplyr::filter(
      rhofunction %in% c("rho.2.1.2.euclidean",
                         "rho.5.1.2.euclidean",
                         "rho.10.1.2.euclidean")
    ),

  aDO =
    affinityDictionaryOrigin %>% dplyr::filter(
      SpeciesAffinities %in% c("rep_0", "rep_1", "runif", "evensplit_01"),
      PatchAffinities %in% c("rep_0", "rep_0.25", "rep_0.5",
                             "rep_0.75", "rep_1")
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
      ID == 1 # averaged central tendency.
    ),

  # interventionDispersalDictionaryChoice
  iDispChoice = "p",
  # interventionDistanceDictionaryChoice
  iDistChoice = "p"
))
