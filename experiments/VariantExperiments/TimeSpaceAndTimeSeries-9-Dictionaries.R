library(dplyr)

# Pool-Patch: #################################################################
poolpatchDictionaryOrigin <- expand.grid(
  BasalConsumerRatio = 1/2,
  NSpecies = c(100, 200),
  PoolFunction = "RMTRCode2::LawMorton1996_species",
  # PoolParameters = c(
  #   paste("Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
  #         "LogBodySize = c(-2, -1, -1, 0)", sep = "; ")
  # ),
  PoolK1InteractionEffectiveness = 10^c(-2, -1.5, -1, -0.5),
  PoolK2ConsumerSizeAdvantage = 10^c(1, 1.5),
  PoolK3ConsumerPredationRange = 10^c(-1, -0.5, log10(0.5), 0),
  PoolK4ConsumerEfficiency = c(0.05, 0.1, 0.15, 0.2),
  PoolK5BasalBiomass = c(30, 100, 300),
  PoolK6CoefOfVariation = c(0, 0.1, 0.2),
  PoolBasalLogBodySize = c("c(-2, -1)", "c(-2, 0)", "c(-3, -1)"),
  PoolConsumerLogBodySize = c("c(-1, -0)", "c(-1, 1)", "c(-1, 2)"),
  PoolDispersalSpeed = 1, # Value divided by DispersalResistance to get current.
  NumberEnvironments = c(1, 2, 6, 10, 12),
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# Dynamics: ###################################################################
dynamicsDictionaryOrigin <- expand.grid(
  InteractionFunction = "RMTRCode2::LawMorton1996_CommunityMat",
  # InteractionParameters = "Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
  InteractionK1InteractionEffectiveness = 10^c(-2, -1.5, -1, -0.5),
  InteractionK2ConsumerSizeAdvantage = 10^c(1, 1.5),
  InteractionK3ConsumerPredationRange = 10^c(-1, -0.5, log10(0.5), 0),
  InteractionK4ConsumerEfficiency = c(0.05, 0.1, 0.15, 0.2),
  InteractionK5BasalBiomass = c(30, 100, 300),
  InteractionK6CoefOfVariation = c(0, 0.1, 0.2),
  InteractionEliminationThreshold = 10^c(-5, -4, -3),
  DynamicsFunction = "RMTRCode2::PerCapitaDynamics_Type1",
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# Events: #####################################################################
eventsDictionaryOrigin <- expand.grid(
  EventsFunction = "defaultEvents", # Takes Number of Environments and Species.
  EventsNumberMultiplier = c(1, 2, 10, 20), # I.e. how many more events
  ImmigrationMultiplier = c(0.1, 1, 10), # Mult. is on rate and number.
  ImmigrationFunction = "RMTRCode2::ArrivalFUN_Example2",
  ColonizationPropaguleSize = c(0.04, 0.4, 4), # Initial abundance on Colonize.
  ExtirpationMultiplier = c(0.1, 1, 10), # Frequency multiplier
  ExtirpationFunction = "RMTRCode2::ExtinctFUN_Example2",
  ExtirpationProportion = c(1, 0.9, 0), # Proportion of population removed.
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# Initial Conditions: #########################################################
initialConditionsDictionaryOrigin <- rbind(
  data.frame(Species = "None", Method = NA, Argument = NA,
             stringsAsFactors = FALSE),
  data.frame(Species = c("Basal", "All"),
             Method = c("Solve",  "InverseSize"), Argument = NA,
             stringsAsFactors = FALSE),
  expand.grid(
    Species = c("Basal", "All"),
    Method = c("Set", "Random"),
    Argument = c(1, 100, 1000),
    stringsAsFactors = FALSE
  )
) %>% dplyr::mutate(ID = dplyr::row_number())

# Dispersal: ##################################################################
# Note: ID is not the row number, because access is via:
# dispersalDictionaryOrigin[ifelse(is.na(dispersalDictionaryChoice),
#                           1, dispersalDictionaryChoice + 2), ]
# i.e. NA -> 1 (None), 0 -> 2 (1e0), 1 -> 3 (1e1), etc.
dispersalDictionaryOrigin <- rbind(
  data.frame(Resistance = Inf, Configuration = "None"),
  expand.grid(
    Resistance = 10^c(0:9),
    Configuration = c("Ring", "Line", "Complete"),
    stringsAsFactors = FALSE
  ))
dispersalDictionaryOrigin <- dispersalDictionaryOrigin %>% dplyr::mutate(
  ID = c(NA, 0:nrow(dispersalDictionaryOrigin))[dplyr::row_number()]
)

# Distance/Intensity/Affinity Effect: #########################################
distanceDictionaryOrigin <- data.frame(
  rhofunction = c( # Take patch and species affinities as vector arguments.
    "rho.noop", # Just returns 1.
    "rho.2.0.1.euclidean",
    "rho.2.1.2.euclidean",
    "rho.5.0.1.euclidean",
    "rho.5.1.2.euclidean",
    "rho.10.0.1.euclidean",
    "rho.10.1.2.euclidean"
  ),
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# Pool-Patch Affinity/Adaptation: #############################################
affinityDictionaryOrigin <- expand.grid(
  SpeciesAffinities = c(
    "rep_0",                 # Pool with {0} affinities.
    "rep_0.5",               # Pool with {0.5} affinities.
    "rep_1",                 # Pool with {1} affinities.
    "sample.int.normalized", # Pool {0, 1} patch affinities at random.
    "sample.int.3",          # Pool {0, 0.5, 1} patch affinities at random.
    "runif",                 # Pool [0, 1] patch affinities at random.
    "evensplit_01"           # Pool {0, 1} alternating affinities.
  ),
  PatchAffinities = c(
    # Detection via if string begins with a numeric or a non-numeric.
    # If numeric, it takes it as a fixed set of affinities.
    # If non-numeric, it attempts to treat the string as a function name.
    # In the latter case, it provides ONLY NumberEnvironments as an argument.
    "rep_0", #               Patch {0} affinities.
    "rep_0.25", #            Patch {0.25} affinities.
    "rep_0.5", #             Patch {0.5} affinities.
    "rep_0.75", #            Patch {0.75} affinities.
    "rep_1", #               Patch {1} affinities.
    "gradientline_01", #     Patch {0, 1} affinities. Gradient Line.
    "evensplit_01", #        Patch {0, 1} affinities. Alternating.
    "gradientline_0half1", # Patch {0, 0.5, 1} affinities. Gradient Line.
    "patchTypes.0.Half.1", # Patch {0, 0.5, 1} affinities. Gradient Ring.
    "runifRing", #           Patch [0, 1] affinities. Gradient Ring at Random.
    "evensplit_0.51" #        Patch {0.5, 1} affinities. Alternating.
  ),
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# Intervention Patch Affinities: ##############################################
interventionPatchDictionaryOrigin <- expand.grid(
  PatchAffinities = c(
    # Detection via if string begins with a numeric or a non-numeric.
    # If numeric, it takes it as a fixed set of affinities.
    # If non-numeric, it attempts to treat the string as a function name.
    # In the latter case, it provides ONLY NumberEnvironments as an argument.

    "rep_0", #               Patch {0} affinities.
    "rep_0.25", #            Patch {0.25} affinities.
    "rep_0.5", #             Patch {0.5} affinities.
    "rep_0.75", #            Patch {0.75} affinities.
    "rep_1", #               Patch {1} affinities.
    "gradientline_01", #     Patch {0, 1} affinities. Gradient Line.
    "sample.int.normalized", # Patches -> {0, 1} Unif @ Random
    "patchTypes.0.Half.1", # Patches -> {0, 0.5, 1} Gradient Ring
    "sample.int.3", # Patches -> {0, 0.5, 1} Unif @ Random
    "runifRing", # Patches -> [0, 1] Gradient Ring
    "runif" # Patches -> [0, 1] Unif @ Random
  ),
  InterventionLocation = c(
    # Percentage seems easiest
    NA, # == Random
    1, # Last
    0 # First
  ),
  InterventionPercentage = c(
    1/3,
    0.5,
    2/3,
    1
  ),
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# InterventionTime: ###########################################################
interventionTimeDictionaryOrigin <- data.frame(
  # Time1, Time2; called by eval(str2lang(X)) where X is the string below
  #               and "loaded" is the file that is loaded.
  Time1 = c(
    "median(loaded$Events$Times)",
    "quantile(loaded$Events$Times, p = 0.25)"
  ),
  Time2 = c(
    "1/2 * max(loaded$Events$Times)",
    "quantile(loaded$Events$Times, p = 0.75)"
  ),
  Method = c(# each needs a custom implementation unfortunately!
    "mean",
    "runif"
  ),
  InterventionTimespan = c(
    0 # Instantaneous => Switch
    # Else: Should be numeric > 0, determines timespan for interpolation.
  )
) %>% dplyr::mutate(ID = dplyr::row_number())
