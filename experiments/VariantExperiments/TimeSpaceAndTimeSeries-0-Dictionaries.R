library(dplyr)

# Pool-Patch: #################################################################
poolpatchDictionaryOrigin <- expand.grid(
  BasalConsumerRatio = 1/2,
  NSpecies = c(100, 200),
  PoolFunction = "RMTRCode2::LawMorton1996_species",
  PoolParameters = c(
    paste("Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
          "LogBodySize = c(-2, -1, -1, 0)", sep = "; ")
  ),
  PoolDispersalSpeed = 1, # Value divided by DispersalResistance to get current.
  NumberEnvironments = c(1, 2, 6, 10, 12),
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# Dynamics: ###################################################################
dynamicsDictionaryOrigin <- data.frame(
  InteractionFunction = "RMTRCode2::LawMorton1996_CommunityMat",
  InteractionParameters = "Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
  DynamicsFunction = "RMTRCode2::PerCapitaDynamics_Type1",
  stringsAsFactors = FALSE
) %>% dplyr::mutate(ID = dplyr::row_number())

# Events: #####################################################################
eventsDictionaryOrigin <- expand.grid(
  EventsFunction = "defaultEvents", # Takes Number of Environments and Species.
  EventsNumberMultiplier = c(1, 2, 10, 20), # I.e. how many more events
  ImmigrationMultiplier = c(0.1, 1, 10), # Mult. is on rate and number.
  ImmigrationFunction = "RMTRCode2::ArrivalFUN_Example2",
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
    # toString(rep(0, NumberOfEnvironments)), # Patches -> {0}
    # toString(rep(0.25, NumberOfEnvironments)), # Patches -> {0.25}
    # toString(rep(0.5, NumberOfEnvironments)), # Patches -> {0.5}
    # toString(rep(0.75, NumberOfEnvironments)), # Patches -> {0.75}
    # toString(rep(1, NumberOfEnvironments)), # Patches -> {1}
    # toString(c(rep(0, NumberOfEnvironments/2),
    #            rep(1, NumberOfEnvironments/2))), # Patches -> {0, 1} Gradient
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
