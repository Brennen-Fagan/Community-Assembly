# Introduction: ###############################################################
# Sequel to TimeSpaceAndTimeSeries-2a-Intervention.R and 6a-Simulation.R.
# We're using the base established by 6a to fork the simulations midway through
# with treatments that change the patch "affinities" (name pending).
#
# As always, please see the previous files for some design choices,
# although we aim to improve design at each stage.

# This does not show up as a parameter/ID choice.
runDictionaryChoice <-
  # 1 # "TSTS_Simulations_1-1_1-1_2024-02-13"
  # 2 # "TSTS_Simulations_1-1_2-2_2024-02-14"
  # 3 # "TSTS_Simulations_2-1_2-2_2024-02-14"
  # 4 # "TSTS_Simulations_6-1_2-2_2024-02-15"
  # 5 # "TSTS_Simulations_10-1_2-2_2024-02-15"
  6 # "TSTS_Simulations_11-1_3-3_2024-02-23"
  # 7 # "TSTS_Simulations_11-1_4-4_2024-02-23"

# Parameters: #################################################################
# Note that many of our options here vary between a deterministic mode
# and a stochastic mode. As such, we assign seeds still, but may not use them.

interventionPatchDictionaryChoice <-
  1 # Patches -> {0}
  # 2 # Patches -> {0.5}
  # 3 # Patches -> {1}
  # 4 # Patches -> {0, 1} Gradient
  # 5 # Patches -> {0, 1} Unif @ Random
  # 6 # Patches -> {0, 0.5, 1} Gradient
  # 7 # Patches -> {0, 0.5, 1} Unif @ Random
  # 8 # Patches -> [0, 1] Gradient
  # 9 # Patches -> [0, 1] Unif @ Random
interventionPatchSeedChoice <-
  1 # Used on ...

interventionTimeDictionaryChoice <-
  1 # Deterministic; set to 1/4 * max(Events$Time) + 1/2 median(Events$Time)
  # 2 # Stochastic; anywhere in [1/4 max(Events$Time), 3/4 max(Events$Time)]
  #   # For a double length run, there's at least half a run either side.
interventionTimeSeedChoice <-
  1 # Used on ...

# Not modifying dynamics at the moment.
# dynamicsDictionaryChoice <-
#   1 # Law and Morton 1996, Size-Structured Lotka-Volterra, Default Parameters
# #
# dynamicsSeedChoice <-
#   # 1 # Used on 2024-02-13
#   # 2 # Used on 2024-02-14
#   # 3 # Used on 2024-02-23
#   4 # Used on 2024-02-23 for the 2 patch system

# Probably shouldn't change.
interventionDispersalDictionaryChoice <-
  15 # c(NA, 5, 0)
# Index: Ones place is resistance to Dispersal on a log scale.
#      : Tens place is configuration: 0* = Ring, 1* = Line, 2* = Complete.
#      : Special: NA corresponds to no dispersal.
# Note : No randomness, so we don't need a seed.

# Probably shouldn't change.
# choose r' = r * rho ^ (sign(r)), but what rho?
interventionDistanceDictionaryChoice <- # for m, n in [0, 1], rho(m, n) = ...
  # 1 # 2 ^ (- euclid(m, n)) => rho in [1/2, 1] for 1-D
  2 # 2 ^ (1 - 2 euclid(m, n)) => rho in [1/2, 2] for 1-D

## Other Parameters: ##########################################################
# Most should be pulled from the data already.

directory <- "." # Should be "VariantExperiments"

source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Interventions.R"))

# Libraries: ##################################################################
library(RMTRCode2)
library(dplyr)
library(Matrix)


# Dictionaries: ###############################################################
# > runif(3)*1e8
# [1] 10515098 55871737 11522135
seedsMain <- data.frame(
  "patches"  = 10515098,
  "times"    = 55871737,
  "dynamics" = 11522135
)

runDictionary <- data.frame(
  Directories = c(
    "TSTS_Simulations_1-1_1-1_2024-02-13",
    "TSTS_Simulations_1-1_2-2_2024-02-14",
    "TSTS_Simulations_2-1_2-2_2024-02-14",
    "TSTS_Simulations_6-1_2-2_2024-02-15",
    "TSTS_Simulations_10-1_2-2_2024-02-15",
    "TSTS_Simulations_11-1_3-3_2024-02-23",
    "TSTS_Simulations_11-1_4-4_2024-02-23"
  )
)[runDictionaryChoice, ]

interventionPatchDictionary <- data.frame(
  PatchAffinities = c(
    # Detection via if string begins with a numeric or a non-numeric.
    # If numeric, it takes it as a fixed set of affinities.
    # If non-numeric, it attempts to treat the string as a function name.
    # In the latter case, it provides ONLY NumberEnvironments as an argument.

    toString(rep(0, NumberOfEnvironments)), # Patches -> {0}
    toString(rep(0.5, NumberOfEnvironments)), # Patches -> {0.5}
    toString(rep(1, NumberOfEnvironments)), # Patches -> {1}
    toString(c(rep(0, NumberOfEnvironments/2),
               rep(1, NumberOfEnvironments/2))), # Patches -> {0, 1} Gradient
    "sample.int.normalized", # Patches -> {0, 1} Unif @ Random
    "patchTypes.0.Half.1", # Patches -> {0, 0.5, 1} Gradient
    "sample.int.3", # Patches -> {0, 0.5, 1} Unif @ Random
    "runifRing", # Patches -> [0, 1] Gradient
    "runif" # Patches -> [0, 1] Unif @ Random
  )
)[interventionPatchDictionaryChoice, ]
interventionPatchSeed <- withRandom(
  runif(interventionPatchSeedChoice)[interventionPatchSeedChoice] * 1e8,
  seed = seedsMain$patches
)

interventionTimeDictionary <- data.frame(
  # Time1, Time2; called by eval(str2lang(X)) where X is the string below
  #               and "result" is the file that is loaded.
  Time1 = c(
    "median(result$Events$Times)",
    "quantile(result$Events$Times, p = 0.25)"
  ),
  Time2 = c(
    "1/2 * max(result$Events$Times)",
    "quantile(result$Events$Times, p = 0.75)"
  ),
  Method = c(# each needs a custom implementation unfortunately!
    "mean",
    "runif"
  )
)[interventionTimeDictionaryChoice, ]
interventionTimeSeed <- withRandom(
  runif(interventionTimeSeedChoice)[interventionTimeSeedChoice] * 1e8,
  seed = seedsMain$times
)

interventionDispersalDictionary <- rbind(
  data.frame(Resistance = Inf, Configuration = "None"),
  expand.grid(
    Resistance = 10^c(0:9),
    Configuration = c("Ring", "Line", "Complete")
  ))[ifelse(is.na(interventiondispersalDictionaryChoice),
            1, interventiondispersalDictionaryChoice + 2), ]

interventionDistanceDictionary <- data.frame(
  rhofunction = c( # Take patch
    "rho.2.0.1.euclidean",
    "rho.2.1.2.euclidean"
  )
)[interventionDistanceDictionaryChoice, ]

# Files: ######################################################################

appendID <- paste0(
  # PARAMETERS:
  interventionPatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
  interventionTimeDictionaryChoice, "-",
  interventionDistanceDictionaryChoice, "-", # Sometimes want to change.
  interventionDispersalDictionaryChoice, "-"
  , "_",
  # SEEDS:
  interventionPatchSeedChoice, "-",
  interventionTimeSeedChoice
)

# Note runDictionary is only 1 column, so it provides only a singleton.
datfiles <- dir(path = runDictionary, pattern = "_Simulation_", full.names = T)
