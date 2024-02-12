# Introduction: ###############################################################
# Sequel to TimeSpaceAndTimeSeries-4a-Bootstrap.R.
# We are still introducing a simulation based intervention.
# Here, we're going to set-up some new base case simulations.
# There will be no intervention (we'll save that for the next step this time).
# That makes this a halfway house between 4a (above) and 
# MNA-Image-ExampleOutcome-Create.R, which I've been abusing for base cases.
#
# As always, please see the previous files for some design choices, 
# although we aim to improve design at each stage.

# Parameters: #################################################################
# We're going to group our simulation parameters, as in 4a, but for the purpose
# of replicating the behaviour of MNA-...-Create.R.
# Note there is NO intervention == nothing will be changing during simulation.

simulationsNumber <- 1
poolpatchDictionaryChoice <-
  1  # Pool with no patch affinity.       Patch with no affinities.
  # 2  # Pool {0, 1} patch affinities.      Patch {0, 1} affinities.
  # 3  # Pool {0, 0.5, 1} patch affinities. Patch {0, 1} affinities.
  # 4  # Pool [0, 1] patch affinities.      Patch {0, 1} affinities.
  # 5  # Pool {0, 1} patch affinities.      Patch {0, 0.5, 1} affinities.
  # 6  # Pool {0, 0.5, 1} patch affinities. Patch {0, 0.5, 1} affinities.
  # 7  # Pool [0, 1] patch affinities.      Patch {0, 0.5, 1} affinities.
  # 8  # Pool {0, 1} patch affinities.      Patch [0, 1] affinities.
  # 9  # Pool {0, 0.5, 1} patch affinities. Patch [0, 1] affinities.
  # 10 # Pool [0, 1] patch affinities.      Patch [0, 1] affinities.
poolpatchSeed <- 1

dispersalDictionaryChoice <- 
  -1
  # Index: Ones place is resistance to Dispersal
  #      : Tens place is configuration: 0* = Ring, 1* = Line, 2* = Complete.
  #      : Special: -1 corresponds to no dispersal.
  # Note : No randomness, so we don't need a seed. 
  
eventsDictionaryChoice <- 
  #   Multipliers:
  # 1 # Immigration: 1, Extirpation: 1, NumberOfEvents: 1 # Default
   2 # Immigration: 1, Extirpation: 1, NumberOfEvents: 2 # For Interventions.
eventsSeed <- 1
  
dynamicsDictionaryChoice <-
   1 # Law and Morton 1996, Size-Structured Lotka-Volterra, Default Parameters
  #
dynamicsSeed <- 1

# Libraries: ##################################################################

# Dictionaries: ###############################################################
# > runif(3)*1e8
# [1] 21622193 73825470 83066253
seedsMain <- c(
  "pools"    = 21622193,
  "events"   = 73825470,
  "dynamics" = 83066253
)

runifRing <- function(n, ...) {
  indices <- if (n %% 2) {
    # Odd (why?)
    c(1, seq(from = 2, by = 2, to = n), seq(from = n, by = -2, to = 2))
  } else {
    # Even.
    c(1, seq(from = 2, by = 2, to = n), seq(from = n - 1, by = -2, to = 2))
  }
  sort(runif(n, ...))[indices]
}
sample.int.normalized <- function(n, slots = 2) {
  (sample.int(slots, size = n, replace = TRUE) - 1) / (slots - 1)
}
sample.int.3 <- purrr::partial(sample.int.normalized, slots = 3)

poolpatchDictionary <- data.frame(
  Basals = 34,
  Consumers = 66,
  NumberEnvironments = 10,
  SpeciesAffinities = c(
    # Pool with no patch affinity.
    toString(rep(0, 100)),
    # 2  # Pool {0, 1} patch affinities.      Patch {0, 1} affinities.
    "sample.int.normalized",
    # 3  # Pool {0, 0.5, 1} patch affinities. Patch {0, 1} affinities.
    "sample.int.3",
    # 4  # Pool [0, 1] patch affinities.      Patch {0, 1} affinities.
    "runif",
    # 5  # Pool {0, 1} patch affinities.      Patch {0, 0.5, 1} affinities.
    "sample.int.normalized",
    # 6  # Pool {0, 0.5, 1} patch affinities. Patch {0, 0.5, 1} affinities.
    "sample.int.3",
    # 7  # Pool [0, 1] patch affinities.      Patch {0, 0.5, 1} affinities.
    "runif",
    # 8  # Pool {0, 1} patch affinities.      Patch [0, 1] affinities.
    "sample.int.normalized",
    # 9  # Pool {0, 0.5, 1} patch affinities. Patch [0, 1] affinities.
    "sample.int.3",
    # 10 # Pool [0, 1] patch affinities.      Patch [0, 1] affinities.
    "runif"
  ),
  PatchAffinities = c(
    # Detection via if string begins with a numeric or a non-numeric.
    # If numeric, it takes it as a fixed set of affinities.
    # If non-numeric, it attempts to treat the string as a function name.
    # In the latter case, it provides ONLY NumberEnvironments as an argument.
    toString(rep(0, 10)), #                   Patch with no affinities.
    toString(rep(0, 5), rep(1, 5)), #         Patch {0, 1} affinities.
    toString(rep(0, 5), rep(1, 5)), #         Patch {0, 1} affinities.
    toString(rep(0, 5), rep(1, 5)), #         Patch {0, 1} affinities.
    "0, 0, 0, 0.5, 0.5, 1, 1, 1, 0.5, 0.5", # Patch {0, 0.5, 1} affinities.
    "0, 0, 0, 0.5, 0.5, 1, 1, 1, 0.5, 0.5", # Patch {0, 0.5, 1} affinities.
    "0, 0, 0, 0.5, 0.5, 1, 1, 1, 0.5, 0.5", # Patch {0, 0.5, 1} affinities.
    "runifRing", #                            Patch [0, 1] affinities.
    "runifRing", #                            Patch [0, 1] affinities.
    "runifRing"  #                            Patch [0, 1] affinities.
  )
)[poolpatchDictionaryChoice, ]

dispersalDictionary <- rbind(
  data.frame(Resistance = Inf, Configuration = "None"),
  expand.grid(
    Resistance = 10^c(0:9), 
    Configuration = c("Ring", "Line", "Complete")
  ))[dispersalDictionaryChoice + 2, ]

eventsDictionary <- data.frame(
  
)[eventsDictionaryChoice, ]

dynamicsDictionary <- data.frame(
  
)[dynamicsDictionaryChoice, ]