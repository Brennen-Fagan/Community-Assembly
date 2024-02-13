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

# Easy to convert to cargs <- as.numeric(commandArgs(TRUE)) for parallel.

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
poolpatchSeedChoice <- 1

dynamicsDictionaryChoice <-
  1 # Law and Morton 1996, Size-Structured Lotka-Volterra, Default Parameters
#
dynamicsSeedChoice <- 1

eventsDictionaryChoice <-
  #   Multipliers:
  # 1 # Immigration: 1, Extirpation: 1, NumberOfEvents: 1 # Default
  2 # Immigration: 1, Extirpation: 1, NumberOfEvents: 2 # For Interventions.
eventsSeedChoice <- 1

dispersalDictionaryChoice <-
  NA # c(NA, 0, 5)
  # Index: Ones place is resistance to Dispersal on a log scale.
  #      : Tens place is configuration: 0* = Ring, 1* = Line, 2* = Complete.
  #      : Special: NA corresponds to no dispersal.
  # Note : No randomness, so we don't need a seed.

# choose r' = r * rho ^ (sign(r)), but what rho?
distanceDictionaryChoice <- # for m, n in [0, 1], rho(m, n) = ...
  1 # 2 ^ (- euclid(m, n)) => rho in [1/2, 1] for 1-D
  # 2 # 2 ^ (1 - 2 euclid(m, n)) => rho in [1/2, 2] for 1-D

## Other Parameters: ##########################################################
EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007
MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 10 # Number of steps to reach next event to smooth.


directory <- "." # Should be "VariantExperiments"

loadPoolPatchDynamicsIfAble <- TRUE

# Libraries: ##################################################################
library(RMTRCode2)
library(dplyr)
library(Matrix)

# Functions: ##################################################################

source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
# Defines: retrieveFunction.

# We need a function to generate control random values.
# https://stackoverflow.com/a/59875367 Gwang-Jin Kim and
# https://stackoverflow.com/a/14324316 Romain Francois (same question).
withRandom <- function(expr, seed) {
  if (exists(".Random.seed")) {
    oldSeed <- .Random.seed
    on.exit({.Random.seed <<- oldSeed})
  }
  set.seed(seed)
  expr
}

# Run runif and organise in a smooth-ish ring.
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

# Discrete niche samplers.
sample.int.normalized <- function(n, slots = 2) {
  (sample.int(slots, size = n, replace = TRUE) - 1) / (slots - 1)
}
sample.int.3 <- purrr::partial(sample.int.normalized, slots = 3)

# Coupon Collector's Problem
# I think this is probably higher accuracy than the previous version.
defaultEvents <- function(
  NumberOfEnvironments, NumberOfSpecies, constant = 3
  ) {
  ceiling(
    NumberOfEnvironments * NumberOfSpecies * (
      log(NumberOfEnvironments * NumberOfSpecies) + constant
    )
  )
}

rhofunction <- function(
  base = 2, offset = 0, multiplier = 1, metric = "euclidean"
) {
  force(base);force(offset);force(multiplier)
  function(m, n) {
    base ^ (offset - multiplier * dist(
      matrix(c(m, n), byrow = TRUE, nrow = 2), method = metric)
    )
  }
}

rho.2.0.1.euclidean <- rhofunction()
rho.2.1.2.euclidean <- rhofunction(2, 1, 2)



# Dictionaries: ###############################################################
# > runif(3)*1e8
# [1] 21622193 73825470 83066253
seedsMain <- data.frame(
  "pools"    = 21622193,
  "events"   = 73825470,
  "dynamics" = 83066253
)

poolpatchDictionary <- data.frame(
  Basals = 34,
  Consumers = 66,
  PoolFunction = "RMTRCode2::LawMorton1996_species",
  PoolParameters = c(
    paste("Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
          "LogBodySize = c(-2, -1, -1, 0)", sep = "; ")
  ),
  PoolDispersalSpeed = 1, # Value divided by DispersalResistance to get current.
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
    toString(c(rep(0, 5), rep(1, 5))), #      Patch {0, 1} affinities.
    toString(c(rep(0, 5), rep(1, 5))), #      Patch {0, 1} affinities.
    toString(c(rep(0, 5), rep(1, 5))), #      Patch {0, 1} affinities.
    "0, 0, 0, 0.5, 0.5, 1, 1, 1, 0.5, 0.5", # Patch {0, 0.5, 1} affinities.
    "0, 0, 0, 0.5, 0.5, 1, 1, 1, 0.5, 0.5", # Patch {0, 0.5, 1} affinities.
    "0, 0, 0, 0.5, 0.5, 1, 1, 1, 0.5, 0.5", # Patch {0, 0.5, 1} affinities.
    "runifRing", #                            Patch [0, 1] affinities.
    "runifRing", #                            Patch [0, 1] affinities.
    "runifRing"  #                            Patch [0, 1] affinities.
  )
)[poolpatchDictionaryChoice, ]
poolpatchSeed <- withRandom(
  runif(poolpatchSeedChoice)[poolpatchSeedChoice] * 1e8,
  seed = seedsMain$pools
  )

dynamicsDictionary <- data.frame(
  InteractionFunction = "RMTRCode2::LawMorton1996_CommunityMat",
  InteractionParameters = "Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
  DynamicsFunction = "RMTRCode2::PerCapitaDynamics_Type1"
)[dynamicsDictionaryChoice, ]
dynamicsSeed <- withRandom(
  runif(dynamicsSeedChoice)[dynamicsSeedChoice] * 1e8,
  seed = seedsMain$dynamics
)

eventsDictionary <- data.frame(
  ImmigrationMultiplier = 1,
  ImmigrationFunction = "RMTRCode2::ArrivalFUN_Example2",
  ExtirpationMultiplier = 1,
  ExtirpationFunction = "RMTRCode2::ExtinctFUN_Example2",
  ExtirpationPercentage = 1,
  EventsFunction = "defaultEvents", # Takes Number of Environments and Species.
  EventsNumberMultiplier = c(1, 2)
)[eventsDictionaryChoice, ]
eventsSeed <- withRandom(
  runif(eventsSeedChoice)[eventsSeedChoice] * 1e8,
  seed = seedsMain$events
)

dispersalDictionary <- rbind(
  data.frame(Resistance = Inf, Configuration = "None"),
  expand.grid(
    Resistance = 10^c(0:9),
    Configuration = c("Ring", "Line", "Complete")
  ))[ifelse(is.na(dispersalDictionaryChoice),
            1, dispersalDictionaryChoice + 2), ]

distanceDictionary <- data.frame(
  rhofunction = c( # Take patch
    "rho.2.0.1.euclidean",
    "rho.2.1.2.euclidean"
  )
)[distanceDictionaryChoice, ]

# Files: ######################################################################
datfolder <- file.path(
  directory,
  paste0(
    "TSTS_Simulations_", # Separate the Name (TSTS) and Type (Simulations)
    poolpatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
    dynamicsDictionaryChoice, # "-",
    # eventsDictionaryChoice, "-", # Sometimes want to change.
    # dispersalDictionaryChoice, "-", # Commonly changed.
    # distanceDictionaryChoice, #"-", # Sometimes changed.
    "_", # Separate the Date
    Sys.Date())
)
if (!dir.exists(datfolder)) {
  dir.create(datfolder)
}

datfile <- file.path(
  datfolder,
  paste0(
    "TSTS_Simulation_",
    poolpatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
    dynamicsDictionaryChoice, "-",
    eventsDictionaryChoice, "-", # Sometimes want to change.
    dispersalDictionaryChoice, "-", # Commonly changed.
    distanceDictionaryChoice, #"-", # Sometimes changed.
    ".RData"
  )
)
if (file.exists(datfile)) {
  warning(paste(datfile, "already exists and will be overwritten."))
}

datfile_ppd <- file.path(
  datfolder,
  paste0(
    "TSTS_PoolPatchDynamics_",
    poolpatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
    dynamicsDictionaryChoice, #"-",
    ".RData"
  )
)
if (loadPoolPatchDynamicsIfAble) {
  if (file.exists(datfile_ppd)) {
    datfile_ppd_envir <- attach(datfile_ppd)
  } else {
    datfile_ppd_write <- TRUE
  }
}

# Pools and Interaction Matrices: #############################################
if (exists("datfile_ppd_envir")) {
  # Pool, InteractionMatrices, DynamicsFunction, CharacteristicRate,
  Pool <- datfile_ppd_envir$Pool
  InteractionMatrices <- datfile_ppd_envir$InteractionMatrices
  DynamicsFunction <- datfile_ppd_envir$DynamicsFunction
  CharacteristicRate <- datfile_ppd_envir$CharacteristicRate
} else {
  Pool <- with(poolpatchDictionary, {
    do.call(what = retrieveFunction(PoolFunction),
            # TODO: the ; handling seems excessive, and was not suggested.
          args = c(callMeMaybe2(as.list(strsplit(PoolParameters, "; ")[[1]])),
                   "Basal" = Basals,
                   "Consumer" = Consumers,
                   "seed" = withRandom(runif(1)[1] * 1e8, seed = poolpatchSeed)
          )
    )
  })

  Affinities <- with(poolpatchDictionary, {
    if(!is.na(as.numeric(substr(SpeciesAffinities, 1, 1)))) {
      # Treat as numbers
      as.numeric(unlist(strsplit(SpeciesAffinities, split = ", ")))
    } else {
      # Treat as function
      withRandom(
        retrieveFunction(SpeciesAffinities)(Basals + Consumers),
        seed = withRandom(runif(2)[2] * 1e8, seed = poolpatchSeed)
      )
    }
  })

  Speeds <- with(poolpatchDictionary, {
    if(is.numeric(PoolDispersalSpeed)) {
      rep(PoolDispersalSpeed, Basals + Consumers)
    } else if(!is.na(as.numeric(substr(PoolDispersalSpeed, 1, 1)))) {
      # Treat as numbers
      as.numeric(unlist(strsplit(PoolDispersalSpeed, split = ", ")))
    } else {
      # Treat as function
      withRandom(
        retrieveFunction(PoolDispersalSpeed)(Basals + Consumers),
        seed = withRandom(runif(3)[3] * 1e8, seed = poolpatchSeed)
      )
    }
  })

  Pool <- cbind(
    Pool,
    Speed = Speeds,
    Affinity = Affinities
  )

  PatchAffinities <- matrix(with(poolpatchDictionary, {
    if(is.numeric(PatchAffinities)) {
      rep(PatchAffinities, Basals + Consumers)
    } else if(!is.na(as.numeric(substr(PatchAffinities, 1, 1)))) {
      # Treat as numbers
      as.numeric(unlist(strsplit(PatchAffinities, split = ", ")))
    } else {
      # Treat as function
      withRandom(
        retrieveFunction(PatchAffinities)(NumberEnvironments),
        seed = withRandom(runif(4)[4] * 1e8, seed = poolpatchSeed)
      )
    }
  }), nrow = poolpatchDictionary$NumberEnvironments)

  # NOT THE FINAL PERCAPITADYNAMICS FUNCTION.
  DynamicsFunction <- with(dynamicsDictionary, {
    withRandom(
      retrieveFunction(DynamicsFunction),
      seed = withRandom(runif(1)[1] * 1e8, seed = dynamicsSeed)
    )}
  )

  IntMatFunc <- with(dynamicsDictionary, {
    withRandom(
      purrr::partial(
        retrieveFunction(InteractionFunction),
        ... =, # otherwise, partialised arguments occur first.
        !!!callMeMaybe2(InteractionParameters)
        # !!! suggested by Bing Chat. Not obvious @me how it evaluates the list.
        # Prompt 'R purrr, using a list of arguments to partialize a function'.
        # 2024/01/19
      ),
      seed = withRandom(runif(2)[2] * 1e8, seed = dynamicsSeed)
    )}
  )

  InteractionMatrices <- RMTRCode2::CreateEnvironmentInteractions(
    Pool = Pool,
    NumEnvironments = poolpatchDictionary$NumberEnvironments,
    ComputeInteractionMatrix = IntMatFunc,
    EnvironmentSeeds = withRandom(runif(3)[3] * 1e8, seed = dynamicsSeed)
  )

  CharacteristicRate <- max(unlist(lapply(
    InteractionMatrices$Mats, function(m) {abs(eigen(m)$values)}
  )))

  if (exists("datfile_ppd_write") && datfile_ppd_write) {
    save(Pool, PatchAffinities,
         InteractionMatrices, DynamicsFunction, CharacteristicRate,
         poolpatchSeed, dynamicsSeed,
         file = datfile_ppd)
  }
}

# Events: #####################################################################
EventsEach <- with(poolpatchDictionary, {
  retrieveFunction(eventsDictionary$EventsFunction)(
    NumberEnvironments, Basals + Consumers
  )
  })

Events <- with(eventsDictionary, {
  RMTRCode2::CreateAssemblySequence(
    Species = with(poolpatchDictionary, Basals + Consumers),
    NumEnvironments = poolpatchDictionary$NumberEnvironments,
    ArrivalEvents = EventsEach * ImmigrationMultiplier * EventsNumberMultiplier,
    ArrivalRate = CharacteristicRate * ImmigrationMultiplier,
    ArrivalFUN = retrieveFunction(ImmigrationFunction),
    ExtinctEvents = EventsEach * ExtirpationMultiplier * EventsNumberMultiplier,
    ExtinctRate = CharacteristicRate * ExtirpationMultiplier,
    ExtinctFUN = retrieveFunction(ExtirpationFunction),
    HistorySeed = eventsSeed
  )}
)

print(combinations <-
        table(Events$Events$Species,
              Events$Events$Environment,
              Events$Events$Type))
if(any(combinations == 0)) {warning(
  "Exists a species which doesn't immigrate to an environment."
)}

# Instantiate PerCapitaDynamics: ##############################################
# We've already built the "functional" interactions matrix, but we now need
# to apply the transformation r' = r rho^(sign(r)) and then combine.

#TODO Be careful here. If we add columns to distanceDictionary, we need to
# adjust the call here (since a 1x1 dataframe is returned as just the entry).
# Not a function, so we don't need to on the fly.
# We'll be a bit lazy here, but hopefully readable and clear.
grid <- expand.grid(
  pool = 1:nrow(Pool), # Fastest Varying
  patch = 1:length(PatchAffinities) # Slower Varying.
)
rprime <-
  rep(Pool$ReproductionRate, poolpatchDictionary$NumberEnvironments) *
  mapply(
    grid$pool,
    grid$patch,
    FUN = function(i, j) {
      retrieveFunction(distanceDictionary)(
        Pool[i, grepl("Affinity", colnames(Pool), fixed = TRUE)],
        PatchAffinities[j, ] # Forced to be a matrix.
      )[1]
    }
  ) ^ sign(rep(Pool$ReproductionRate, poolpatchDictionary$NumberEnvironments))

if (is.function(InteractionMatrices$Mats[[1]])) {
  PerCapitaDynamics <- DynamicsFunction(
    rprime,
    function(t, y, parms) {
      Matrix::bdiag(lapply(
        InteractionMatrices$Mats,
        function(matfunc) {matfunc(t, y, parms)}
      ))
      },
    poolpatchDictionary$NumberEnvironments
  )
} else {
  # Treat as constant matrices
  PerCapitaDynamics <- DynamicsFunction(
    rprime,
    Matrix::bdiag(InteractionMatrices$Mats),
    poolpatchDictionary$NumberEnvironments
  )
}

# Instantiate Dispersal Matrix: ###############################################

DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
  EnvironmentDistances = with(c(
    dispersalDictionary,
    Environments = poolpatchDictionary$NumberEnvironments
    ), {
    if (Configuration == "None") {
      DistanceMatrix <- Matrix::sparseMatrix(
        i = Environments, j = Environments, x = 0)
    }
    if (Configuration == "Ring" || Configuration == "Line")
      DistanceMatrix <- Matrix::bandSparse(
        Environments, k = c(-1, 1),
        diagonals = list(rep(Resistance, Environments - 1),
                         rep(Resistance, Environments - 1))
      )
    if (Configuration == "Ring") {
      DistanceMatrix[Environments, 1] <- Resistance
      DistanceMatrix[1, Environments] <- Resistance
    }
    if (Configuration == "Full") {
      DistanceMatrix <- matrix(Resistance,
                               nrow = Environments,
                               ncol = Environments)
      diag(DistanceMatrix) <- 0
    }
    DistanceMatrix
  }
  ),
  SpeciesSpeeds = Pool$Speed
)

# Run Simulation: #############################################################
result <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
  Pool = Pool,
  NumEnvironments = poolpatchDictionary$NumberEnvironments,
  CharacteristicRate = CharacteristicRate,
  Events = Events,
  PerCapitaDynamics = PerCapitaDynamics,
  DispersalMatrix = DispersalMatrix,
  EliminationThreshold = EliminationThreshold,
  ArrivalDensity = ArrivalDensity,
  ExtinctionProportion = eventsDictionary$ExtirpationPercentage,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = TRUE,
  # Using the ellipsis pass through feature:
  Timescale = "Simulation",
  Affinity = list(
    PatchAffinities = PatchAffinities,
    EffectiveReproductionRate = rprime
  )
)

# Save Simulation: ############################################################
save(result,
     file = datfile
)
