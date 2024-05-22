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

poolpatchDictionaryChoice <-
  # 17 # 100 Species, 2 Patches, Species {0, 1} Alternating, Patches {0.5}.
  # 18 # 200 Species, 2 Patches, Species {0, 1} Alternating, Patches {0.5}.
  19 # 100 Species, 10 Patches, Species {0, 1} Alternating, Patches {0.5}.
  # 20 # 200 Species, 10 Patches, Species {0, 1} Alternating, Patches {0.5}.

poolpatchSeedChoice <-
  # 1 # Used on 2024-02-13
  # 2 # Used on 2024-02-14, 2024-02-15
  # 3 # Used on 2024-02-23
  # 4 # Used on 2024-02-23 for the 2 patch system
  # 5 # Used on 2024-03-08 for choice 17 above.
  # 6 # Used on 2024-03-08 for choice 18 above.
  7 # Used on 2024-03-08 for choice 19 above.
  # 8 # Used on 2024-03-08 for choice 20 above.

dynamicsDictionaryChoice <-
  1 # Law and Morton 1996, Size-Structured Lotka-Volterra, Default Parameters
#
dynamicsSeedChoice <-
  # 1 # Used on 2024-02-13
  # 2 # Used on 2024-02-14
  # 3 # Used on 2024-02-23
  # 4 # Used on 2024-02-23 for the 2 patch system
  # 5 # Used on 2024-03-08 for choice 17 above.
  # 6 # Used on 2024-03-08 for choice 18 above.
  7 # Used on 2024-03-08 for choice 19 above.
  # 8 # Used on 2024-03-08 for choice 20 above.

eventsDictionaryChoice <-
  #   Multipliers:
  # 1 # Immigration: 1, Extirpation: 1, NumberOfEvents: 1 # Default
  2 # Immigration: 1, Extirpation: 1, NumberOfEvents: 2 # For Interventions.
eventsSeedChoice <-
  # 1 # Used on 2024-02-13
  # 2 # Used on 2024-02-14 for both 1-1 and 2-1.
  # 3 # Used on 2024-02-23
  # 4 # Used on 2024-02-23 for the 2 patch system
  # 5 # Used on 2024-03-08 for ppDChoice 17 above.
  # 6 # Used on 2024-03-08 for ppDChoice 18 above.
  7 # Used on 2024-03-08 for ppDChoice 19 above.
  # 8 # Used on 2024-03-08 for ppDChoice 20 above.

dispersalDictionaryChoice <-
  15 # c(NA, 5, 0)
  # Index: Ones place is resistance to Dispersal on a log scale.
  #      : Tens place is configuration: 0* = Ring, 1* = Line, 2* = Complete.
  #      : Special: NA corresponds to no dispersal.
  # Note : No randomness, so we don't need a seed.

# choose r' = r * rho ^ (sign(r)), but what rho?
distanceDictionaryChoice <- # for m, n in [0, 1], rho(m, n) = ...
  # 1 # 2 ^ (- euclid(m, n)) => rho in [1/2, 1] for 1-D
  2 # 2 ^ (1 - 2 euclid(m, n)) => rho in [1/2, 2] for 1-D
  # 3 # 10 ^ (1 - 2 euclid(m, n)) => rho in [1/10, 10] for 1-D

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

# Why? so we can have single argument functions with partials.
repFixed <- function(value = 0.5) {
  force(value)
  function(n) {rep(value, n)}
}
rep_0 <- repFixed(0)
rep_0.5 <- repFixed()
rep_1 <- repFixed(1)

evensplit <- function(values = c(0, 1)) {
  force(values)
  function(n) {
    c(rep(values, times = floor(n / length(values))),
      if (n %% length(values) != 0) {
        values[1:(n %% length(values))]
      })
  }
}
evensplit_01 <- evensplit()

gradientline_01 <- function(n) {
  c(rep(0, ceiling(n/2)), rep(1, floor(n/2)))
}
gradientline_0half1 <- function(n) {
  left <- rep(0, floor(n / 3)); right <- rep(1, floor(n/3))
  return(c(left, rep(0.5, n - length(left) - length(right)), right))
}

# Dictionaries: ###############################################################
# > runif(3)*1e8
# [1] 21622193 73825470 83066253
seedsMain <- data.frame(
  "pools"    = 21622193,
  "events"   = 73825470,
  "dynamics" = 83066253
)

poolpatchDictionary <- expand.grid(
  BasalConsumerRatio = 1/2,
  NSpecies = c(100, 200),
  PoolFunction = "RMTRCode2::LawMorton1996_species",
  PoolParameters = c(
    paste("Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
          "LogBodySize = c(-2, -1, -1, 0)", sep = "; ")
  ),
  PoolDispersalSpeed = 1, # Value divided by DispersalResistance to get current.
  NumberEnvironments = c(2, 10),
  SpeciesAffinities = c(
    # Pool with {0.5} affinities.
    "rep_0.5",
    # 2  # Pool {0, 1} patch affinities at random.
    "sample.int.normalized",
    # 3  # Pool {0, 0.5, 1} patch affinities at random.
    "sample.int.3",
    # 4  # Pool [0, 1] patch affinities at random.
    "runif",
    # 5  # Pool {0, 1} alternating affinities.
    "evensplit_01"
  ),
  PatchAffinities = c(
    # Detection via if string begins with a numeric or a non-numeric.
    # If numeric, it takes it as a fixed set of affinities.
    # If non-numeric, it attempts to treat the string as a function name.
    # In the latter case, it provides ONLY NumberEnvironments as an argument.
    "rep_0.5", #             Patch {0.5} affinities.
    "gradientline_01", #     Patch {0, 1} affinities. Gradient Line.
    "evensplit_01", #        Patch {0, 1} affinities. Alternating.
    "gradientline_0half1", # Patch {0, 0.5, 1} affinities. Gradient Line.
    "patchTypes.0.Half.1", # Patch {0, 0.5, 1} affinities. Gradient Ring.
    "runifRing" #            Patch [0, 1] affinities. Gradient Ring at Random.
  ),
  stringsAsFactors = FALSE
)[poolpatchDictionaryChoice, ] %>% dplyr::mutate(
  Basals = ceiling((1 - (1 + BasalConsumerRatio)^(-1)) * NSpecies),
  Consumers = NSpecies - Basals
)
poolpatchSeed <- withRandom(
  runif(poolpatchSeedChoice)[poolpatchSeedChoice] * 1e8,
  seed = seedsMain$pools
  )

dynamicsDictionary <- data.frame(
  InteractionFunction = "RMTRCode2::LawMorton1996_CommunityMat",
  InteractionParameters = "Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
  DynamicsFunction = "RMTRCode2::PerCapitaDynamics_Type1",
  stringsAsFactors = FALSE
)[dynamicsDictionaryChoice, ]
dynamicsSeed <- withRandom(
  runif(dynamicsSeedChoice)[dynamicsSeedChoice] * 1e8,
  seed = seedsMain$dynamics
)

eventsDictionary <- expand.grid(
  EventsFunction = "defaultEvents", # Takes Number of Environments and Species.
  EventsNumberMultiplier = c(1, 2),
  ImmigrationMultiplier = 1,
  ImmigrationFunction = "RMTRCode2::ArrivalFUN_Example2",
  ExtirpationMultiplier = 1, # Frequency multiplier
  ExtirpationFunction = "RMTRCode2::ExtinctFUN_Example2",
  ExtirpationProportion = c(1, 0.9, 0), # Proportion of population removed.
  stringsAsFactors = FALSE
)[eventsDictionaryChoice, ]
eventsSeed <- withRandom(
  runif(eventsSeedChoice)[eventsSeedChoice] * 1e8,
  seed = seedsMain$events
)

dispersalDictionary <- rbind(
  data.frame(Resistance = Inf, Configuration = "None"),
  expand.grid(
    Resistance = 10^c(0:9),
    Configuration = c("Ring", "Line", "Complete"),
    stringsAsFactors = FALSE
  ))[ifelse(is.na(dispersalDictionaryChoice),
            1, dispersalDictionaryChoice + 2), ]

distanceDictionary <- data.frame(
  rhofunction = c( # Take patch
    "rho.2.0.1.euclidean",
    "rho.2.1.2.euclidean",
    "rho.10.1.2.euclidean",
    stringsAsFactors = FALSE
  )
)[distanceDictionaryChoice, ]

# Files: ######################################################################
partialID <- paste0(
  # PARAMETERS:
  poolpatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
  dynamicsDictionaryChoice# , # "-",
  # eventsDictionaryChoice, "-", # Sometimes want to change.
  # dispersalDictionaryChoice, "-", # Commonly changed.
  # distanceDictionaryChoice, #"-", # Sometimes changed.
  , "_",
  # SEEDS:
  poolpatchSeedChoice, "-",
  dynamicsSeedChoice
)
datfolder <- file.path(
  directory,
  paste0(
    "TSTS_Simulations_", # Separate the Name (TSTS) and Type (Simulations)
    partialID,
    "_", # Separate the Date
    Sys.Date())
)
if (!dir.exists(datfolder)) {
  dir.create(datfolder)
}

fullID <- paste0(
  # PARAMETERS:
  poolpatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
  dynamicsDictionaryChoice, "-",
  eventsDictionaryChoice, "-", # Sometimes want to change.
  dispersalDictionaryChoice, "-", # Commonly changed.
  distanceDictionaryChoice #"-", # Sometimes changed.
  , "_",
  # SEEDS:
  poolpatchSeedChoice, "-",
  dynamicsSeedChoice, "-",
  eventsSeedChoice
)

datfile <- file.path(
  datfolder,
  paste0("TSTS_Simulation_", fullID, ".RData")
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
  PatchAffinities <- datfile_ppd_envir$PatchAffinities
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
         poolpatchSeed, dynamicsSeed, ID = partialID,
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

if (is.function(rprime)) {
  # Calculate rprime using Parms$Patch
  if (is.function(InteractionMatrices$Mats[[1]])) {
    # Calculate and combine interaction matrices on the fly.
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
  }
  else {
    # Just combine the interaction matrices.
    PerCapitaDynamics <- DynamicsFunction(
      rprime,
      Matrix::bdiag(InteractionMatrices$Mats),
      poolpatchDictionary$NumberEnvironments
    )
  }
} else {
  # Treat rprime as constant and explicitly calculated.
  if (is.function(InteractionMatrices$Mats[[1]])) {
    # Calculate and combine interaction matrices on the fly.
    PerCapitaDynamics <- DynamicsFunction(
      rprime,
      function(t, y, parms) {
        Matrix::bdiag(lapply(
          InteractionMatrices$Mats,
          function(matfunc) {matfunc(t, y, parms)}
        ))
      }
    )
  }
  else {
    # Just combine the interaction matrices.
    PerCapitaDynamics <- DynamicsFunction(
      rprime,
      Matrix::bdiag(InteractionMatrices$Mats)
    )
  }
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
    if (Configuration == "Complete") {
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
  ExtinctionProportion = eventsDictionary$ExtirpationProportion,
  MaximumTimeStep = MaximumTimeStep,
  BetweenEventSteps = BetweenEventSteps,
  Verbose = TRUE,
  # Using the ellipsis pass through feature:
  Timescale = "Simulation",
  ID = fullID,
  Affinity = list(
    PatchAffinities = PatchAffinities,
    EffectiveReproductionRate = rprime
  )
)

# Save Simulation: ############################################################
save(result,
     file = datfile
)
