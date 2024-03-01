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

# While this code can be run in parallel, I'm generally disinclined.
# I've not written it to suggest mass production and would rather
# embarassingly parallel approached be used before proper foreach.
# Concerns include memory and time commitment.
# Foreach is used instead in order to facilitate coding and consistency.
cores <- 1

overwriteOutput <- TRUE

# Parameters: #################################################################
# Note that many of our options here vary between a deterministic mode
# and a stochastic mode. As such, we assign seeds still, but may not use them.

interventionPatchDictionaryChoice <-
  # 1 # Random 50% Patches -> {0}
  # 2 # Random 50% Patches -> {0.5}
  # 3 # Random 50% Patches -> {1}
  # 4 # Random 50% Patches -> {0, 1} Gradient
  # 5 # Random 50% Patches -> {0, 1} Unif @ Random
  # 6 # Random 50% Patches -> {0, 0.5, 1} Gradient Ring
  # 7 # Random 50% Patches -> {0, 0.5, 1} Unif @ Random
  # 8 # Random 50% Patches -> [0, 1] Gradient Ring
  # 9 # Random 50% Patches -> [0, 1] Unif @ Random
  10 # Last 50% Patches -> {0}
  # 11 # Last 50% Patches -> {0.5}
  # 12 # Last 50% Patches -> {1}
  # 13 # Last 50% Patches -> {0, 1} Gradient
  # 14 # Last 50% Patches -> {0, 1} Unif @ Random
  # 15 # Last 50% Patches -> {0, 0.5, 1} Gradient Ring
  # 16 # Last 50% Patches -> {0, 0.5, 1} Unif @ Random
  # 17 # Last 50% Patches -> [0, 1] Gradient Ring
  # 18 # Last 50% Patches -> [0, 1] Unif @ Random
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


library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

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

# We'll pre-load the pool and patch dynamics. This allows us to infer some
# parameters.
datPoolMats <- dir(runDictionary,
                   "PoolPatchDynamics.+[.]RData$",
                   full.names = T)
poolMats <- new.env()
load(datPoolMats, envir = poolMats)
NumberOfEnvironments <- length(poolMats$InteractionMatrices$Mats)

interventionPatchDictionary <- expand.grid(
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
    "patchTypes.0.Half.1", # Patches -> {0, 0.5, 1} Gradient Ring
    "sample.int.3", # Patches -> {0, 0.5, 1} Unif @ Random
    "runifRing", # Patches -> [0, 1] Gradient Ring
    "runif" # Patches -> [0, 1] Unif @ Random
  ),
  InterventionPercentage = c(
    0.5
  ),
  InterventionLocation = c(
    # Percentage seems easiest
    NA, # == Random
    1, # Last
    0 # First
  ),
  stringsAsFactors = FALSE
)[interventionPatchDictionaryChoice, , drop = FALSE]
interventionPatchSeed <- withRandom(
  runif(interventionPatchSeedChoice)[interventionPatchSeedChoice] * 1e8,
  seed = seedsMain$patches
)

interventionTimeDictionary <- data.frame(
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
)[interventionTimeDictionaryChoice, ]
interventionTimeSeed <- withRandom(
  runif(interventionTimeSeedChoice)[interventionTimeSeedChoice] * 1e8,
  seed = seedsMain$times
)

# TECH DEBT: Copied from 6a-simulations.R. Should be a common resource in final.
# Furthermore, we may want to add an option to use the old dispersal;
# many of my previous runs were the same pool-patches but different dispersals.
# This implementation forces them to experience a new dispersal, and without a
# proper transition (I don't think I have a dynamic dispersal function?).
interventionDispersalDictionary <- rbind(
  data.frame(Resistance = Inf, Configuration = "None"),
  expand.grid(
    Resistance = 10^c(0:9),
    Configuration = c("Ring", "Line", "Complete")
  ))[ifelse(is.na(interventionDispersalDictionaryChoice),
            1, interventionDispersalDictionaryChoice + 2), ]

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
  # Where dynamics would go if necessary.
  interventionTimeDictionaryChoice, "-",
  interventionDispersalDictionaryChoice, "-", # Sometimes want to change.
  interventionDistanceDictionaryChoice
  , "_",
  # SEEDS:
  interventionPatchSeedChoice, "-",
  # Where dynamics would go if necessary.
  interventionTimeSeedChoice
)

# Note runDictionary is only 1 column, so it provides only a singleton.
datfiles <- dir(path = runDictionary,
                pattern = "Simulation.+[.]RData$",
                full.names = T)

# Instantiate Dispersal Matrix: ###############################################
# Copied from 6a-Simulations.R. Should the configuration "switch" be a function
# in order to make it a common resource as well?
DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
  EnvironmentDistances = with(c(
    interventionDispersalDictionary,
    Environments = NumberOfEnvironments
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

# Interventions: ##############################################################
interventionSuccess <- foreach::foreach(
  x = iterators::iter(datfiles)
) %op% {
  # Extract exposed properties.
  x_properties <- strsplit(
    strsplit(basename(x), split = ".", fixed = TRUE)[[1]][1],
    split = "_", fixed = TRUE)
  stopifnot(length(x_properties) == 1)

  # Reformulate file name.
  filename <- file.path(
    dirname(x),
    paste0(x_properties[[1]][1],
           "_Intervention_",
           x_properties[[1]][3],"_", x_properties[[1]][4],
           "_", appendID,
           if (length(x_properties[[1]]) > 4)
             paste0("_", x_properties[[1]][4:length(x_properties[[1]])],
                    collapse = ""),
           ".RData")
  )

  # Do we re-run/skip if the file already exists?
  if (file.exists(filename)) {
    if (overwriteOutput) {
      warning(paste(filename, "already exists and will be overwritten."))
    } else {
      return(1)
    }
  }

  # Load file.
  loaded <- load(x) # names
  stopifnot(length(loaded) == 1)
  loaded <- (get(loaded)) # objects

  ### Determine Interventions: ################################################
  # Number:
  NumberOfInterventions <- round(
    NumberOfEnvironments * interventionPatchDictionary$InterventionPercentage
  )

  # New values:
  # Note: this is technically excessive, but keeping the number of environments
  # in total is easier for programming below.
  interventionPatchAffinities <- matrix(with(interventionPatchDictionary, {
    if(is.numeric(PatchAffinities)) {
      rep(PatchAffinities, nrow(poolMats$Pool))
    } else if(!is.na(as.numeric(substr(PatchAffinities, 1, 1)))) {
      # Treat as numbers
      as.numeric(unlist(strsplit(PatchAffinities, split = ", ")))
    } else {
      # Treat as function
      withRandom(
        retrieveFunction(PatchAffinities)(NumberOfEnvironments),
        seed = withRandom(runif(1)[1] * 1e8, seed = interventionPatchSeed)
      )
    }
  }), nrow = NumberOfEnvironments)

  # Locations:
  interventionPatches <- with(interventionPatchDictionary, {
    if (is.na(InterventionLocation)) {
      # Random: Ring style (fill rightwards)
      interventionPatches <- ((withRandom(
        sample.int(NumberOfEnvironments, 1),
        seed = withRandom(runif(2)[2] * 1e8, seed = interventionPatchSeed)
      ) + 1 : ceiling(NumberOfEnvironments * InterventionPercentage)
      ) %% NumberOfEnvironments) + 1
    } else if (InterventionLocation == 0) {
      # Left: Line style (fill inwards)
      interventionPatches <-
        1 : ceiling(NumberOfEnvironments * InterventionPercentage)
    } else if (InterventionLocation == 1) {
      # Right: Line style (fill inwards)
      interventionPatches <-
        # (0.40, 0.50] of 10 means 10, 9, 8, 7, 6
        (floor(NumberOfEnvironments * (1 - InterventionPercentage)) + 1) :
        NumberOfEnvironments
    } else {
      # We choose to interpret as centred so that someone can have easily
      # understood output for 0.5 rather than skewed one direction or another.
      # Ring style (fill symmetrically), alternate left->right->left->right
      # (e.g. 0.5 with 10 patches 5 -> 4 -> 6 -> 3 -> 7
      #           while 5 patches 3 -> 2 -> 4)
      interventionPatches <-
        round((NumberOfEnvironments - 1) * InterventionLocation) + 1
      interventionPatches <- interventionPatches +
        seq(from = ceiling(-NumberOfEnvironments * InterventionPercentage / 2),
            to = floor(NumberOfEnvironments * InterventionPercentage / 2),
            by = 1)
    }
  })

  # Time:
  interventionTime <- with(interventionTimeDictionary, {
    if (Method == "mean") {
      mean(c(eval(str2lang(Time1)), eval(str2lang(Time2))))
    } else if(Method == "runif") {
      withRandom(
        runif(1, min = eval(str2lang(Time1), max = eval(str2lang(Time2)))),
        seed = withRandom(runif(1)[1] * 1e8, seed = interventionTimeSeed)
      )
    } else {
      stop(paste0("Method (Intervention Time): ", method, " not implemented."))
    }
  })

  ##### Post-intervention adjusted intrinsic growth/decay rates: ##############
  # TECH DEBT: Copied from 6a-Simulation.R.
  rho <- retrieveFunction(data.frame(
    rhofunction = c( # Take patch
      "rho.2.0.1.euclidean",
      "rho.2.1.2.euclidean"
    )
  )[as.numeric(strsplit(x_properties[[1]], "-", fixed = TRUE)[[3]][5]), ])

  grid <- expand.grid(
    pool = 1:nrow(Pool), # Fastest Varying
    patch = 1:length(PatchAffinities) # Slower Varying.
  )
  rprime <- with(poolMats, {
    mapply(
      grid$pool, # Species
      grid$patch, # Location
      FUN = function(i, j) {
        ifelse(
          j %in% interventionPatches, # if in intervention
          Pool$ReproductionRate[i] * rho( # recalculate for the new patches.
            Pool[i, grepl("Affinity", colnames(Pool), fixed = TRUE)],
            interventionPatchAffinities[j, ] # Forced to be a matrix.
          )[1]^sign(Pool$ReproductionRate[i]),
          loaded$Ellipsis$Affinity$EffectiveReproductionRate[ # else use old
            (j - 1) * nrow(Pool) + i
            ]
        )
      }
    )
  })

  if (interventionTimeDictionary$InterventionTimespan == 0) {
    rprime <- switchMatrices(
      loaded$Ellipsis$Affinity$EffectiveReproductionRate,
      rprime,
      switchtime = interventionTime
    )
  } else {
    rprime <- interpolateMatrices(
      loaded$Ellipsis$Affinity$EffectiveReproductionRate,
      rprime,
      switchtime = interventionTime,
      timespan = interventionTimeDictionary$InterventionTimespan
    )
  }

  interventionPerCapitaDynamics <- with(poolMats, {
    # TECH DEBT: Copied from 6a-simulations.R
    if (is.function(rprime)) {
      # Calculate rprime using Parms$Patch
      if (is.function(InteractionMatrices$Mats[[1]])) {
        # Calculate and combine interaction matrices on the fly.
        DynamicsFunction(
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
        DynamicsFunction(
          rprime,
          Matrix::bdiag(InteractionMatrices$Mats),
          poolpatchDictionary$NumberEnvironments
        )
      }
    } else {
      # Treat rprime as constant and explicitly calculated.
      if (is.function(InteractionMatrices$Mats[[1]])) {
        # Calculate and combine interaction matrices on the fly.
        DynamicsFunction(
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
        DynamicsFunction(
          rprime,
          Matrix::bdiag(InteractionMatrices$Mats)
        )
      }
    }
  })

  # Set new initial values based on when we fork the simulation.
  # Note: we don't need to worry about time conversions yet because all
  # calculations so far are in the same time scale as was provided.
  ## Abundance
  timeInterventionRow <- which.max(
    loaded$Abundance[, 1] > interventionTime
  ) - 1

  timeInitial <- loaded$Abundance[timeInterventionRow, 1]
  abundanceInitial <- loaded$Abundance[timeInterventionRow, -1]

  ## Events
  # Note formatting is a list containing a data.frame named Events.
  eventsPostIntervention <- list(
    Events = loaded$Events %>% dplyr::filter(
      Times > s$Simulation$Abundance[timeInterventionRow, 1]
    ))
  # Why not timeIntervention? To make sure that we don't miss out on an event.
  # Possibly unnecessary.
  eventsPostIntervention$Events$Success <- NA

  # Need to worry about time conversion here because we need to be on the
  # simulation scale if we weren't already (to match with the rates).
  # Run Simulation: ###########################################################
  result <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
    # From Intervening
    PopulationInitial = abundanceInitial,
    TimeInitial = timeInitial,
    # From PoolMats
    Pool = poolMats$Pool,
    NumEnvironments = NumberOfEnvironments,
    CharacteristicRate = poolMats$CharacteristicRate,
    # Recalculated
    Events = eventsPostIntervention,
    PerCapitaDynamics = interventionPerCapitaDynamics,
    DispersalMatrix = DispersalMatrix,
    # From Loaded
    EliminationThreshold = loaded$Parameters$EliminationThreshold,
    ArrivalDensity = loaded$Parameters$ArrivalDensity,
    ExtinctionProportion = loaded$Parameters$ExtinctionProportion,
    MaximumTimeStep = loaded$Parameters$MaximumTimeStep,
    BetweenEventSteps = loaded$Parameters$BetweenEventSteps,
    Verbose = TRUE,
    # Using the ellipsis pass through feature:
    Timescale = "Simulation",
    ID = paste0(loaded$Ellipsis$ID, "_", appendID),
    Affinity = list(
      PatchAffinitiesOld = loaded$Ellipsis$Affinity$PatchAffinities,
      PatchAffinitiesIntervention = PatchAffinities,
      PatchInterventions = interventionPatches,
      EffectiveReproductionRate = rprime
    )
  )

  # Save Simulation: ############################################################
  save(result,
       file = filename
  )

  return(0)
}

