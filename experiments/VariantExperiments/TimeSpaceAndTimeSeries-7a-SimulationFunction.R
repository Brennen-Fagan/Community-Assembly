# Libraries: ##################################################################
library(RMTRCode2)
library(dplyr)
library(Matrix)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Dictionaries.R"))

# Function Definition: ########################################################
simulationWrapper <- function(
  poolpatchDictionaryChoice,
  poolpatchSeedChoice,
  dynamicsDictionaryChoice,
  dynamicsSeedChoice,
  eventsDictionaryChoice,
  eventsSeedChoice,
  dispersalDictionaryChoice,
  distanceDictionaryChoice,
  parameters = list(), # Borrowing from stats::optim's control template
  loadPoolPatchDynamicsIfAble = TRUE,
  seedsMain = data.frame(
    "pools"    = 21622193,
    "events"   = 73825470,
    "dynamics" = 83066253
  ),
  returnResults = FALSE,
  saveResults = TRUE,
  skipIfSaveExists = TRUE, # Precedence over the next argument.
  errorIfSaveExists = FALSE
) {
  #  Parameters: ##############################################################
  params <- list(
    EliminationThreshold = 10^-4, # Below which species are removed.
    ArrivalDensity =  4 * 10 ^ 3, # Traill et al. 2007: MULTIPLIER ON ABOVE.
    MaximumTimeStep = 1, # Maximum time solver can proceed without elimination.
    BetweenEventSteps = 10 # Number of steps to reach next event to smooth.
  )
  namespold <- names(params)
  params[(namespnew <- names(parameters))] <- parameters
  if (length(namespNotFound <- namespnew[!namespnew %in% namespold])) {
    warning(
      "Unknown names in parameters: ", paste(namespNotFound, collapse = ", ")
    )
  }

  # Dictionaries: #############################################################
  poolpatchDictionary <-
    poolpatchDictionaryOrigin[poolpatchDictionaryChoice, ] %>% dplyr::mutate(
      Basals = ceiling((1 - (1 + BasalConsumerRatio)^(-1)) * NSpecies),
      Consumers = NSpecies - Basals
    )
  poolpatchSeed <- withRandom(
    runif(poolpatchSeedChoice)[poolpatchSeedChoice] * 1e8,
    seed = seedsMain$pools
  )

  dynamicsDictionary <-
    dynamicsDictionaryOrigin[dynamicsDictionaryChoice, ]
  dynamicsSeed <- withRandom(
    runif(dynamicsSeedChoice)[dynamicsSeedChoice] * 1e8,
    seed = seedsMain$dynamics
  )

  eventsDictionary <-
    eventsDictionaryOrigin[eventsDictionaryChoice, ]
  eventsSeed <- withRandom(
    runif(eventsSeedChoice)[eventsSeedChoice] * 1e8,
    seed = seedsMain$events
  )

  dispersalDictionary <-
    dispersalDictionaryOrigin[ifelse(is.na(dispersalDictionaryChoice),
                                     1, dispersalDictionaryChoice + 2), ]

  distanceDictionary <-
    distanceDictionaryOrigin[distanceDictionaryChoice, ]

  # Files: ######################################################################
  partialID <- paste0(
    # PARAMETERS:
    poolpatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
    dynamicsDictionaryChoice, "_",
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
    distanceDictionaryChoice, "_",
    # SEEDS:
    poolpatchSeedChoice, "-",
    dynamicsSeedChoice, "-",
    eventsSeedChoice
  )


  if(saveResults) {
    datfile <- file.path(
      datfolder,
      paste0("TSTS_Simulation_", fullID, ".RData")
    )
    if (file.exists(datfile)) {
      if (skipIfSaveExists) {
        warning(paste(datfile, "already exists. Returning NULL."))
        return(-1)
      } else if (errorIfSaveExists) {
        stop(paste(datfile, "already exists. Throwing error."))
      } else {
        warning(paste(datfile, "already exists and will be overwritten."))
      }
    }
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

  # Pools and Interaction Matrices: ###########################################
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
              args = c(
                callMeMaybe2(as.list(strsplit(PoolParameters, "; ")[[1]])),
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
    "Exists a species which doesn't immigrate to/extirpate from an environment."
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
        retrieveFunction(distanceDictionary$rhofunction)(
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
  if (poolpatchDictionary$NumberEnvironments > 1) {
    DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
      EnvironmentDistances = convertDispersalDictToDistMatrix(
        dispersalDictionary,
        nEnv = poolpatchDictionary$NumberEnvironments
      ),
      SpeciesSpeeds = Pool$Speed
    )
  } else {
    DispersalMatrix <- Matrix::sparseMatrix(
      i = {}, j = {}, # From documentation
      dims = c(nrow(Pool), nrow(Pool))
      )
  }

  # Run Simulation: ###########################################################
  result <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
    Pool = Pool,
    NumEnvironments = poolpatchDictionary$NumberEnvironments,
    CharacteristicRate = CharacteristicRate,
    Events = Events,
    PerCapitaDynamics = PerCapitaDynamics,
    DispersalMatrix = DispersalMatrix,
    EliminationThreshold = params$EliminationThreshold,
    ArrivalDensity = params$ArrivalDensity,
    ExtinctionProportion = eventsDictionary$ExtirpationProportion,
    MaximumTimeStep = params$MaximumTimeStep,
    BetweenEventSteps = params$BetweenEventSteps,
    Verbose = TRUE,
    # Using the ellipsis pass through feature:
    Timescale = "Simulation",
    ID = fullID,
    Affinity = list(
      PatchAffinities = PatchAffinities,
      EffectiveReproductionRate = rprime
    )
  )

  # Save Simulation: ##########################################################
  if(saveResults)
    save(result, file = datfile)

  if(returnResults) {
    return(results)
  } else {
    return(0)
  }
}
