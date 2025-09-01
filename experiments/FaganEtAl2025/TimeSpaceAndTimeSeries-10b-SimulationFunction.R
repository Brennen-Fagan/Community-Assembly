# Define a wrapper function to more easily keep track of the inputs from the
# dictionary combinations.

# Libraries: ##################################################################
library(RMTRCode2)
library(dplyr)
library(Matrix)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "FaganEtAl2025"
source(file.path(directory, "TimeSpaceAndTimeSeries-10-Dictionaries.R"))

# Function Definition: ########################################################
simulationWrapper <- function(
  poolpatchDictionaryChoice,
  poolpatchSeedChoice,
  dynamicsDictionaryChoice,
  dynamicsSeedChoice,
  eventsDictionaryChoice,
  eventsSeedChoice,
  initialConditionsDictionaryChoice,
  initialConditionsSeedChoice,
  dispersalDictionaryChoice,
  distanceDictionaryChoice,
  speciesAffinityDictionaryChoice,
  speciesAffinitySeedChoice,
  patchAffinityDictionaryChoice,
  patchAffinitySeedChoice,
  parameters = list(), # Borrowing from stats::optim's control template
  loadPoolPatchDynamicsIfAble = TRUE,
  seedsMain = data.frame(
    # > runif(6)*1e8
    # [1] 45980437 17222244 64343803  3361566 79927625  8554899
    "pools"             = 45980437,
    "events"            = 17222244,
    "dynamics"          = 64343803,
    "speciesAffinity"   = 3361566,
    "patchAffinity"     = 79927625,
    "initialConditions" = 8554899
  ),
  logisticCarryingCapacity = NULL, # list(Basal=x, Consumer=y) or (Total=z)
  returnResults = FALSE,
  saveResults = TRUE,
  skipIfSaveExists = TRUE, # Precedence over the next argument.
  errorIfSaveExists = FALSE
) {
  #  Parameters: ##############################################################
  # Intelligent defaults and downstream parameter handling.
  params <- list(
    EliminationThreshold = 10^-4, # Below which species are removed.
    ArrivalDensity =  4 * 10 ^ 3 * 10^-4, # Traill et al. 2007: MULTIPLIED BY ET.
    MaximumTimeStep = 5, # Maximum time solver can proceed without elimination.
    BetweenEventSteps = 2, # Number of steps to reach next event to smooth.
    Date = Sys.Date() # Label for the folders
  )
  namespold <- names(params)
  params[(namespnew <- names(parameters))] <- parameters
  if (length(namespNotFound <- namespnew[!namespnew %in% namespold])) {
    warning(
      "Unknown names in parameters: ", paste(namespNotFound, collapse = ", ")
    )
  }

  # Dictionaries: #############################################################
  # Match the provided entry values to the actual dictionary entries.
  # Initialise the random number generator factories while we are at it.
  poolpatchDictionary <-
    poolpatchDictionaryOrigin[poolpatchDictionaryChoice, ] %>% dplyr::mutate(
      Basals = ceiling((1 - (1 + BasalConsumerRatio)^(-1)) * NSpecies),
      Consumers = NSpecies - Basals
    )
  poolpatchSeed <- withRandom(
    runif(poolpatchSeedChoice)[poolpatchSeedChoice] * 1e8,
    seed = seedsMain$pools
  )
  poolpatchSeedIndex <- indexFactory()

  dynamicsDictionary <-
    dynamicsDictionaryOrigin[dynamicsDictionaryChoice, ]
  dynamicsSeed <- withRandom(
    runif(dynamicsSeedChoice)[dynamicsSeedChoice] * 1e8,
    seed = seedsMain$dynamics
  )
  dynamicsSeedIndex <- indexFactory()

  eventsDictionary <-
    eventsDictionaryOrigin[eventsDictionaryChoice, ]
  eventsSeed <- withRandom(
    runif(eventsSeedChoice)[eventsSeedChoice] * 1e8,
    seed = seedsMain$events
  )
  eventsSeedIndex <- indexFactory()

  initialConditionsDictionary <-
    initialConditionsDictionaryOrigin[initialConditionsDictionaryChoice,]
  initialConditionsSeed <- withRandom(
    runif(initialConditionsSeedChoice)[initialConditionsSeedChoice] * 1e8,
    seed = seedsMain$initialConditions
  )
  initialConditionsSeedIndex <- indexFactory()

  dispersalDictionary <-
    dispersalDictionaryOrigin[ifelse(is.na(dispersalDictionaryChoice),
                                     1, dispersalDictionaryChoice + 2), ]

  distanceDictionary <-
    distanceDictionaryOrigin[distanceDictionaryChoice, ]

  speciesAffinityDictionary <-
    speciesAffinityDictionaryOrigin[speciesAffinityDictionaryChoice, ]
  speciesAffinitySeed <- withRandom(
    runif(speciesAffinitySeedChoice)[speciesAffinitySeedChoice] * 1e8,
    seed = seedsMain$speciesAffinity
  )
  speciesAffinitySeedIndex <- indexFactory()

  patchAffinityDictionary <-
    patchAffinityDictionaryOrigin[patchAffinityDictionaryChoice, ]
  patchAffinitySeed <- withRandom(
    runif(patchAffinitySeedChoice)[patchAffinitySeedChoice] * 1e8,
    seed = seedsMain$patchAffinity
  )
  patchAffinitySeedIndex <- indexFactory()

  # Edit params: ##############################################################
  # Brief the user if they sent in values that overwrite the old defaults
  # without invoking the parameters arguments.
  if ("InteractionEliminationThreshold" %in% names(dynamicsDictionary)) {
    warning(paste0("Overwriting Elimination Threshold with: ",
                   dynamicsDictionary$InteractionEliminationThreshold))
    params$EliminationThreshold <-
      dynamicsDictionary$InteractionEliminationThreshold
  }
  if ("ColonizationPropaguleSize" %in% names(eventsDictionary)) {
    warning(paste0("Overwriting Elimination Threshold with: ",
                   eventsDictionary$ColonizationPropaguleSize))
    params$ArrivalDensity <-
      eventsDictionary$ColonizationPropaguleSize
  }

  # Files: ####################################################################
  # Organise the file names and folders to make tracking of outputs easier.
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
      params$Date)
  )
  if (!dir.exists(datfolder)) {
    dir.create(datfolder)
  }

  fullID <- paste0(
    # PARAMETERS:
    poolpatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
    dynamicsDictionaryChoice, "-",
    eventsDictionaryChoice, "-", # Sometimes want to change.
    initialConditionsDictionaryChoice, "-",
    dispersalDictionaryChoice, "-", # Commonly changed.
    distanceDictionaryChoice, "-", # Commonly changed.
    speciesAffinityDictionaryChoice, "-", # Of Experimental interest.
    patchAffinityDictionaryChoice, "_",
    # SEEDS:
    poolpatchSeedChoice, "-",
    dynamicsSeedChoice, "-",
    eventsSeedChoice, "-",
    initialConditionsSeedChoice, "-",
    speciesAffinitySeedChoice, "-",
    patchAffinitySeedChoice
  )

  # Make sure that the call has the intended consequences.
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

  # Save some time and processing power for creating the shared resources.
  # (This seems to default to parallel and uncapped on research servers???)
  datfile_ppd <- file.path(
    datfolder,
    paste0(
      "TSTS_PoolPatchDynamics_", partialID, ".RData"
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
  # Safely load or create and save the shared resources.
  if (exists("datfile_ppd_envir")) {
    # Pool, InteractionMatrices, DynamicsFunction, CharacteristicRate,
    Pool <- datfile_ppd_envir$Pool
    InteractionMatrices <- datfile_ppd_envir$InteractionMatrices
    DynamicsFunction <- datfile_ppd_envir$DynamicsFunction
    CharacteristicRate <- datfile_ppd_envir$CharacteristicRate
  } else {
    # For each: set a random seed, establish a temporary environment, and run
    # the appropriate function.
    id <- poolpatchSeedIndex()
    Pool <- with(poolpatchDictionary, {
      do.call(what = retrieveFunction(PoolFunction),
              args = c(
                "Basal" = Basals,
                "Consumer" = Consumers,
                "Parameters" = callMeMaybe2(paste0(
                  "c(",
                  toString(c(PoolK1InteractionEffectiveness,
                             PoolK2ConsumerSizeAdvantage,
                             PoolK3ConsumerPredationRange,
                             PoolK4ConsumerEfficiency,
                             PoolK5BasalBiomass,
                             PoolK6CoefOfVariation)),
                  ")"
                )),
                "LogBodySize" = callMeMaybe2(paste0(
                  "c(", PoolBasalLogBodySize, ",", PoolConsumerLogBodySize, ")"
                )),
                "seed" = withRandom(runif(id)[id] * 1e8, seed = poolpatchSeed)
              )
      )
    })

    id <- poolpatchSeedIndex()
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
          seed = withRandom(runif(id)[id] * 1e8, seed = poolpatchSeed)
        )
      }
    })

    Pool <- cbind(
      Pool,
      Speed = Speeds
    )

    id <- dynamicsSeedIndex()
    # NOT THE FINAL PERCAPITADYNAMICS FUNCTION.
    DynamicsFunction <- with(dynamicsDictionary, {
      withRandom(
        retrieveFunction(DynamicsFunction),
        seed = withRandom(runif(id)[id] * 1e8, seed = dynamicsSeed)
      )}
    )

    id <- dynamicsSeedIndex()
    IntMatFunc <- with(dynamicsDictionary, {
      withRandom(
        purrr::partial(
          retrieveFunction(InteractionFunction),
          ... =, # otherwise, partialised arguments occur first.
          "Parameters" = c(InteractionK1InteractionEffectiveness,
                           InteractionK2ConsumerSizeAdvantage,
                           InteractionK3ConsumerPredationRange,
                           InteractionK4ConsumerEfficiency,
                           InteractionK5BasalBiomass,
                           InteractionK6CoefOfVariation)
        ),
        seed = withRandom(runif(id)[id] * 1e8, seed = dynamicsSeed)
      )}
    )

    id <- dynamicsSeedIndex()
    InteractionMatrices <- RMTRCode2::CreateEnvironmentInteractions(
      Pool = Pool,
      NumEnvironments = poolpatchDictionary$NumberEnvironments,
      ComputeInteractionMatrix = IntMatFunc,
      EnvironmentSeeds = withRandom(runif(id)[id] * 1e8, seed = dynamicsSeed)
    )

    CharacteristicRate <- max(unlist(lapply(
      InteractionMatrices$Mats, function(m) {abs(eigen(m)$values)}
    )))

    if (exists("datfile_ppd_write") && datfile_ppd_write) {
      save(Pool,
           InteractionMatrices, DynamicsFunction, CharacteristicRate,
           poolpatchSeed, dynamicsSeed, ID = partialID,
           file = datfile_ppd)
    }
  }

  # Affinities are not a part of the pool-patch framework.
  # Instead, we assign afterwards to allow us to vary them independently.
  id <- speciesAffinitySeedIndex()
  SpeciesAffinities <- with(speciesAffinityDictionary, {
    if(!is.na(as.numeric(substr(SpeciesAffinities, 1, 1)))) {
      # Treat as numbers
      as.numeric(unlist(strsplit(SpeciesAffinities, split = ", ")))
    } else {
      # Treat as function
      withRandom(
        retrieveFunction(SpeciesAffinities)(nrow(Pool)),
        seed = withRandom(runif(id)[id] * 1e8, seed = speciesAffinitySeed)
      )
    }
  })

  Pool <- cbind(Pool, Affinity = SpeciesAffinities)

  id <- patchAffinitySeedIndex()
  PatchAffinities <- matrix(with(patchAffinityDictionary, {
    if(is.numeric(PatchAffinities)) {
      rep(PatchAffinities, Basals + Consumers)
    } else if(!is.na(as.numeric(substr(PatchAffinities, 1, 1)))) {
      # Treat as numbers
      as.numeric(unlist(strsplit(PatchAffinities, split = ", ")))
    } else {
      # Treat as function
      withRandom(
        retrieveFunction(PatchAffinities)(poolpatchDictionary$NumberEnvironments),
        seed = withRandom(runif(id)[id] * 1e8, seed = patchAffinitySeed)
      )
    }
  }), nrow = poolpatchDictionary$NumberEnvironments)

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

  # Instantiate PerCapitaDynamics: ############################################
  # We've already built the "functional" interactions matrix, but we now need
  # to apply the transformation r' = r rho^(sign(r)) and then combine.

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

  if (!is.null(logisticCarryingCapacity)) {
    if ("basal" %in% tolower(names(logisticCarryingCapacity))) {
      if ("consumer" %in% tolower(names(logisticCarryingCapacity))) {
        # Both Basal and Consumer, but separately. See below for descriptions.
        stop("logisticCarryingCapacity term not implemented.")
      } else {
        # Only Basal
        rprimeMax <- rprime
        basalVec <- (Pool$Type == "Basal")
        sizeVec <- Pool$Size * basalVec
        sizeVecLen <- length(sizeVec)
        rprime <- function(t, y, parms, ...) {
          #TODO but sizes are in the pool, and we need only for this patch, and we need only basals.
          # Automatically evaluated per patch (parms$Patch == i, unlist-lapply)
          # RMTRCode2::PerCapitaDynamics_Type1, but the whole y is provided.
          rprimeMax * (
            1 - basalVec * sum(
              y[1:sizeVecLen + sizeVecLen*(parms$Patch - 1)] * sizeVec
            ) / logisticCarryingCapacity$Basal
          )
        }
      }
    } else if ("consumer" %in% tolower(names(logisticCarryingCapacity))) {
      # Only Consumer
      # NOTE: Not clear about implementation if r(prime) is negative.
      #       Presumably, it would instead act on the consumption term.
      stop("logisticCarryingCapacity term not implemented.")
    } else if ("total" %in% tolower(names(logisticCarryingCapacity))) {
      # Both Basal and Consumer, together.
      # NOTE: Implementing the consumer side of this is not obvious.
      #       For Basals, this is a sum over all species, not only Basals.
      stop("logisticCarryingCapacity term not implemented.")
    } else {
      stop("logisticCarryingCapacity term not recognised.")
    }
  }

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

  # Instantiate Dispersal Matrix: #############################################
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

  # Set initial population conditions: ########################################
  # For each environment, decide what to instantiate, apply the method, then
  # supply the remainders (within the environment for correct order) with 0s.
  # Following calculations assume that the pool is ordered basal first.
  id <- initialConditionsSeedIndex()
  popInitial <-
    unlist(lapply(
      1:poolpatchDictionary$NumberEnvironments,
      function(i) {
        pI <- rep(0, nrow(Pool))

        pICalc <- switch(
          initialConditionsDictionary$Species,
          "None" = NA,
          "Basal" = 1:sum(Pool$Type == "Basal"),
          "All" = 1:nrow(Pool))

        if (!is.na(pICalc[1])) {
          pI[pICalc] <- switch(
            initialConditionsDictionary$Method,
            "Solve" = solve(InteractionMatrices$Mats[[i]][pICalc, pICalc],
                            -rprime[pICalc]), # TODO rprime might be a function!
            "InverseSize" = 1/Pool$Size[pICalc],
            "Set" = rep(initialConditionsDictionary$Argument, length(pICalc)),
            "Random" = withRandom(
              runif(length(pICalc),
                    min = params$EliminationThreshold,
                    max = initialConditionsDictionary$Argument),
              seed = withRandom(runif(id)[id] * 1e8,
                                seed = initialConditionsSeed)
            )
          )
        }
        return(pI)
      }))

  # Run Simulation: ###########################################################
  result <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
    Pool = Pool,
    PopulationInitial = popInitial,
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
    Verbose = FALSE,
    # Using the ellipsis pass through feature:
    Timescale = "Simulation",
    ID = fullID,
    Affinity = list(
      SpeciesAffinities = SpeciesAffinities,
      PatchAffinities = PatchAffinities,
      EffectiveReproductionRate = rprime
    )#,
    # Can't not be passed through without a separate call that doesn't pass it.
    # LogisticCarryingCapacity = logisticCarryingCapacity
    # But it is implicitly stored by the function inside
    # Ellipsis$Affinity$EffectiveReproductionRate, i.e.,
    # environment(result$Ellipsis$Affinity$EffectiveReproductionRate)$logisticCarryingCapacity
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
