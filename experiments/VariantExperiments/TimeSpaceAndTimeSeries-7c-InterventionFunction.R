# Libraries: ##################################################################
library(RMTRCode2)
library(dplyr)
library(Matrix)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Dictionaries.R"))

# Function Definition: ########################################################
handlerID <- function(
  ID # list containing: Tag (optional) +
  # poolpatchDictionaryChoice, poolpatchSeedChoice, dynamicsDictionaryChoice,
  # dynamicsSeedChoice, eventsDictionaryChoice, eventsSeedChoice,
  # dispersalDictionaryChoice, distanceDictionaryChoice, Date (YYYY-MM-DD)
) {
  # ID Handling: ##############################################################
  stopmsg <- "ID not interpreted correctly. Try providing explicit names."
  if (!is.null(names(ID))) { # List/data.frame of components.
    if (toupper("Tag") %in% toupper(names(ID))) {
      runDictionaryTag <- ID[which(toupper("Tag") == toupper(names(ID)))]
    } else {
      runDictionaryTag <- "TSTS_Simulations"
    }

    runDictionary <- with(ID, {
      paste0(runDictionaryTag, "_",
             # PARAMETERS:
             poolpatchDictionaryChoice, "-", dynamicsDictionaryChoice, "_",
             # SEEDS:
             poolpatchSeedChoice, "-", dynamicsSeedChoice, "_",
             Date)
    })

    datfile <- with(
      ID, paste0(
        runDictionaryTag, "_",
        # PARAMETERS:
        poolpatchDictionaryChoice, "-", dynamicsDictionaryChoice, "-",
        eventsDictionaryChoice, "-", dispersalDictionaryChoice, "-",
        distanceDictionaryChoice, "_",
        # SEEDS:
        poolpatchSeedChoice, "-", dynamicsSeedChoice, "-", eventsSeedChoice,
        ".RData")
    )

  } else if (length(ID) == 1) {# The file.path(dirname, basename) itself.
    if (is.list(ID)) {
      ID <- unlist(ID)
    }
    runDictionary <- dirname(ID)
    datfile <- basename(ID)

  } else if (length(ID) <= 10 && length(ID) >= 8) { # Unnamed list of components
    runDictionary0 <- 0
    runDictionaryDate <- NULL
    if (length(ID) > 8) {
      # Maybe Date
      runDictionaryDateTest <- tryCatch(
        {
          runDictionaryDate <- as.Date(ID[[length(ID)]],
                                       tryFormats = c("%Y-%m-%d"));
          TRUE
        },
        error = function(e) return(FALSE))
      if(length(ID) == 10 && !runDictionaryDateTest) {
        stop("ID has 10 entries, but last not parsed as %Y-%m-%d date.")
      }
      if (length(ID) == 10 || !runDictionaryDateTest) {
        # == 9 entries and last is not a date (X)or have 10 entries.
        runDictionaryTag <- ID[[1]]
        runDictionary0 <- 1
      } else {
        runDictionaryTag <- "TSTS_Simulations"
      }
    }

    tryCatch({
      runDictionary <- # Note awkward order.
        paste0(runDictionaryTag, "_",
               # PARAMETERS:
               # poolpatchDictionaryChoice, "-", dynamicsDictionaryChoice, "_",
               ID[[runDictionary0 + 1]], "-", ID[[runDictionary0 + 3]], "_",
               # SEEDS:
               # poolpatchSeedChoice, "-", dynamicsSeedChoice, "_",
               ID[[runDictionary0 + 2]], "-", ID[[runDictionary0 + 4]], "_",
               runDictionaryDate)

      if (is.null(runDictionaryDate)) {
        if (!exists("directory")) directory <- "."
        candidates <- grep(dir(directory), pattern = runDictionary,
                           value = TRUE)
        if (length(candidates) != 1) {
          stop("Provided ID ambiguous. Try providing explicit names.")
        }
        runDictionary <- candidates
      }

      datfile <- with(
        ID, paste0(
          runDictionaryTag, "_",
          # PARAMETERS:
          # poolpatchDictionaryChoice, "-", dynamicsDictionaryChoice, "-",
          ID[[runDictionary0 + 1]], "-", ID[[runDictionary0 + 3]], "_",
          # eventsDictionaryChoice, "-", dispersalDictionaryChoice, "-",
          ID[[runDictionary0 + 5]], "-", ID[[runDictionary0 + 7]], "_",
          # distanceDictionaryChoice, "_",
          ID[[runDictionary0 + 8]], "_",
          # SEEDS:
          # poolpatchSeedChoice, "-", dynamicsSeedChoice, "-", eventsSeedChoice,
          ID[[runDictionary0 + 2]], "-", ID[[runDictionary0 + 4]], "-",
          ID[[runDictionary0 + 6]],
          ".RData")
      )

    }, error = function(e) {
      stop(stopmsg)
    }
    )
  } else {
    stop(stopmsg)
  }

  return(list(
    runDictionary = runDictionary, # dirname
    datfile = file.path(runDictionary, datfile) # file.path(dirname, basename)
  ))
}

interventionWrapper <- function(
  ID, # See handlerID.
  interventionPatchDictionaryChoice,
  interventionPatchSeedChoice,
  interventionTimeDictionaryChoice,
  interventionTimeSeedChoice,
  interventionDispersalDictionaryChoice,
  interventionDistanceDictionaryChoice,
  seedsMain = data.frame(
    "patches"  = 10515098,
    "times"    = 55871737,
    "dynamics" = 11522135
  ),
  returnResults = FALSE,
  saveResults = TRUE,
  skipIfSaveExists = TRUE, # Precedence over the next argument.
  errorIfSaveExists = FALSE
) {
  files <- handlerID(ID)
  runDictionary <- files$runDictionary
  datfile <- files$datfile

  # Source Files: #############################################################

  # Extract exposed properties.
  datfile_properties <- strsplit(
    strsplit(basename(datfile), split = ".", fixed = TRUE)[[1]][1],
    split = "_", fixed = TRUE)
  stopifnot(length(datfile_properties) == 1)

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

  # Put failure point before loads.
  if(saveResults) {
    # Reformulate file name.
    filename <- file.path(
      dirname(datfile),
      paste0(datfile_properties[[1]][1],
             "_Intervention_",
             datfile_properties[[1]][3],"_", datfile_properties[[1]][4],
             "_", appendID,
             if (length(datfile_properties[[1]]) > 4)
               paste0("_", datfile_properties[[1]][
                 4:length(datfile_properties[[1]])
                 ],
                 collapse = ""),
             ".RData")
    )

    if (file.exists(filename)) {
      if (skipIfSaveExists) {
        warning(paste(filename, "already exists. Returning NULL."))
        return(-1)
      } else if (errorIfSaveExists) {
        stop(paste(filename, "already exists. Throwing error."))
      } else {
        warning(paste(filename, "already exists and will be overwritten."))
      }
    }
  }

  datPoolMats <- dir(runDictionary,
                     "PoolPatchDynamics.+[.]RData$",
                     full.names = T)
  poolMats <- new.env()
  load(datPoolMats, envir = poolMats)
  NumberOfEnvironments <- length(poolMats$InteractionMatrices$Mats)

  # Load file.
  loaded <- load(datfile) # names
  stopifnot(length(loaded) == 1)
  loaded <- get(loaded) # objects

  # Regular function: #########################################################

  interventionPatchDictionary <-
    interventionPatchDictionaryOrigin[
      interventionPatchDictionaryChoice, , drop = FALSE
      ]
  interventionPatchSeed <- withRandom(
    runif(interventionPatchSeedChoice)[interventionPatchSeedChoice] * 1e8,
    seed = seedsMain$patches
  )

  interventionTimeDictionary <- interventionTimeDictionaryOrigin[
    interventionTimeDictionaryChoice,
    ]
  interventionTimeSeed <- withRandom(
    runif(interventionTimeSeedChoice)[interventionTimeSeedChoice] * 1e8,
    seed = seedsMain$times
  )

  if (toupper(interventionDispersalDictionaryChoice) %in%
      c("P", "PRE", "PREV", "PREVIOUS") # "PREVIOUS" is meaning. "p" for storage.
  ) {
    interventionDispersalDictionaryChoice <- as.numeric(strsplit(
      datfile_properties[[1]][3], split = '-'
    )[[1]][4])
  }
  interventionDispersalDictionary <- dispersalDictionaryOrigin[
    ifelse(is.na(interventionDispersalDictionaryChoice),
           1, interventionDispersalDictionaryChoice + 2),
    ]

  if (toupper(interventionDistanceDictionaryChoice) %in%
      c("P", "PRE", "PREV", "PREVIOUS") # "PREVIOUS" is meaning. "p" for storage.
  ) {
    interventionDistanceDictionaryChoice <- as.numeric(strsplit(
      datfile_properties[[1]][3], split = '-'
    )[[1]][5])
  }
  interventionDistanceDictionary <-
    distanceDictionaryOrigin[interventionDistanceDictionaryChoice, ]

  # Instantiate Dispersal Matrix: #############################################
  if (NumberOfEnvironments > 1) {
    DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
      EnvironmentDistances = convertDispersalDictToDistMatrix(
        interventionDispersalDictionary,
        nEnv = NumberOfEnvironments
      ),
      SpeciesSpeeds = poolMats$Pool$Speed
    )
  } else {
    DispersalMatrix <- Matrix::sparseMatrix(
      i = {}, j = {}, # From documentation
      dims = c(nrow(poolMats$Pool), nrow(poolMats$Pool))
    )
  }

  # Interventions: ############################################################
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
  rho <- retrieveFunction(interventionDistanceDictionary$rhofunction)

  grid <- expand.grid(
    pool = 1:nrow(poolMats$Pool), # Fastest Varying
    patch = 1:length(interventionPatchAffinities) # Slower Varying.
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

  # Need to specify each patch separately with how I've implemented the
  # DynamicsFunction. (I recall wanting to have per patch stats.)
  # Note, this is inefficient.
  if (interventionTimeDictionary$InterventionTimespan == 0) {
    rprimeSwitches <- lapply(1:NumberOfEnvironments, function(i) {
      if (i %in% interventionPatches) {
        switchMatrices(
          loaded$Ellipsis$Affinity$EffectiveReproductionRate[
            (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
            ],
          rprime[
            (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
            ],
          switchtime = interventionTime
        )
      } else {
        function(t,...) rprime[
          (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
          ]
      }
    })
    rprimef <- function(t, parms, ...) {
      return(rprimeSwitches[[parms$Patch]](t, parms, ...))
    }
  } else {
    rprimeSwitches <- lapply(1:NumberOfEnvironments, function(i) {
      if (i %in% interventionPatches) {
        interpolateMatrices(
          loaded$Ellipsis$Affinity$EffectiveReproductionRate[
            (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
            ],
          rprime[
            (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
            ],
          switchtime = interventionTime,
          timespan = interventionTimeDictionary$InterventionTimespan
        )
      } else {
        function(t,...) rprime[
          (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
          ]
      }
    })
    rprimef <- function(t, parms, ...) {
      return(rprimeSwitches[[parms$Patch]](t, parms, ...))
    }
  }

  interventionPerCapitaDynamics <- with(poolMats, {
    # TECH DEBT: Copied from 6a-simulations.R
    if (is.function(rprimef)) {
      # Calculate rprimef using Parms$Patch
      if (is.function(InteractionMatrices$Mats[[1]])) {
        # Calculate and combine interaction matrices on the fly.
        DynamicsFunction(
          rprimef,
          function(t, y, parms) {
            Matrix::bdiag(lapply(
              InteractionMatrices$Mats,
              function(matfunc) {matfunc(t, y, parms)}
            ))
          },
          NumberOfEnvironments
        )
      }
      else {
        # Just combine the interaction matrices.
        DynamicsFunction(
          rprimef,
          Matrix::bdiag(InteractionMatrices$Mats),
          NumberOfEnvironments
        )
      }
    } else {
      # Treat rprimef as constant and explicitly calculated.
      if (is.function(InteractionMatrices$Mats[[1]])) {
        # Calculate and combine interaction matrices on the fly.
        DynamicsFunction(
          rprimef,
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
          rprimef,
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
      Times > timeInitial
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
    ParentRun = x,
    ID = paste0(loaded$Ellipsis$ID, "_", appendID),
    Affinity = list(
      PatchAffinitiesOld = loaded$Ellipsis$Affinity$PatchAffinities,
      PatchAffinitiesIntervention = interventionPatchAffinities,
      PatchInterventions = interventionPatches,
      TimeIntervention = interventionTime,
      EffectiveReproductionRateOld =
        loaded$Ellipsis$Affinity$EffectiveReproductionRate,
      EffectiveReproductionRateIntervention = rprime
    )
  )

  # Save Simulation: ##########################################################
  if(saveResults)
    save(result, file = filename)

  if(returnResults) {
    return(results)
  } else {
    return(0)
  }
}
