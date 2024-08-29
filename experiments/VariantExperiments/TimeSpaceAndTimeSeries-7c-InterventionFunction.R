# Libraries: ##################################################################
library(RMTRCode2)
library(dplyr)
library(Matrix)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Dictionaries.R"))

# Function Definition: ########################################################
interventionWrapper <- function(
  ID, # list containing: Tag (optional) + 
  # poolpatchDictionaryChoice, poolpatchSeedChoice, dynamicsDictionaryChoice,
  # dynamicsSeedChoice, eventsDictionaryChoice, eventsSeedChoice,
  # dispersalDictionaryChoice, distanceDictionaryChoice, Date (YYYY-MM-DD)
  interventionPatchDictionaryChoice,
  interventionPatchSeedChoice,
  interventionTimeDictionaryChoice,
  interventionTimeSeedChoice,
  interventionDispersalDictionaryChoice,
  interventionDistanceDictionaryChoice,
  parameters = list(), # Borrowing from stats::optim's control template
  loadPoolPatchDynamicsIfAble = TRUE,
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
      runDictionary <- dirname(unlist(ID))
      datfile <- unlist(ID)
    } else {
      runDictionary <- dirname(ID)
      datfile <- ID
    }
    
  } else if (length(ID) <= 10 && length(ID) >= 8) { # Unnamed list of components
    if (length(ID) > 8) {
      # Maybe Date
      runDictionaryDateTest <- tryCatch(
        {as.Date(ID[[length(ID)]], tryFormats = c("%Y-%m-%d")); TRUE},
        error = function(e) return(FALSE))
      if(length(ID) == 10 && !runDictionaryDateTest) {
        stop("ID has 10 entries, but last not parsed as %Y-%m-%d date.")
      }
      if (length(ID) == 10 || !runDictionaryDateTest) {
        # == 9 entries and last is not a date (X)or have 10 entries.
        runDictionaryTag <- ID[[1]]
      } else {
        runDictionaryTag <- "TSTS_Simulations"
      }
    }
    
    runDictionary <- tryCatch({
      file.path(
        paste0()
      )
    }, error = function(e) {
      stop(stopmsg)
      }
    )
  } else {
    stop(stopmsg)
  }
  
  
  # Source Files: #############################################################
  
  # Extract exposed properties.
  datfile_properties <- strsplit(
    strsplit(basename(datfile), split = ".", fixed = TRUE)[[1]][1],
    split = "_", fixed = TRUE)
  stopifnot(length(datfile_properties) == 1)
  
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
  
  if (!toupper(interventionDispersalDictionaryChoice) %in%
      c("P", "PRE", "PREV", "PREVIOUS") # "PREVIOUS" is meaning. "p" for storage.
  ) {
    interventionDispersalDictionary <- dispersalDictionaryOrigin[
      ifelse(is.na(interventionDispersalDictionaryChoice),
             1, interventionDispersalDictionaryChoice + 2),
      ]
  } else {
    previousDispersalDictionaryChoice <- as.numeric(strsplit(
      datfile_properties[[1]][3], split = '-'
    )[[1]][4])
    interventionDispersalDictionary <- dispersalDictionaryOrigin[
      ifelse(is.na(previousDispersalDictionaryChoice),
             1, previousDispersalDictionaryChoice + 2), 
      ]
  }
  
  if (!toupper(interventionDistanceDictionaryChoice) %in%
      c("P", "PRE", "PREV", "PREVIOUS") # "PREVIOUS" is meaning. "p" for storage.
  ) {
    interventionDistanceDictionary <-
      distanceDictionaryOrigin[interventionDistanceDictionaryChoice, ]
  } else {
    interventionDistanceDictionaryChoice <- "p"
  }
  
  # Instantiate Dispersal Matrix: #############################################
  if (NumberOfEnvironments > 1) {
    DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
      EnvironmentDistances = convertDispersalDictToDistMatrix(
        interventionDispersalDictionary,
        nEnv = NumberOfEnvironments
      ),
      SpeciesSpeeds = Pool$Speed
    )
  } else {
    DispersalMatrix <- Matrix::sparseMatrix(
      i = {}, j = {}, # From documentation
      dims = c(nrow(Pool), nrow(Pool))
    )
  }
  
  # Interventions: ############################################################
  
}