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
  files <- handlerID(ID)
  runDictionary <- files$runDictionary
  datfile <- files$datfile
  
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