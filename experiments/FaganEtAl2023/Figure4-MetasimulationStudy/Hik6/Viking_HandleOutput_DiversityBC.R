# Script Goals: ###############################################################
#   0. In Parallel,
#   1. Load the Data,
#   2. Extract Alpha, Beta, Gamma diversities over time,
#   3. Record Data as Present,
#   4. Determine if we can compress the Data further.
#   5. Save the Data in a more compressed format if so,
#   6. Save the biodiversity values separately, but alongside its parameters.
#   7. Determine Data that is missing,
#   8. Save list of missing data to a file.

# Parameters: #################################################################
divide_time_by <- 1E4
preferred_rows_per_event <- 1.5

outputLocation <- file.path(".", paste0("SaveDiversityBC_", Sys.Date()))
if (!dir.exists(outputLocation)) {
  dir.create(outputLocation, showWarnings = FALSE)
}

# Libraries: ##################################################################
print("Loading libraries")
librarypath <- file.path(".", "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

packages <- c(
  "dplyr",      # Data Manipulation
  "Matrix",     # Common Format
  "readr",
  "vegan",      # Biodiversity measurements
  "parallel",   # Base parallel dependency
  "doParallel", # For compatible cluster
  "foreach",    # For parallel evaluation
  "iterators"   # For parallel evaluation
)

# Does not work because it tries to use the system-wide libraries... Oops.
# update.packages(repos = 'https://cloud.r-project.org',
#                 oldPkgs = packages, ask = FALSE)

for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, lib = librarypath,
                     repos = 'https://cloud.r-project.org',
                     dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

if (!require("RMTRCode2", character.only = TRUE)) {
  install.packages(
    "RMTRCode2_0.2.tar.gz", lib = librarypath,
    repos = NULL, type = "source"
  )
}
library(RMTRCode2)#, lib.loc = librarypath) # lib.loc shouldn't be necessary.

source("HelperFunctions.R") # local override with the patch so don't rebuild.

# # Receive the Selection: ######################################################
print("Receiving number of cores:")

# cargs <- 4
cargs <- as.numeric(commandArgs(trailingOnly = TRUE))

clust <- parallel::makeCluster(cargs[1], outfile = "")
doParallel::registerDoParallel(clust)

# Files: ######################################################################
print("Identifying files.")

directories <- dir(path = ".", pattern = "SaveOutput",
                   full.names = TRUE, include.dirs = TRUE)
files <- dir(path = directories,
             pattern = "^MNA[-]Hik6.+[.]RData$",
             full.names = TRUE, recursive = TRUE,
             include.dirs = TRUE, no.. = TRUE)

# Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
filesParameters <- dir(path = "Prepared_2023-03-03", pattern = "^MNA[-].+Cases[.]csv$",
                       full.names = TRUE, recursive = TRUE,
                       include.dirs = TRUE, no.. = TRUE)
parametersMaster <- load(file.path("Prepared_2023-03-03", "MNA-Hik6Parameters.RData"))
parametersCSVs <- lapply(filesParameters, utils::read.csv)

# 0 In Parallel: ##############################################################
results <- foreach::foreach(
  file = iterators::iter(files), .packages = "dplyr"
) %dopar% {
  print(file)
  # Retrieve the trailing id numbers from before the file extension.
  # 1st: Set, 2nd: CaseNumber, 3rd: History, 4th: Part
  # Note, if last two are not present, all histories are bundled together.
  # If last two are present, each file has a single part of a single history.
  idNums <- suppressWarnings(na.omit(
    as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                             split = "-", fixed = TRUE)[[1]],
                    n = 4))),
    classes = "simpleWarning")

  print(idNums)

  fileName <- file.path(
    outputLocation,
    paste0("MNA-Hik6", LETTERS[idNums[1]], "-Cases-DiversityBC-",
           paste(idNums, collapse = "-"), ".RData")
  )

  if (file.exists(fileName)) {
    print("Already done.")

    if (length(idNums) == 2){
      return(c(idNums, NA, NA))
    } else {
      return(idNums)
    }
  }

  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    fileContents <- get(fileContents)

    print("loaded")

    # 2 Extract Alpha, Beta, Gamma diversities over time: ######################
    if (length(idNums) == 2) { # fileContents is a list of results
      Diversity <- list(
        Diversities = lapply(
          fileContents,
          FUN = thinAndCalculateDiversities,
          nspecies = systemBase$species,
          preferred_rows_per_event = preferred_rows_per_event
        ))
    } else if (length(idNums) == 4) {# fileContents is a result
      Diversity <- list(
        Diversities = list(thinAndCalculateDiversities(
          fileContents,
          nspecies = systemBase$species,
          preferred_rows_per_event = preferred_rows_per_event
        ))
      )
    } else stop("idNums is not 2 or 4.")
    print("diversity")
    # 6. Save the biodiversity values...: #####################################
    # separately, but alongside its parameters.
    # idNums indicates which variation we are dealing with, as well as the
    # associated parameters.
    # This is split between MNA-MasterParameters.RData and the .csv files.
    theParameters <- parametersCSVs[[idNums[1]]][((idNums[2] - 1) %/% 10) + 1,]
    Diversity$PoolMod <- systemMods$PoolBodySizes[theParameters$Framework, ]
    Diversity$NoiseMod <- systemMods$InteractionsNoiseMultiplier[theParameters$Framework]
    Diversity$NeutralMod<- systemMods$NeutralRateMultipliers[theParameters$Neutral,]
    Diversity$SpaceMod <- systemMods$SpaceDistanceMultiplier[theParameters$Space]

    save(Diversity,
         file = fileName
    )
    print("save")

    # Ran into memory issue. We'll see if this fixes it.
    rm(Diversity)
    rm(fileContents)
    gc()

    ###   3. Record Data as Present: ##########################################
    # This is already in the file list, but we can make it easier by taking the
    # id numbers. If there is an error at any point, record the relevant id
    # numbers and the error to know that something went wrong.
    if (length(idNums) == 2){
      return(c(idNums, NA, NA))
    } else {
      return(idNums)
    }
  }, error = function(e, id = idNums) {
    return(list(ID = id,
                error = e))
  })
}

# Cleanup: ####################################################################
parallel::stopCluster(clust)

# 7. Determine Data that is missing: ##########################################
results_iflist <- unlist(lapply(results, is.list))
results_success <- do.call(rbind, results[!results_iflist])
results_error <- do.call(rbind,
                         lapply(results[results_iflist], function(x) x$ID))
results_present <- dplyr::bind_rows(
  data.frame(results_success),
  data.frame(results_error)
)
names(results_present) <- c("Set", "CaseNumber", "History", "Part")
results_present <- results_present %>% dplyr::arrange(Set, CaseNumber)
results_missing <- data.frame(
  Set = rep(1:3, unlist(lapply(parametersCSVs[1:3], nrow)) * 10),
  CaseNumber = unlist(lapply(
    unlist(lapply(parametersCSVs[1:3], nrow)) * 10,
    seq, from = 1, by = 1
  ))
)
results_missing <- dplyr::anti_join(
  results_missing, results_present, by = c("Set", "CaseNumber")
)

# 8. Save list of missing data to a file: #####################################
write.csv(results_missing,
          file = file.path(outputLocation,
                           paste0("MNA-Hik6MissingRuns.csv"))
)
write.csv(results_error,
          file = file.path(outputLocation,
                           paste0("MNA-Hik6ErrorRuns.csv"))
)
