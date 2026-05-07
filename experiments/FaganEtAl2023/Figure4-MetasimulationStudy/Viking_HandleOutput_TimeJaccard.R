# Script Goals: ###############################################################
#   0. In Parallel,
#   1. Load the Data,
#   2. Extract temporal Jaccard for each environment over time,
#   3. Save the beta values separately, but alongside their parameters.

# Parameters: #################################################################
divide_time_by <- 1E4
#preferred_rows_per_event <- 8

outputLocation <- file.path(".", paste0("SaveTimeJaccard_", Sys.Date()))
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
    "RMTRCode2_2022-02-12_1.0.tar.gz", lib = librarypath,
    repos = NULL, type = "source"
  )
}
library(RMTRCode2)#, lib.loc = librarypath) # lib.loc shouldn't be necessary.

if (!exists("Calculate_TimeJaccard",
            where = asNamespace("RMTRCode2"),
            mode = "function")) {
  source(file.path(".", "Calculate_TimeJaccard.R"))
  if (!exists("Calculate_TimeJaccard")) {
    stop("Could not find Calculate_TimeJaccard.")
  }
}

# Receive the Selection: ######################################################
print("Receiving number of cores:")

cargs <- 2
# cargs <- as.numeric(commandArgs(trailingOnly = TRUE))

clust <- parallel::makeCluster(cargs[1], outfile = "")
doParallel::registerDoParallel(clust)

# Send locally installed packages to cluster.
parallel::clusterExport(
  cl = clust, list("allLibraryPaths")
)

parallel::clusterEvalQ(
  cl = clust,
  expr = {
    .libPaths(allLibraryPaths)
    library(RMTRCode2)
    library(vegan)
  }
)

# Files: ######################################################################
print("Identifying files.")

directories <- dir(
  path = ".", pattern = "SaveOutput[_]",
  #pattern = "Viking[_]SaveOutput[_]2021[-]12[-]28[_]2022[-]03[-]01",
  full.names = TRUE, include.dirs = TRUE)
files <- dir(path = directories,
             pattern = "^MNA[-]Master.+[.]RData$",
             full.names = TRUE, recursive = TRUE,
             include.dirs = TRUE, no.. = TRUE)

# Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
filesParameters <- dir(path = "Prepared_2021-12-20",
                       pattern = "^MNA[-]Master.+Cases[.]csv$",
                       full.names = TRUE, recursive = TRUE,
                       include.dirs = TRUE, no.. = TRUE)
parametersMaster <- load(file.path("Prepared_2021-12-20",
                                   "MNA-MasterParameters.RData"))
parametersCSVs <- lapply(filesParameters, utils::read.csv)

# Note the problem observed 31/01/2022, 1400 in the calendar.
# To accomodate, reassign inf to the *end* of the Space list.
systemMods$SpaceDistanceMultiplier <- systemMods$SpaceDistanceMultiplier[c(2:9, 1)]

# 0 In Parallel: ##############################################################
results <- foreach::foreach(
  file = iterators::iter(files),
  .packages = "dplyr"
) %dopar% {
  gc() # Try to reduce the memory usage.
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
    paste0("MNA-Master", LETTERS[idNums[1]], "-Cases-TimeJaccard-",
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

    # 2 Extract Beta diversities over time: ######################
    if (length(idNums) == 2) { # fileContents is a list of results
      TimeJaccard <- list(
        Betas = lapply(
          fileContents,
          FUN = Calculate_TimeJaccard,
          nspecies = systemBase$species
        ))
    } else if (length(idNums) == 4) {# fileContents is a result
      TimeJaccard <- list(
        Betas = list(Calculate_TimeJaccard(
          fileContents,
          nspecies = systemBase$species,
          minTime = max(fileContents$Events$Times)/101
        ))
      )
    } else stop("idNums is not 2 or 4.")

    # Formatting is now the same, so we can apply a single function to
    # make sure the time units are the same as in other functions.

    TimeJaccard$Betas <- lapply(TimeJaccard$Betas, function(b) {
      b %>% dplyr::mutate(
        Time  = Time / divide_time_by
      )
    })

    print("diversity")
    # 6. Save the biodiversity values...: #####################################
    # separately, but alongside its parameters.
    # idNums indicates which variation we are dealing with, as well as the
    # associated parameters.
    # This is split between MNA-MasterParameters.RData and the .csv files.
    theParameters <- parametersCSVs[[idNums[1]]][((idNums[2] - 1) %/% 10) + 1,]
    TimeJaccard$PoolMod <- systemMods$PoolBodySizes[theParameters$Framework, ]
    TimeJaccard$NoiseMod <- systemMods$InteractionsNoiseMultiplier[theParameters$Framework]
    TimeJaccard$NeutralMod<- systemMods$NeutralRateMultipliers[theParameters$Neutral,]
    TimeJaccard$SpaceMod <- systemMods$SpaceDistanceMultiplier[theParameters$Space]

    save(TimeJaccard,
         file = fileName
    )
    print("save")

    # Ran into memory issue. We'll see if this fixes it.
    rm(TimeJaccard)
    rm(fileContents)

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
