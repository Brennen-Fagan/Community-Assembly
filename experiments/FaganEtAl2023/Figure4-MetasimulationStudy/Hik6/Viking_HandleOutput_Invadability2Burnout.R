# Script Goals: ###############################################################
#   0. In Parallel,
#   1. Load the Data,
#   2. Extract invadabilities (by occupancy (site-species pair)) over time,
#   3. Save the invadability values separately, but alongside their parameters.

# Parameters: #################################################################
divide_time_by <- 1E4
preferred_rows_per_event <- 1.5

burnout <- 6.5 * 1E4 # In original units.

outputLocation <- file.path(".", paste0("SaveInvadability_", Sys.Date()))
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
  "parallel",   # Base parallel dependency
  "doParallel", # For compatible cluster
  "foreach",    # For parallel evaluation
  "iterators"#,   # For parallel evaluation
  #"WeightedCluster" # For Medoids to reduce computational load.
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

if (!exists("CalculateLocalInvadables_KnockOn",
            where = asNamespace("RMTRCode2"),
            mode = "function")) {
  if (!exists("CalculateLocalInvadables_KnockOn")) {
    stop("Could not find invadability calculation.")
  }
}

# # Receive the Selection: ######################################################
print("Receiving number of cores:")

# cargs <- 4
cargs <- as.numeric(commandArgs(trailingOnly = TRUE))

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
    library(Matrix)
  }
)

# Files: ######################################################################
print("Identifying files.")

directories <- dir(
  path = ".", pattern = "SaveOutput",
  #pattern = "Viking[_]SaveOutput[_]2021[-]12[-]28[_]2022[-]03[-]01",
  full.names = TRUE, include.dirs = TRUE)
files <- dir(path = directories,
             pattern = "^MNA[-]Hik6.+[.]RData$",
             full.names = TRUE, recursive = TRUE,
             include.dirs = TRUE, no.. = TRUE)

# Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
filesParameters <- dir(path = "Prepared_2023-03-03",
                       pattern = "^MNA[-]Hik6.+Cases[.]csv$",
                       full.names = TRUE, recursive = TRUE,
                       include.dirs = TRUE, no.. = TRUE)
parametersMaster <- load(file.path("Prepared_2023-03-03",
                                   "MNA-Hik6Parameters.RData"))
parametersCSVs <- lapply(filesParameters, utils::read.csv)

# We'll need the per capita dynamics and the dispersal matrices.
filesPrepared <- dir(path = "Prepared_2023-03-03",
                     pattern = "^MNA[-]Hik6.+[-]Cases[-]Prepared[.]RData$",
                     full.names = TRUE, recursive = TRUE,
                     include.dirs = TRUE, no.. = TRUE)
preparedCases <- lapply(filesPrepared, function(fil) {
  tempenv <- new.env()
  load(fil, envir = tempenv)
  # We only need the matrices and the reproductive rates, not the events.
  # This should cut down on the memory requirements.
  rm(list = c("evnts"), envir = tempenv)
  return(tempenv)
})


# 0 In Parallel: ##############################################################
# start_time <- Sys.time()
results <- foreach::foreach(
  file = iterators::iter(files),
  .packages = c("dplyr")
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
    paste0("MNA-Hik6", LETTERS[idNums[1]], "-Cases-Invadability-",
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

  if (length(idNums) == 4 && idNums[4] == 1) {
    print("Not end of simulation.")
    return(idNums)
  }

  # Build up Required Objects: ################################################
  # See, e.g., VikingRunHalves.R

  IntMats <- Matrix::bdiag(
    preparedCases[[idNums[1]]]$matrs[[ # Case
      ((idNums[2] - 1) %/% 10) + 1 # Parameters
      ]][[((idNums[2] - 1) %% 10) + 1 # System
          ]]$Mats
  )

  pool <- preparedCases[[idNums[1]]]$pools[[ # Case
    ((idNums[2] - 1) %/% 10) + 1 # Parameters
    ]][[((idNums[2] - 1) %% 10) + 1 # System
        ]]

  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    fileContents <- get(fileContents)

    print("loaded")

    if (length(idNums) == 2) {
      numEnvs <- fileContents[[1]]$NumEnvironments
    } else if (length(idNums) == 4) {
      numEnvs <- fileContents$NumEnvironments
    }

    # Create objects dependent on loaded contents.
    PerCapitaDynamics <- RMTRCode2::PerCapitaDynamics_Type1(
      pool$ReproductionRate, IntMats,
      NumEnvironments = numEnvs
    )

    # Move this earlier than in _Diversity.R.
    # We need it to know the amount of dispersal.
    theParameters <- parametersCSVs[[idNums[1]]][((idNums[2] - 1) %/% 10) + 1,]

    PerIslandDistance <- systemMods$SpaceDistanceMultiplier[theParameters$Space]

    # Ring Dynamics
    DistanceMatrix <- Matrix::bandSparse(
      numEnvs, k = c(-1, 1),
      diagonals = list(rep(PerIslandDistance, numEnvs - 1),
                       rep(PerIslandDistance, numEnvs - 1))
    )
    DistanceMatrix[numEnvs, 1] <- PerIslandDistance
    DistanceMatrix[1, numEnvs] <- PerIslandDistance

    DispersalMatrices <- RMTRCode2::CreateDispersalMatrix(
      EnvironmentDistances = DistanceMatrix,
      SpeciesSpeeds = rep(1, nrow(pool))
    )

    # 2 Extract Invadabilities over time: ######################
    if (length(idNums) == 2) { # fileContents is a list of results
      Invadability <- list(
        Invadabilities = lapply(
          fileContents,
          FUN = thinAndCalculateInvadabilities,
          dyn = PerCapitaDynamics,
          dis = DispersalMatrices,
          preferred_rows_per_event = preferred_rows_per_event
        ))
    } else if (length(idNums) == 4) {# fileContents is a result
      Invadability <- list(
        Invadabilities = list(thinAndCalculateInvadabilities(
          fileContents,
          dyn = PerCapitaDynamics,
          dis = DispersalMatrices,
          preferred_rows_per_event = preferred_rows_per_event
        ))
      )
    } else stop("idNums is not 2 or 4.")
    print("Invadability")
    # 3. Save the invadability values...: #####################################
    # separately, but alongside its parameters.
    # idNums indicates which variation we are dealing with, as well as the
    # associated parameters.
    # This is split between MNA-MasterParameters.RData and the .csv files.
    Invadability$PoolMod <- systemMods$PoolBodySizes[theParameters$Framework, ]
    Invadability$NoiseMod <- systemMods$InteractionsNoiseMultiplier[theParameters$Framework]
    Invadability$NeutralMod<- systemMods$NeutralRateMultipliers[theParameters$Neutral,]
    Invadability$SpaceMod <- PerIslandDistance

    save(Invadability,
         file = fileName
    )
    print("save")

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
    print(paste(toString(id), e))
    return(list(ID = id,
                error = e))
  })
}
# end_time <- Sys.time()

# Cleanup: ####################################################################
parallel::stopCluster(clust)
