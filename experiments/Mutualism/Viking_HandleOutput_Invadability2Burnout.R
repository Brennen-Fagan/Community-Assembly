# Script Goals: ###############################################################
#   0. In Parallel,
#   1. Load the Data,
#   2. Extract invadabilities (by occupancy (site-species pair)) over time,
#   3. Save the invadability values separately, but alongside their parameters.

# Parameters: #################################################################
SpeciesTypes <- c(34, 66) * 2
divide_time_by <- 1E2
preferred_rows_per_event <- 1.5

outputLocation <- file.path(".", paste0("SaveInvadability_", Sys.Date()))
if (!dir.exists(outputLocation)) {
  dir.create(outputLocation, showWarnings = FALSE)
}

# Libraries: ##################################################################
# print("Loading libraries")
# librarypath <- file.path(".", "Rlibs")
# if (!dir.exists(librarypath)) {
#   dir.create(librarypath, showWarnings = FALSE)
# }
# .libPaths(c(librarypath, .libPaths()))
#
# allLibraryPaths <- .libPaths()

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
  # if (!require(package, character.only = TRUE)) {
  #   install.packages(package, lib = librarypath,
  #                    repos = 'https://cloud.r-project.org',
  #                    dependencies = TRUE)
  # }
  library(package, character.only = TRUE)
}

# if (!require("RMTRCode2", character.only = TRUE)) {
#   install.packages(
#     "RMTRCode2_2022-02-12_1.0.tar.gz", lib = librarypath,
#     repos = NULL, type = "source"
#   )
# }
library(RMTRCode2)#, lib.loc = librarypath) # lib.loc shouldn't be necessary.

# if (!exists("CalculateLocalInvadables_KnockOn",
#             where = asNamespace("RMTRCode2"),
#             mode = "function")) {
#   source(file.path(".", "Invadability.R"))
#   if (!exists("CalculateLocalInvadables_KnockOn")) {
#     stop("Could not find invadability calculation.")
#   }
# }

# # Receive the Selection: ######################################################
# print("Receiving number of cores:")
#
cargs <- 4
# cargs <- as.numeric(commandArgs(trailingOnly = TRUE))
#
clust <- parallel::makeCluster(cargs[1], outfile = "")
doParallel::registerDoParallel(clust)
#
# # Send locally installed packages to cluster.
# parallel::clusterExport(
#   cl = clust, list("allLibraryPaths")
# )
#
parallel::clusterEvalQ(
  cl = clust,
  expr = {
    #     .libPaths(allLibraryPaths)
    library(RMTRCode2)
    library(Matrix)
  }
)

# Files: ######################################################################
print("Identifying files.")

# directories <- dir(
#   path = ".", pattern = "SaveOutput[_]MissingInv",
#   #pattern = "Viking[_]SaveOutput[_]2021[-]12[-]28[_]2022[-]03[-]01",
#   full.names = TRUE, include.dirs = TRUE)
directories <- c("Data_2023-06-26", "Data_2023-06-27", "Data_2023-07-06",
                 "Data2_2023-06-30")
files <- dir(
  path = directories,
  pattern = "^.+[-]ExampleExtProp[-].+[.]RData$",
  full.names = TRUE, recursive = TRUE,
  include.dirs = TRUE, no.. = TRUE)

# # Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
# filesParameters <- dir(path = "Prepared_2021-12-20",
#                        pattern = "^MNA[-]Master.+Cases[.]csv$",
#                        full.names = TRUE, recursive = TRUE,
#                        include.dirs = TRUE, no.. = TRUE)
# parametersMaster <- load(file.path("Prepared_2021-12-20",
#                                    "MNA-MasterParameters.RData"))
# parametersCSVs <- lapply(filesParameters, utils::read.csv)

# We'll need the per capita dynamics and the dispersal matrices.
# filesPrepared <- dir(path = directories,
#                      pattern = "^.+[-]ExampleOutcome[-]PoolMats[-].+[.]RData$",
#                      full.names = TRUE, recursive = TRUE,
#                      include.dirs = TRUE, no.. = TRUE)
# preparedCases <- lapply(filesPrepared, function(fil) {
#   tempenv <- new.env()
#   load(fil, envir = tempenv)
#   # We only need the matrices and the reproductive rates, not the events.
#   # This should cut down on the memory requirements.
#   rm(list = c("evnts"), envir = tempenv)
#   return(tempenv)
# })


# 0 In Parallel: ##############################################################
# start_time <- Sys.time()
results <- foreach::foreach(
  file = iterators::iter(files),
  .packages = c("dplyr")#, "RMTRCode2")
) %dopar% {
  gc() # Try to reduce the memory usage.
  print(file)
  # # Retrieve the trailing id numbers from before the file extension.
  # # 1st: Set, 2nd: CaseNumber, 3rd: History, 4th: Part
  # # Note, if last two are not present, all histories are bundled together.
  # # If last two are present, each file has a single part of a single history.
  # idNums <- suppressWarnings(na.omit(
  #   as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
  #                            split = "-", fixed = TRUE)[[1]],
  #                   n = 4))),
  #   classes = "simpleWarning")
  #
  #   print(idNums)

  ids <- strsplit(tools::file_path_sans_ext(basename(file)),
                  split = "-", fixed = TRUE)[[1]][c(1, 6)]
  print(ids)


  fileName <- file.path(
    outputLocation,
    paste0(ids[1], "-Invadability-", ids[2], ".RData")
  )

  if (file.exists(fileName)) {
    print("Already done.")
    return(ids)
  }

  # Build up Required Objects: ################################################
  # See, e.g., VikingRunHalves.R

  load(file.path(dirname(file),
                 paste0(ids[1], "-ExampleOutcome-PoolMats-Env10.RData")))

  IntMat <- Matrix::bdiag(InteractionMatrices$Mats)


  # IntMats <- Matrix::bdiag(
  #   preparedCases[[idNums[1]]]$matrs[[ # Case
  #     ((idNums[2] - 1) %/% 10) + 1 # Parameters
  #     ]][[((idNums[2] - 1) %% 10) + 1 # System
  #         ]]$Mats
  # )

  # pool <- preparedCases[[idNums[1]]]$pools[[ # Case
  #   ((idNums[2] - 1) %/% 10) + 1 # Parameters
  #   ]][[((idNums[2] - 1) %% 10) + 1 # System
  #       ]]

  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    fileContents <- get(fileContents)

    print("loaded")

    numEnvs <- fileContents$NumEnvironments

    # Create objects dependent on loaded contents.
    if (ids[1] %in% c(
      "MNA", "LM1996PermuteWithin", "LM1996PermuteWithinPool")) {
      PerCapitaDynamics <- RMTRCode2::PerCapitaDynamics_Type1(
        Pool$ReproductionRate, IntMat,
        NumEnvironments = numEnvs
      )
    } else if (ids[1] %in% c("Mutualism")) {
      PerCapitaDynamics <- PerCapitaDynamics_Mutualistic1(
        reprate, IntMat,
        NumEnvironments = numEnvs,
        SpeciesTypes = SpeciesTypes, Saturations = Saturation
      )
    } else if (ids[1] %in% c("MutTF2010", "PrPrTF2010")) {
      PerCapitaDynamics <- PerCapitaDynamics_Mutualistic3(
        reprate, IntMat,
        NumEnvironments = numEnvs,
        SpeciesTypes = SpeciesTypes, Saturations = Saturation
      )
    } else {
      stop("Ids[1] not recognised; typo in conditions?")
    }

    # Move this earlier than in _Diversity.R.
    # We need it to know the amount of dispersal.
    # theParameters <- parametersCSVs[[idNums[1]]][((idNums[2] - 1) %/% 10) + 1,]

    # PerIslandDistance <- systemMods$SpaceDistanceMultiplier[theParameters$Space]
    PerIslandDistance <- 10^as.numeric(sub("_", "-", ids[2], fixed = T))

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
      SpeciesSpeeds = rep(1, nrow(Pool))
    )

    # 2 Extract Invadabilities over time: ######################
    Invadability <- list(
      Invadabilities = list(thinAndCalculateInvadabilities(
        fileContents,
        dyn = PerCapitaDynamics,
        dis = DispersalMatrices,
        preferred_rows_per_event = preferred_rows_per_event,
        divide_time_by = divide_time_by,
        burnout =
          with(fileContents$Events, Times[length(Times)]) +
          fileContents$ReactionTime
      ))
    )
    print("Invadability")
    # 3. Save the invadability values...: #####################################
    # separately, but alongside its parameters.
    # idNums indicates which variation we are dealing with, as well as the
    # associated parameters.
    # This is split between MNA-MasterParameters.RData and the .csv files.
    Diversity$SysType <- ids[1]
    Diversity$SpaceMod <- ids[2]

    save(Invadability,
         file = fileName
    )
    print("save")

    ###   3. Record Data as Present: ##########################################
    # This is already in the file list, but we can make it easier by taking the
    # id numbers. If there is an error at any point, record the relevant id
    # numbers and the error to know that something went wrong.
    return(ids)
  }, error = function(e, id = ids) {
    print(paste(toString(id), e))
    return(list(ID = id,
                error = e))
  })
}
# end_time <- Sys.time()

# Cleanup: ####################################################################
parallel::stopCluster(clust)
