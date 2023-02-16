# See Note Wednesday 02/03/2022, 1930 for explanation.
# Roughly, attributes were misassigned.
# This file loads, reassigns attributes, and saves the
# file UNDER THE SAME NAME IT WAS LOADED (overwriting).

# Note the script is not smart.
# It can't detect when the attributes are wrong.
# We exploit the fact that we'll be making a list
# and then can start generating files with the patched
# diversity calculator without changing the list
# to continue working in parallel.

# For job script:
# The largest file is 15 MB, requiring 1.3 GB to open.
# This suggests that 4GB is probably enough to perform
# our operations. (Per Core/Task of course.)

# Strategy:
# Create a list of the files to rectify.
# Assign each task to a core (somehow).
# Each task reads its file name,
# Reads its file and correct parameters,
# Corrects the file, and
# Saves and overwrites the original file.

# Two philosophies:
# We could create a manager, which drafts the list
# and assigns workers using foreach.
# Alternatively, we could create a script that drafts
# the list and then submits a large set of workers
# to process the list.

# It is probably best for recording purposes to save
# the list as a file either way.

# Parameters: #################################################################

directory <- '.'

fileLocation <- file.path(
  directory, "SaveDiversity_2022-03-05"
  )

files <- dir(fileLocation,
             pattern = "[.]RData$",
             full.names = TRUE)

utils::write.csv(
  x = data.frame(files = files, success = NA),
  file = paste0("RectifiedFiles_Attempt_", Sys.Date(), ".csv")
)

# Libraries: ##################################################################

print("Loading libraries")
librarypath <- file.path(".", "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

packages <- c(
  "parallel",   # Base parallel dependency
  "doParallel", # For compatible cluster
  "foreach",    # For parallel evaluation
  "iterators"   # For parallel evaluation
)

for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, lib = librarypath,
                     repos = 'https://cloud.r-project.org',
                     dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

# Receive the Selection: ######################################################
print("Receiving number of cores:")

cargs <- 2
# cargs <- as.numeric(commandArgs(trailingOnly = TRUE))

clust <- parallel::makeCluster(cargs[1], outfile = "")
doParallel::registerDoParallel(clust)

# Parameter Files: ############################################################
# Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
filesParameters <- dir(
  path = file.path(
    directory,
    "Prepared_2021-12-20"
  ),
  pattern = "^MNA[-].+Cases[.]csv$",
  full.names = TRUE, recursive = TRUE,
  include.dirs = TRUE, no.. = TRUE
)
parametersMaster <- load(file.path(directory, "Prepared_2021-12-20", "MNA-MasterParameters.RData"))
parametersCSVs <- lapply(filesParameters, utils::read.csv)

# Note the problem observed 31/01/2022, 1400 in the calendar.
# To accomodate, reassign inf to the *end* of the Space list.
systemMods$SpaceDistanceMultiplier <- systemMods$SpaceDistanceMultiplier[c(2:9, 1)]

parallel::clusterExport(clust, parametersMaster)
parallel::clusterExport(clust, "parametersCSVs")

# In Parallel: ################################################################
results <- foreach::foreach(
  file = iterators::iter(files), .packages = "dplyr"
) %dopar% {
  idNums <- suppressWarnings(na.omit(
    as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                             split = "-", fixed = TRUE)[[1]],
                    n = 4))),
    classes = "simpleWarning")

  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    if (!exists("Diversity")) {
      stop(paste("Diversity not found in file", file))
    }

    print("loaded")

    theParameters <- parametersCSVs[[idNums[1]]][((idNums[2] - 1) %/% 10) + 1,]
    Diversity$PoolMod <- systemMods$PoolBodySizes[theParameters$Framework, ]
    Diversity$NoiseMod <- systemMods$InteractionsNoiseMultiplier[theParameters$Framework]
    Diversity$NeutralMod<- systemMods$NeutralRateMultipliers[theParameters$Neutral,]
    Diversity$SpaceMod <- systemMods$SpaceDistanceMultiplier[theParameters$Space]

    save(Diversity,
         file = file
    )
    print("save")

    return(data.frame(
      file = file,
      success = TRUE
    ))
  }, error = function(e, id = file) {
    return(data.frame(
      file = id,
      success = e
    ))
  })
}

# Cleanup: ####################################################################
parallel::stopCluster(clust)

# Save for checking: ##########################################################
utils::write.csv(
  x = do.call(rbind, results),
  file = paste0("RectifiedFiles_", Sys.Date(), ".csv")
)
