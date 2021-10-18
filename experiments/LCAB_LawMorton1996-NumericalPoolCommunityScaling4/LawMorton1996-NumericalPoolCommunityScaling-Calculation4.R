# Load all libraries. #########################################################
# Location:
directory <- '.'
versionNum <- 4

print("Loading libraries")
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

packages <- c(
  "foreach",
  "parallel",
  "doParallel",
  "dplyr",
  "tidyr",
  "Rcpp"
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
    "RMTRCode2_1.0.tar.gz", lib = librarypath,
    repos = NULL, type = "source"
  )
}
library(RMTRCode2)#, lib.loc = librarypath) # lib.loc shouldn't be necessary.

# Create folder for saving.
directory_save <- file.path(directory, paste0("save-", Sys.Date()))
# if (dir.exists(directory_save)) {
#   directory_save <- paste0(directory_save, "-", length(dir(directory)) + 1)
# }
if (!dir.exists(directory_save)) {
  dir.create(directory_save, recursive = TRUE)
}

# Run sanity checks. ##########################################################

print(paste("Cores detected: ", parallel::detectCores()))

cargs <- as.numeric(commandArgs(TRUE))

print('Args provided:')
print(cargs)

# Below does not work: mismatch between assigned cores and cores detectable.
# cores <- max(c(parallel::detectCores(), 1))
cores <- cargs[1]
clust <- parallel::makeCluster(cores, outfile = "")
doParallel::registerDoParallel(clust)

if (!is.na(cargs[2])) {
  start <- cargs[2]
} else {
  start <- 1
}

# Perform calculation. ########################################################
source(paste0(
    "LawMorton1996-NumericalPoolCommunityScaling-Settings", versionNum, ".R"
  ))

load(
  file = paste0(
    "LawMorton1996-NumericalPoolCommunityScaling-PoolMats", versionNum, ".RDS"
    )
)

end <- length(pools)

results <- foreach::foreach(
  i = start : end,
  ourPool = iterators::iter(pools[start:end]),
  ourMat = iterators::iter(communityMats[start:end])#,
  #.packages = c("RMTRCode2", packages),
  #.export = c("allLibraryPaths")
) %:% foreach::foreach(
  seed = iterators::iter(seedsRun[((i - 1) * runs + 1) : (i * runs)]),
  .combine = "rbind"
) %dopar% {
  .libPaths(allLibraryPaths)
  library(RMTRCode2)
  
  savepath <- file.path(
    directory_save,
    paste0("run-", i, "-", seed, ".RDS")
  )

  print(paste(i, seed, Sys.time(), "checking", savepath))
  # Don't bother if we already have done this calculation.
  # Don't load it either yet, we likely have a large number of other calcs to do
  if (file.exists(savepath)) {
    print(paste(i, seed, Sys.time(), "killing"))
    return(
      data.frame(
        Combination = i,
        Result = savepath
      )
    )
  }
  
  print(paste(i, seed, Sys.time(), "running"))

  aRun <- LawMorton1996_NumericalAssembly(
    Pool = ourPool,
    CommunityMat = ourMat,
    IntegratorTimeStep = 10000,
    InnerTimeStepSize = 5000,
    ArrivalEvents = events,
    ReturnValues = c("Abundance", "Sequence"),
    seed = seed
  )
  
  print(paste(i, seed, Sys.time(), "summarizing"))

  lastRow <- which(is.na(aRun$Sequence$Community[-1]))[1]
  if (is.na(lastRow)) {
    lastRow <- nrow(aRun$Sequence)
  }
  aRun$Sequence <- aRun$Sequence[1:lastRow, ]

  retval <- if (LawMorton1996_CheckUninvadable(
    AbundanceRow = aRun$Abundance[nrow(aRun$Abundance), ],
    Pool = ourPool,
    CommunityMatrix = ourMat
  )) {
    # and, if so, store the result.
    #return(
      list(
        Combination = i,
        Result = aRun$Sequence,
        Abund = aRun$Abundance[nrow(aRun$Abundance), ]
      )
    #)
  } else {
    #return(
      list(
        Combination = i,
        Result = toString(NULL),
        Abund = toString(NULL)
      )
    #)
  }
  
  print(paste(i, seed, Sys.time(), "saving"))

  save(
    file = savepath,
    retval
  )
  
  print(paste(i, seed, Sys.time(), "returning"))

  return(retval)
}

# save(
 #  file = "LawMorton1996-NumericalPoolCommunityScaling-results.RDS",
  # results
# )

parallel::stopCluster(clust)

