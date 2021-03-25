# Load all libraries. #########################################################
# Location:
directory <- '.'

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

# Run sanity checks. ##########################################################

print(paste("Cores detected: ", parallel::detectCores()))

cargs <- as.numeric(commandArgs(TRUE))

print('Args provided:')
print(cargs)

# Below does not work: mismatch between assigned cores and cores detectable.
# cores <- max(c(parallel::detectCores(), 1))
cores <- cargs[1]
clust <- parallel::makeCluster(cores)
doParallel::registerDoParallel(clust)

# Perform calculation. ########################################################
set.seed(38427042)

basal <- c(3, 10, 30, 100, 300, 1000)
consumer <- c(3, 10, 30, 100, 300, 1000) * 2
events <- (max(basal) + max(consumer)) * 2
runs <- 100

logBodySize <- c(-2, -1, -1, 1) # Morton and Law 1997 version.
parameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)

# Need to rerun seedsPrep to get the random number generation right for seedsRun
seedsPrep <- runif(2 * length(basal) * length(consumer)) * 1E8
seedsRun <- runif(runs * length(basal) * length(consumer)) * 1E8

load(
  file = "LawMorton1996-NumericalPoolCommunityScaling-PoolMats.RDS"
)

results <- foreach::foreach(
  i = 1 : (length(basal) * length(consumer)),
  ourPool = iterators::iter(pools),
  ourMat = iterators::iter(communityMats)#,
  #.packages = c("RMTRCode2", packages),
  #.export = c("allLibraryPaths")
) %:% foreach::foreach(
  seed = iterators::iter(seedsRun[((i - 1) * runs + 1) : (i * runs)])
) %dopar% {
  .libPaths(allLibraryPaths)
  library(RMTRCode2)

  aRun <- LawMorton1996_NumericalAssembly(
    Pool = ourPool,
    CommunityMat = ourMat,
    IntegratorTimeStep = 10000,
    InnerTimeStepSize = 5000,
    ArrivalEvents = events,
    ReturnValues = c("Abundance", "Sequence"),
    seed = seed
  )

  lastRow <- which(is.na(aRun$Sequence$Community[-1]))[1]
  if (is.na(lastRow)) {
    lastRow <- nrow(aRun$Sequence)
  }
  aRun$Sequence <- aRun$Sequence[1:lastRow, ]

  if (LawMorton1996_CheckUninvadable(
    AbundanceRow = aRun$Abundance[nrow(aRun$Abundance), ],
    Pool = ourPool,
    CommunityMatrix = ourMat
  )) {
    # and, if so, store the result.
    return(aRun$Sequence)
  } else {
    return(toString(NULL))
  }
}

save(
  file = "LawMorton1996-NumericalPoolCommunityScaling-results.RDS",
  results
)

parallel::stopCluster(clust)

