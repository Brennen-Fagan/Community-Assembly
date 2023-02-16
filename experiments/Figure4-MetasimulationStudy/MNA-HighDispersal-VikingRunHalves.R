# For running on Viking in an embarrassingly parallel manner.
# We have the cases, pools, matr(ice)s, and ev(e)nts already generated.
# We now need
#   to receive the selection to use,
#   to load the appropriate system,
#   to generate the master interaction matrix,
#   to generate the per capita dynamics,
#   to generate the spatial dynamics,
#   to run the system, and
#   finally to save the system.

# Furthermore, there is a danger of a segfault during evaluation.
# To avoid this, we'll be doing two things.
# First, run the code in halves. We'll run 1-1, 1-2, 2-1, 2-2, 3-1, 3-2, etc.
# Second, we'll be saving after each half in order to be on the safe side.
# We can then skip halves that have already been computed.

# Universal Settings: #########################################################
EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 30 # Number of steps to reach next event to smooth.

directory <- '.'

outputLocation <- file.path(directory, paste0("SaveOutput_", Sys.Date()))
if (!dir.exists(outputLocation)) {
  dir.create(outputLocation, showWarnings = FALSE)
}

# Set-up: #####################################################################
print("Loading libraries")
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

packages <- c(
  "dplyr",
  "Matrix"
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

# Receive the Selection: ######################################################
print("Receive")
# Two CArgs:
#   Case:       A: 1,
#   Row Number: A: 1000, (i.e. `nrow`s)
cargs <- as.numeric(commandArgs(trailingOnly = TRUE))
print(cargs)

# Load the System: ############################################################
print("Load")
cases <- NULL
pools <- NULL
matrs <- NULL
evnts <- NULL

candidates <- dir(pattern = "MNA[-]HiDisp.+[-]Cases[-]Prepared[.]RData$",
                  path = "Prepared_2022-05-03",
                  full.names = TRUE)

stopifnot(length(candidates) == 1)

loaded <- load(candidates[cargs[1]])
outputName <- strsplit(basename(candidates[cargs[1]]),
                       split = ".", fixed = TRUE)[[1]]
outputName <- paste0(outputName[-length(outputName)],
                     "-", cargs[1], "-", cargs[2], "-")#, ".RData")

stopifnot(
  c("cases", "pools", "matrs", "evnts") %in% loaded,
  unlist(lapply(list(cases, pools, matrs, evnts),
                function(i) !is.null(i)))
)

case <- cases[cargs[2], ]
pool <- pools[[case$Parameters]][[case$System]]
matr <- matrs[[case$Parameters]][[case$System]]
evnt <- evnts[[case$Parameters]][[case$System]]

Environments <- length(matr$Mats)

# Generate the master interaction matrix: #####################################
print("Master")
IntMats <- Matrix::bdiag(matr$Mats)

# Generate the per capita dynamics: ###########################################
print("PerCap")
PerCapitaDynamics <- RMTRCode2::PerCapitaDynamics_Type1(
  pool$ReproductionRate, IntMats,
  NumEnvironments = Environments
)

# Generate the spatial dynamics: ##############################################
print("Space")
PerIslandDistance <- switch(
  case$Space,
  Inf, 1e9, 1e8, 1e7, 1e6,
  1e5, 1e4, 1e3, 1e2, 1e1, 1e0
)

# Ring Dynamics
DistanceMatrix <- Matrix::bandSparse(
  Environments, k = c(-1, 1),
  diagonals = list(rep(PerIslandDistance, Environments - 1),
                   rep(PerIslandDistance, Environments - 1))
)
DistanceMatrix[Environments, 1] <- PerIslandDistance
DistanceMatrix[1, Environments] <- PerIslandDistance

DispersalMatrices <- RMTRCode2::CreateDispersalMatrix(
  EnvironmentDistances = DistanceMatrix,
  SpeciesSpeeds = rep(1, nrow(pool))
)

# Run the system: #############################################################
print("Run")
# Run times appear to be short enough that we can just run each set of histories
# (per pool-environ combination) on the same node without hitting time-out.
#results <- list()
folderCandidates <- dir(pattern = "SaveOutput[_][0-9]{4}[-][0-9]{2}[-][0-9]{2}")
fileCandidates_full <- dir(folderCandidates, full.names = TRUE)
fileCandidates <- basename(fileCandidates_full)
for (i in seq_along(evnt)){
  print(i)
  fileName1 <- file.path(
    outputLocation, paste0(
      outputName, i, "-", 1, ".RData"
    )
  )

  ev <- evnt[[i]]
  stopTime <- max(ev$Events$Times) / 2

  if (
    #!file.exists(fileName1)
    !any(indices <- (paste0(
      outputName, i, "-", 1, ".RData"
    ) == fileCandidates))
    ) {
    print(1)
    ev$Events <- ev$Events %>% dplyr::filter(Times < stopTime)
    results <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
      Pool = pool, NumEnvironments = Environments,
      InteractionMatrices = matr,
      Events = ev,
      PerCapitaDynamics = PerCapitaDynamics,
      DispersalMatrix = DispersalMatrices,
      EliminationThreshold = EliminationThreshold,
      ArrivalDensity = ArrivalDensity,
      MaximumTimeStep = MaximumTimeStep,
      BetweenEventSteps = BetweenEventSteps,
      Verbose = FALSE
    )

    print("Save")

    save(
      results, file = fileName1
    )
  }

  fileName2 <- file.path(
    outputLocation, paste0(
      outputName, i, "-", 2, ".RData"
    )
  )

  if (
    #!file.exists(fileName2)
    !any((paste0(
      outputName, i, "-", 2, ".RData"
    ) == fileCandidates))
    ) {

    if(!exists("results")) {
      load(tail(fileCandidates_full[indices], n = 1))
    }

    print(2)

    ev <- evnt[[i]]
    # We'll stop immediately before the last event, then make sure to include
    # it as the starting point for the next event.
    lastEventTime <- max(results$Events$Times)
    ev$Events <- ev$Events %>% dplyr::filter(Times >= lastEventTime)
    stopRow <- which.max(results$Abundance[, 1] >= lastEventTime) - 1

    results <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
      PopulationInitial = results$Abundance[stopRow, -1],
      TimeInitial = results$Abundance[stopRow, 1],
      Pool = pool, NumEnvironments = Environments,
      InteractionMatrices = matr,
      Events = ev,
      PerCapitaDynamics = PerCapitaDynamics,
      DispersalMatrix = DispersalMatrices,
      EliminationThreshold = EliminationThreshold,
      ArrivalDensity = ArrivalDensity,
      MaximumTimeStep = MaximumTimeStep,
      BetweenEventSteps = BetweenEventSteps,
      Verbose = FALSE
    )

    print("Save")

    save(
      results, file = fileName2
    )
  }

  rm(results)
}


# Save the system: ############################################################
# print("Save")
# save(
#   results, file = file.path(
#     outputLocation, outputName
#   )
# )
