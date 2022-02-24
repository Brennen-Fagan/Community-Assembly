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
    "RMTRCode2_1.0.tar.gz", lib = librarypath,
    repos = NULL, type = "source"
  )
}
library(RMTRCode2)#, lib.loc = librarypath) # lib.loc shouldn't be necessary.

# Receive the Selection: ######################################################
print("Receive")
# Two CArgs:
#   Case:       A: 1,   B: 2,   C: 3
#   Row Number: A: 400, B: 160, C: 1600 (i.e. `nrow`s)
cargs <- c(1, 350) #as.numeric(commandArgs(trailingOnly = TRUE))

# Load the System: ############################################################
print("Load")
cases <- NULL
pools <- NULL
matrs <- NULL
evnts <- NULL

candidates <- dir(pattern = "MNA[-]Master.+[-]Cases[-]Prepared[.]RData$",
                  path = "Prepared_2021-12-20",
                  full.names = TRUE)

stopifnot(length(candidates) == 3)

loaded <- load(candidates[cargs[1]])
outputName <- strsplit(basename(candidates[cargs[1]]),
                       split = ".", fixed = TRUE)[[1]]
outputName <- paste0(outputName[-length(outputName)],
                     "-", cargs[1], "-", cargs[2], "-stop", ".RData")

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
  1e9, 1e8, 1e7, 1e6,
  1e5, 1e4, 1e3, 1e0, Inf # See Note, 31/01/2022, Calendar, 1400
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
results <- list()
for (i in 5){#seq_along(evnt)) {
  ev <- evnt[[i]]

  #lapply(evnt, function(ev) {
  # Actual Event 56417, need to go earlier
  ev$Events <- ev$Events %>% dplyr::filter(Times < 53000)
  results[[i]] <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
    Pool = pool, NumEnvironments = Environments,
    InteractionMatrices = matr,
    Events = ev,
    PerCapitaDynamics = PerCapitaDynamics,
    DispersalMatrix = DispersalMatrices,
    EliminationThreshold = EliminationThreshold,
    ArrivalDensity = ArrivalDensity,
    MaximumTimeStep = MaximumTimeStep,
    BetweenEventSteps = BetweenEventSteps,
    Verbose = TRUE # FALSE
  )#})
}


# Save the system: ############################################################
print("Save")
save(
  results, file = file.path(
    outputLocation, outputName
  )
)
