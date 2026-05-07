# Introduction: ###############################################################
# Discovered a need for the specific intervention times.
# Borrowing from 10f-Diversity to collect all of the intervention times
# using the "now standard" substituteTimeIntervention.
# In the future, we might need to consider dynamically grabbing it from
# each set of, or each individual!, intervention file, which would be messier.

# Parameters: #################################################################
# If TimeIntervention is not present, should we calculate it?
substituteTimeIntervention <- function(times) {
  mean(c(median(times), 1/2 * max(times)))
}
# substituteTimeIntervention <- NULL

datfolders <- dir(pattern = "TSTS_Simulations_.+2025-07-30$")
# datfolders <- dir(pattern = "TSTS_Simulations_.+2025-01-2[0-9]$")

cargs <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
if (!exists("cores")) {
  if (is.null(cargs) || is.na(cargs)) {
    cores <- 1
  } else {
    cores <- cargs
  }
}

# Libraries: ##################################################################
directory <- '.'
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

library(tidytable)

library(parallel)
library(iterators)
library(doParallel)
library(foreach)

labelInterventionTime <- function(lo) {
  id <- strsplit(
    strsplit(
      # Remove .RData.
      strsplit(basename(lo$Ellipsis$ID), ".", fixed = TRUE)[[1]][1],
      # Split seeds off.
      "_", fixed = TRUE)[[1]],
    # Separate out the id values.
    "-", fixed = TRUE
  )

  if (length(id) < 3) {
    # I.e., no intervention.
    id[[3]] <- rep(NA, 4)
    id[[4]] <- rep(NA, 2)
  }

  tidytable::data.table(
    TimeIntervention = lo$Ellipsis$Affinity$TimeIntervention,
    PoolPatch = id[[1]][1],
    PoolPatchSeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    InitialConditions = id[[1]][4],
    InitialConditionsSeed = id[[2]][4],
    Dispersal = id[[1]][5],
    NicheDistance = id[[1]][6],
    SpeciesAffinity = id[[1]][7],
    SpeciesAffinitySeed = id[[2]][5],
    PatchAffinity = id[[1]][8],
    PatchAffinitySeed = id[[2]][6],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
  )
}

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

# Load Data: ##################################################################
datfolders_properties <- strsplit(datfolders, split = "_")
if ( length(datfolders_properties) > 1 &&
     with(list(x = unlist(lapply(datfolders_properties, length))),
          any(x[1] != x)) ){
  stop("Differing folder types implies differing file types.")
}

flag <- datfolders_properties[[1]][1]
if (flag == "TSTS") {
  splitchar <- "_"
} else if (flag == "Data") {
  splitchar <- "-"
} else {
  stop("Folder type not recognized.")
}

# Calculations: ###############################################################
allfiles <- dir(datfolders, full.names = TRUE,
                pattern = "Simulation")#"(Simulation|Result|Intervention)")

Interventions <- foreach::foreach(
  id = iterators::iter(1:length(allfiles)),
  x = iterators::iter(allfiles)
  #, .packages = c("dplyr", "RMTRCode2")
) %op% {
  directory <- '.'
  librarypath <- file.path(directory, "Rlibs")
  .libPaths(c(librarypath, .libPaths()))
  library(tidytable)

  x_properties <- strsplit(basename(x), split = splitchar)
  stopifnot(length(x_properties) == 1#,
            #x_properties[[1]][1] == "TSTS",
            #x_properties[[1]][2] == "Simulation"
  )

  print(paste(id, x))
  x_dir <- dirname(x)

  # Load result to analyse.
  loaded <- load(x) # names
  stopifnot(length(loaded) == 1)
  loaded <- (get(loaded)) # objects

  # Unify format, double check time scale and make sure on same time scale.
  if (!"ReactionTime" %in% names(loaded$Ellipsis)) {
    loaded$Ellipsis$ReactionTime <- loaded$ReactionTime
  }
  if (!"Timescale" %in% names(loaded$Ellipsis) ||
      loaded$Ellipsis$Timescale == "Simulation") {
    loaded$Events$Times <-
      loaded$Events$Times / loaded$Ellipsis$ReactionTime
    loaded$Abundance[, 1] <-
      loaded$Abundance[, 1] / loaded$Ellipsis$ReactionTime
    if ("TimeIntervention" %in% names(loaded$Ellipsis$Affinity)) {
      loaded$Ellipsis$Affinity$TimeIntervention <-
        loaded$Ellipsis$Affinity$TimeIntervention /
        loaded$Ellipsis$ReactionTime
    }
    loaded$Ellipsis$Timescale <- "Characteristic"
  }

  if (!"TimeIntervention" %in% names(loaded$Ellipsis$Affinity) &&
      !is.null(substituteTimeIntervention)) {
    loaded$Ellipsis$Affinity$TimeIntervention <-
      substituteTimeIntervention(loaded$Events$Times)
  }

  labelInterventionTime(loaded)
}

if (cores > 1)
  parallel::stopCluster(clust)


InterventionTimes <- tidytable::bind_rows(Interventions)

save(InterventionTimes, file = "TSTS_Interventions_10a1.RData")
