# Introduction: ###############################################################
# Sequel to TimeSpaceAndTimeSeries-6c-Diversities.R.
# We're customising the previous version to the new spec of 8.
# We're also going to try to incorporate the various diversity metrics that
# have been additionally requested.

# Parameters: #################################################################
alsoload <- FALSE # if TRUE, try to load all diversity files encountered.
# if FALSE, only try to create new diversity files (and return the outputs).
overwrite <- FALSE

# If TimeIntervention is not present, should we calculate it?
substituteTimeIntervention <- function(times) {
  mean(c(median(times), 1/2 * max(times)))
}
# substituteTimeIntervention <- NULL

datfolders <- dir(pattern = "TSTS_Simulations_.+2025-07-30$")

cargs <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
if (!exists("cores")) {
  if (is.null(cargs) || is.na(cargs)) {
    cores <- 1
  } else {
    cores <- cargs
  }
}

preferredTimestep <- 10 # Characteristic Time Scale Units
# Previously Event rate was ~1/CTU, now it's more like ~0.1/CTU in theory.
# (This seems reasonably close in practice looking at 1 example and observing
#  that we are talking about events of a single type at a time.)
# Reminder, smaller values mean more data, mean more storage.

# Minimum number of rows per event is 2 (was 10).
# Default timesteps per event is 5 (was 1) in Simulation Time Scale Units.
# In the example I looked at, we've generally slowed things down so that the
# default/maximum timestep of 5 is hit consistently, ~97% of steps.
# In the example I looked at, the STU = 38 * CTU.
#
# The preferred Timestep should be larger than this number.
# The current value is about 1 recording per 2 events (of ext. and/or imm.).

# Libraries: ##################################################################
directory <- '.'
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

library(dplyr)
library(RMTRCode2)
library(betapart) # Dependency for the diversity calculation

# source("TimeSpaceAndTimeSeries-0-Functions.R") # Abundance metrics.
source(file.path("R", "Calculate_Species.R"))

library(parallel)
library(iterators)
library(doParallel)
library(foreach)

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

poolmats <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Pool"), function(x) {
    names <- load(x)
    return(c(mget(names), "Dir" = dirname(x), "File" = basename(x)))
  })

# Calculations: ###############################################################
allfiles <- dir(datfolders, full.names = TRUE,
                pattern = "(Simulation|Result|Intervention)")

Diversity <- foreach::foreach(
  id = iterators::iter(1:length(allfiles)),
  x = iterators::iter(allfiles)
  #, .packages = c("dplyr", "RMTRCode2")
) %op% {
  directory <- '.'
  librarypath <- file.path(directory, "Rlibs")
  .libPaths(c(librarypath, .libPaths()))
  library(RMTRCode2); library(dplyr)

  x_properties <- strsplit(basename(x), split = splitchar)
  stopifnot(length(x_properties) == 1#,
            #x_properties[[1]][1] == "TSTS",
            #x_properties[[1]][2] == "Simulation"
  )

  filename <- file.path(
    dirname(x),
    if (flag == "TSTS") {
      paste0(c(x_properties[[1]][1],
               "Diversity",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_Diversity_",
             paste0(x_properties[[1]][5:length(x_properties[[1]])],
                    collapse = "_"))
    } else {
      paste0("Diversity_", x)
    }
  )

  if(!overwrite && file.exists(filename)) {
    if (alsoload) {
      load(filename)
    }
  } else {
    print(paste(id, filename))
    x_dir <- dirname(x)
    x_poolind <- which(unlist(lapply(poolmats, function(y) y$Dir == x_dir)))
    if(any(x_poolind)) {
      x_pool <- poolmats[[x_poolind]]$Pool
    } else {
      x_pool <- NULL
    }

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


    loaded$Abundance <- thinAbundanceTimes(
      abundance = loaded$Abundance,
      threshold = loaded$Parameters$EliminationThreshold,
      times = sort(unique(c(
        # Start time
        min(loaded$Abundance[, 1]),
        # Base sequence, aligned to CTU time scale.
        seq(
          from = ceiling(loaded$Abundance[1, 1]/preferredTimestep)*preferredTimestep,
          to = max(loaded$Abundance[, 1]),
          by = max(diff(loaded$Abundance[, 1]), preferredTimestep)
        ),
        # Intervention times, aligned s.t. 0 = intervention.
        loaded$Ellipsis$Affinity$TimeIntervention +
          c(0:10, 2*(6:10), 20+3*(1:10)) # by 1s til 10, 2s til 20, 3s til 50.
      )))
    )

    if (exists("x_pool") && !is.null(x_pool)) {
      numberOfSpecies <- nrow(x_pool)
    } else {
      numberOfSpecies <- (ncol(loaded$Abundance) - 1) / loaded$NumEnvironments
    }


    # Bad implementation (in a few ways!); bins presences in 0.1s of time units.
    Presence <- RMTRCode2::Calculate_Species(loaded, bintimes = FALSE)
    if (nrow(Presence) == 0) {
      warning("No species presences detected at this resolution.")
    }

    # Major edit that is somewhat of a backslide.
    # Instead of programmatically using the named different columns of the
    # pool, we're going to be working off of a list entry. For speed of coding,
    # I'll be treating it as a nx1 vector, rather than a possibly named matrix.
    if ("SpeciesAffinities" %in% names(loaded$Ellipsis$Affinity)) {
      # Identify Niche Cuts. If discrete, this is by value. If continuous, or
      # there are many bins, then this is by binning.
      AffinitiesBinned <-
        if (length(unique(loaded$Ellipsis$Affinity$SpeciesAffinities)) >= 5) {
          cut(loaded$Ellipsis$Affinity$SpeciesAffinities,
              breaks = max(ceiling((loaded$NumEnvironments + 1)/2), 5))
        } else {
          loaded$Ellipsis$Affinity$SpeciesAffinities
        }

      # If Basal/Consumer is present, we can break it into the constituents.
      if ("Type" %in% names(x_pool)) {
        AffinitiesBinned <- paste0(x_pool$Type, "_", AffinitiesBinned)
      }

      DiversityNiche <- lapply(
        unique(AffinitiesBinned),
        function(AffinityType) {
          # Identify the subset of the abundance matrix that has the type.
          idcolumns <- c(
            1,
            rep(which(AffinitiesBinned == AffinityType),
                loaded$NumEnvironments) + 1 +  # Time Column
              numberOfSpecies * rep(
                (1:loaded$NumEnvironments) - 1,
                each = sum(AffinitiesBinned == AffinityType)) # Env Indexing
          )
          loaded_subset <- loaded
          loaded_subset$Abundance <- loaded_subset$Abundance[, idcolumns]

          if ("Size" %in% names(x_pool))
            # Note the (-)1 is there to preserve (resp., remove) the times.
            sizes_subset <- x_pool$Size[idcolumns[-1]]

          Diversities <- calculateDiversityMetrics(
            abundance = loaded_subset$Abundance,
            nspecies = length(idcolumns) - 1,
            nenvironments = loaded$NumEnvironments,
            sizes = if ("Size" %in% names(x_pool)) sizes_subset
          )

          Diversities$Subset <- AffinityType
          return(Diversities)
        }) %>% dplyr::bind_rows()

      # Add in the relevant trait information.
      if (nrow(Presence) > 0) {
        Presence <- Presence %>% dplyr::left_join(
          y = data.frame(
            Species = 1:length(loaded$Ellipsis$Affinity$SpeciesAffinities),
            Size = if ("Size" %in% names(x_pool)) x_pool$Size,
            Type = if ("Type" %in% names(x_pool)) x_pool$Type,
            Affinity = loaded$Ellipsis$Affinity$SpeciesAffinities
          ), by = "Species"
        )
      }
    }

    DiversityAll <- calculateDiversityMetrics(
      abundance = loaded$Abundance,
      nspecies = numberOfSpecies,
      nenvironments = loaded$NumEnvironments,
      sizes = if ("Size" %in% names(x_pool)) x_pool$Size
    )
    DiversityAll$Subset <- NA

    Diversity <-
      list(
        Diversity =
          dplyr::bind_rows(DiversityAll,
                           if (exists("DiversityNiche")) DiversityNiche),
        Presence = Presence,
        Ellipsis = loaded$Ellipsis
      )

    if ("ParentRun" %in% names(Diversity$Ellipsis))
      Diversity$Ellipsis$GrandparentRun <- Diversity$Ellipsis$ParentRun
    Diversity$Ellipsis$ParentRun <- x

    # So now Diversity contains summary statistics, presence/absence values,
    # and simulation metadata.
    save(Diversity, file = filename)
  }
}

if (cores > 1)
  parallel::stopCluster(clust)
