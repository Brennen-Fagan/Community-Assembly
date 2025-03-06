# Introduction: ###############################################################
# Branch of TimeSpaceAndTimeSeries-9f-Diversities.R.
# This file intends to keep track of colonization, dynamic extirpations, and
# neutral extirpations by consumer/basal and affinity type.

# Parameters: #################################################################
alsoload <- FALSE # if TRUE, try to load all ColExt files encountered.
# if FALSE, only try to create new ColExt files (and return the outputs).

#datfolders <- dir(pattern = "TSTS_Simulations_")#.+2024-11-19$")
datfolders <- dir(pattern = "TSTS_Simulations_.+2025-01-2.$")
# datfolders <- dir(pattern = "CompareEliminationThresholds$")
cores <- 8 # Parallelization?

# Libraries: ##################################################################
library(dplyr)
library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R") # Abundance metrics.

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
ColExt <- foreach::foreach(
  x = iterators::iter(
    dir(datfolders, full.names = TRUE,
        pattern = "(Simulation|Result|Intervention)")
  ), .packages = c("dplyr", "RMTRCode2")
) %op% {
  x_properties <- strsplit(basename(x), split = splitchar)
  stopifnot(length(x_properties) == 1#,
            #x_properties[[1]][1] == "TSTS",
            #x_properties[[1]][2] == "Simulation"
  )

  filename <- file.path(
    dirname(x),
    if (flag == "TSTS") {
      paste0(c(x_properties[[1]][1],
               "ColExt",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_ColExt_",
             paste0(x_properties[[1]][5:length(x_properties[[1]])],
                    collapse = "_"))
    } else {
      paste0("ColExt_", x)
    }
  )

  if(file.exists(filename)) {
    if (alsoload) {
      load(filename)
    }
  } else {
    print(filename)
    x_dir <- dirname(x)
    x_poolind <- which(unlist(lapply(poolmats, function(y) y$Dir == x_dir)))
    if(length(x_poolind) == 1) {
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
      loaded$Ellipsis$Timescale <- "Characteristic"
    }


    ColExt <- calculateColExtMetrics(loaded)

    # Add in Traits.
    if (!is.null(x_pool)) {
      ColExt <- ColExt %>% dplyr::left_join(x_pool, by = c("Species" = "ID"))
    }

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

      ColExt <- ColExt %>% dplyr::left_join(
        data.frame(Species = 1:length(result$Ellipsis$Affinity$SpeciesAffinities),
                   Affinity = result$Ellipsis$Affinity$SpeciesAffinities,
                   AffinityBins = AffinitiesBinned),
        by = "Species")
    }

    ColExt <- list(Events = ColExt, Ellipsis = list())
    if ("ParentRun" %in% names(loaded$Ellipsis))
      ColExt$Ellipsis$GrandparentRun <- loaded$Ellipsis$ParentRun
    ColExt$Ellipsis$ParentRun <- x

    # So now ColExt contains both neutral and dynamic observed events through
    # time, the appropriate bins to ease plotting, and simulation metadata.
    save(ColExt, file = filename)
  }
}

if (cores > 1)
  parallel::stopCluster(clust)
