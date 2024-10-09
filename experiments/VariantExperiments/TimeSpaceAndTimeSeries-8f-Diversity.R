# Introduction: ###############################################################
# Sequel to TimeSpaceAndTimeSeries-6c-Diversities.R.
# We're customising the previous version to the new spec of 8.
# We're also going to try to incorporate the various diversity metrics that
# have been additionally requested.

# Parameters: #################################################################
datfolders <- dir(pattern = "TSTS_Simulations_")
cores <- 4 # Parallelization?
rows_per_event <- 1 # Maximum rows per event (base resolution) is 10.

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
Diversity <- foreach::foreach(
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
               "Diversity",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_Diversity_",
             paste0(x_properties[[1]][5:length(x_properties[[1]])],
                    collapse = "_"))
    }
  )

  if(file.exists(filename)) {
    load(filename)
  } else {
    print(filename)
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
    if (loaded$Ellipsis$Timescale == "Simulation") {
      loaded$Events$Times <-
        loaded$Events$Times / loaded$Ellipsis$ReactionTime
      loaded$Abundance[, 1] <-
        loaded$Abundance[, 1] / loaded$Ellipsis$ReactionTime
      loaded$Ellipsis$Timescale <- "Characteristic"
    }

    loaded$Abundance <- thinAbundanceEqualTimeSteps(
      abundance = loaded$Abundance,
      # events = loaded$Events, # unnecessary here.
      threshold = loaded$Parameters$EliminationThreshold,
      preferredTimeStep = rows_per_event, # Abusing the characteristic t scale.
      preferredStart = ceiling(loaded$Abundance[1, 1]) # Must be at or after.
    )

    if (!is.null(x_pool)) {
      numberOfSpecies <- nrow(x_pool)
    } else {
      numberOfSpecies <- (ncol(loaded$Abundance) - 1) / loaded$NumEnvironments
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

          Diversities <- calculateDiversityMetrics(
            loaded_subset$Abundance,
            length(idcolumns) - 1, loaded$NumEnvironments
          )

          Diversities$Subset <- AffinityType
          return(Diversities)
        }) %>% dplyr::bind_rows()
    }
    DiversityAll <- calculateDiversityMetrics(
      loaded$Abundance, numberOfSpecies, loaded$NumEnvironments
    )
    DiversityAll$Subset <- NA

    Diversity <-
      list(
        Diversity =
          dplyr::bind_rows(DiversityAll,
                           if (exists("DiversityNiche")) DiversityNiche),
        Ellipsis = loaded$Ellipsis
      )

    if ("ParentRun" %in% names(Diversity$Ellipsis))
      Diversity$Ellipsis$GrandparentRun <- Diversity$Ellipsis$ParentRun
    Diversity$Ellipsis$ParentRun <- x

    save(Diversity, file = filename)
  }
}

if (cores > 1)
  parallel::stopCluster(clust)
