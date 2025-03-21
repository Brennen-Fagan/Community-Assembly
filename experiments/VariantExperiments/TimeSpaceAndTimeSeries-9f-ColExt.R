# Introduction: ###############################################################
# Branch of TimeSpaceAndTimeSeries-9f-Diversities.R.
# This file intends to keep track of colonization, dynamic extirpations, and
# neutral extirpations by consumer/basal and affinity type.

# Parameters: #################################################################
alsoload <- TRUE # if TRUE, try to load all ColExt files encountered.
# if FALSE, only try to create new ColExt files (and return the outputs).

#datfolders <- dir(pattern = "TSTS_Simulations_")#.+2024-11-19$")
datfolders <- dir(pattern = "TSTS_Simulations_.+2025-01-2.$")
# datfolders <- dir(pattern = "CompareEliminationThresholds$")
cores <- 8 # Parallelization?

# Libraries: ##################################################################
library(dplyr)
library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R") # Abundance metrics.
# debug(calculateColExtMetrics)

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
      loaded <- load(filename)
      stopifnot(length(loaded) == 1)
      loaded <- (get(loaded)) # objects
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

    if (alsoload) {
      ColExt # return the object to the foreach loop.
    }
  }
}



# Now to process into a compact whole.
# We have two types: base and intervention.
# A nice distinguishing feature is whether they have a GrandparentRun
# listed in the Ellipsis argument. If they do, set them aside.
# If they don't, set them to the other side and index them by ParentRun.
# Then we can connect the ones with GrandparentRun attributes to the appropriate
# ParentRun attributes.
# This leaves a problem in the form of the intervention time, however.
# The most correct solution is likely to load up the intervention run to check.
# (Or else, re-run the entire evaluation while including this information.)

# Might need parallel
# Separate out and label, grabbing intervention times.
# ColExtBase <- vector("list")
# ColExtIntervention <- vector("list")
# for (CE in ColExt) {
#   if ("GrandparentRun" %in% names(CE$Ellipsis)) {
#
#     loaded <- load(CE$Ellipsis$ParentRun) # names
#     stopifnot(length(loaded) == 1)
#     loaded <- (get(loaded)) # objects
#
#     CE$Ellipsis$TimeIntervention <-
#       loaded$Ellipsis$Affinity$TimeIntervention / loaded$ReactionTime
#
#     ColExt2$ColExtIntervention <- c(ColExt2$ColExtIntervention, list(CE))
#
#   } else {
#     ColExt2$ColExtBase <- c(ColExt2$ColExtBase, list(CE))
#     names(ColExt2$ColExtBase)[length(ColExt2$ColExtBase)] <-
#       CE$Ellipsis$ParentRun
#   }
# }


ColExtBase <- vector("list")
ColExtIntervention <- vector("list")
CEIs <- unlist(lapply(
  ColExt, function(CE) "GrandparentRun" %in% names(CE$Ellipsis)
))
CEBs <- which(!CEIs)
CEIs <- which(CEIs)

ColExtIntervention <- foreach::foreach(
  x = iterators::iter(
    CEIs
  ), .packages = c("dplyr", "RMTRCode2")
) %op% {
  CE <- ColExt[[x]]
  loaded <- load(CE$Ellipsis$ParentRun) # names
  stopifnot(length(loaded) == 1)
  loaded <- (get(loaded)) # objects

  CE$Ellipsis$TimeIntervention <-
    loaded$Ellipsis$Affinity$TimeIntervention / loaded$ReactionTime

  CE
}

ColExtBase <- ColExt[CEBs]
names(ColExtBase) <-
  unlist(lapply(ColExtBase, function(CE) CE$Ellipsis$ParentRun))


# ColExt <- vector("list")
# Process Intervention CEs, deposit into the CE results.
# for (CE in ColExtIntervention) {
#   CEBase <- ColExtBase[CE$Ellipsis$GrandparentRun]
#   CE$Events <- rbind(
#     CEBase$Events %>% dplyr::filter(
#       Times < CE$Ellipsis$TimeIntervention
#     ),
#     CE$Events
#   )
#   ColExt <- c(ColExt, list(CE))
# }
ColExtIntervention <- foreach::foreach(
  CE = iterators::iter(
    ColExtIntervention
  ), .packages = c("dplyr", "RMTRCode2")
) %op% {
  CEBase <- ColExtBase[[CE$Ellipsis$GrandparentRun]]
  CE$Events <- rbind(
    CEBase$Events %>% dplyr::filter(
      Times < CE$Ellipsis$TimeIntervention
    ),
    CE$Events
  )
  CE
}

# Recombine into a single set.
ColExt <- c(ColExtBase, ColExtIntervention)

# Save this almost processed object so we don't miss out.
save(ColExt, file = "ColExt9a9_full.RData")

# Flatten the object to facilitate plotting.

# Save the flat object for combination with the flattened diversities.

if (cores > 1)
  parallel::stopCluster(clust)
