# Introduction: ###############################################################
# Branch of TimeSpaceAndTimeSeries-9f-Diversities.R.
# This file intends to keep track of colonization, dynamic extirpations, and
# neutral extirpations by consumer/basal and affinity type.

# Parameters: #################################################################
alsoload <- TRUE # if TRUE, try to load all ColExt files encountered.
# if FALSE, only try to create new ColExt files (and return the outputs).
overwrite <- FALSE #TRUE # if TRUE, ignore whether a previous file exists.

datfolders <- dir(pattern = "TSTS_Simulations_.+2025-07-30$")

cargs <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
if (!exists("cores")) {
  if (is.null(cargs) || is.na(cargs)) {
    cores <- 1
  } else {
    cores <- cargs
  }
}

print(paste("Cores", cores))

# Libraries: ##################################################################
librarypath <- file.path(".", "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

library(dplyr)
library(RMTRCode2)

source(file.path("R", "calculateColExtMetrics.R"))

library(parallel)
library(iterators)
library(doParallel)
library(foreach)

# Functions: ##################################################################

flattenCEs <- function(CE) {
  id <- strsplit(
    strsplit(
      # Remove .RData.
      strsplit(basename(CE$Ellipsis$ParentRun), ".", fixed = TRUE)[[1]][1],
      # Remove TSTS_Type and split seeds off.
      "_", fixed = TRUE)[[1]][-c(1:2)],
    # Separate out the id values.
    "-", fixed = TRUE
  )

  if (length(id) < 3) {
    # I.e., no intervention.
    id[[3]] <- rep(NA, 4)
    id[[4]] <- rep(NA, 2)
  }

  tidytable::data.table(CE$Events) %>% tidytable::rename(
    EventType = Type.x,
    SpeciesType = Type.y
  ) %>% tidytable::mutate(
    PostIntervention = if("TimeIntervention" %in% names(CE$Ellipsis)) {
      Times > CE$Ellipsis$TimeIntervention
    } else {
      NA
    },
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

poolmats <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Pool"), function(x) {
    names <- load(x)
    return(c(mget(names), "Dir" = dirname(x), "File" = basename(x)))
  })

# Calculations: ###############################################################
allfiles <- dir(datfolders, full.names = TRUE,
                pattern = "(Simulation|Result|Intervention)")

ColExt <- foreach::foreach(
  id = iterators::iter(1:length(allfiles)),
  x = iterators::iter(allfiles)
) %op% {
  .libPaths(c(librarypath, .libPaths()))
  library("dplyr")
  library("RMTRCode2")

  x_properties <- strsplit(basename(x), split = splitchar)
  stopifnot(length(x_properties) == 1)

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

  if(!overwrite && file.exists(filename)) {
    if (alsoload) {
      loaded <- load(filename)
      stopifnot(length(loaded) == 1)
      loaded <- (get(loaded)) # objects
    }
  } else {
    print(paste(id, filename))
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
        data.frame(
          Species = 1:length(result$Ellipsis$Affinity$SpeciesAffinities),
          Affinity = result$Ellipsis$Affinity$SpeciesAffinities,
          AffinityBins = AffinitiesBinned),
        by = "Species")
    }

    ColExt <- list(Events = ColExt, Ellipsis = list())
    if ("ParentRun" %in% names(loaded$Ellipsis)) {
      # Must be an Intervention Run
      ColExt$Ellipsis$GrandparentRun <- loaded$Ellipsis$ParentRun
      ColExt$Ellipsis$TimeIntervention <-
        loaded$Ellipsis$Affinity$TimeIntervention / loaded$ReactionTime
      ColExt$Ellipsis$TimeStart <-
        loaded$Abundance[1, 1] # / loaded$ReactionTime # Already scaled above.
    }
    ColExt$Ellipsis$ParentRun <- x

    # So now ColExt contains both neutral and dynamic observed events through
    # time, the appropriate bins to ease plotting, and simulation metadata.
    save(ColExt, file = filename)

    if (alsoload) {
      ColExt # return the object to the foreach loop.
    }
  }
}


if (alsoload) {
  # Now to process into a compact whole.
  # We have two types: base and intervention.
  # A nice distinguishing feature is whether they have a GrandparentRun
  # listed in the Ellipsis argument. If they do, set them aside.
  # If they don't, set them to the other side and index them by ParentRun.
  # Then we connect the ones with GrandparentRun attributes to the appropriate
  # ParentRun attributes.
  # This leaves a problem in the form of the intervention time, however.
  # The most correct solution is likely to load the intervention run to check.
  # (Or else, re-run the entire evaluation while including this information.)

  ColExtBase <- vector("list")
  ColExtIntervention <- vector("list")
  CEIs <- unlist(lapply(
    ColExt, function(CE) "GrandparentRun" %in% names(CE$Ellipsis)
  ))
  CEBs <- which(!CEIs)
  CEIs <- which(CEIs)

  ColExtBase <- ColExt[CEBs]
  names(ColExtBase) <-
    unlist(lapply(ColExtBase, function(CE) CE$Ellipsis$ParentRun))

  # Process Intervention CEs, deposit into the CE results.
  ColExtIntervention <- foreach::foreach(
    CEindex = iterators::iter(
      CEIs
    ), .packages = c("dplyr", "RMTRCode2")
  ) %op% {
    CE <- ColExt[[CEindex]]
    CEBase <- ColExtBase[[CE$Ellipsis$GrandparentRun]]
    CE$Events <- rbind(
      CEBase$Events %>% dplyr::filter(
        Times < CE$Ellipsis$TimeStart
      ),
      CE$Events
    )
    CE
  }

  # Recombine into a single set.
  ColExt <- c(ColExtBase, ColExtIntervention)

  # Save this almost processed object so we don't miss out.
  save(ColExt, file = "ColExt10a1_full.RData")

  # Flatten the object to facilitate plotting. Fairly fast believe it or not.
  ColExt <- tidytable::bind_rows(lapply(ColExt, flattenCEs))

  # Save the flat object for combination with the flattened diversities.
  save(ColExt, file = "ColExt10a1_flat.RData")
}

if (cores > 1)
  parallel::stopCluster(clust)
