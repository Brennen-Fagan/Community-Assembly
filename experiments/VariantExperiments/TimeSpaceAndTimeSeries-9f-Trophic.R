# Introduction: ###############################################################
# Branch of TimeSpaceAndTimeSeries-9f-Diversities.R.
# This file intends to keep track of changes in the trophic structure through
# time.

# Parameters: #################################################################
alsoload <- TRUE # if TRUE, try to load all ColExt files encountered.
# if FALSE, only try to create new ColExt files (and return the outputs).
overwrite <- TRUE # if TRUE, ignore whether a previous file exists.

#datfolders <- dir(pattern = "TSTS_Simulations_")#.+2024-11-19$")
datfolders <- dir(pattern = "TSTS_Simulations_.+2025-01-2.$")
# datfolders <- dir(pattern = "CompareEliminationThresholds$")
cores <- 1 # Parallelization?
#cargs <- as.numeric(commandArgs(TRUE))
#cores <- cargs[1]

# Libraries: ##################################################################
#librarypath <- file.path(".", "Rlibs")
#if (!dir.exists(librarypath)) {
#  dir.create(librarypath, showWarnings = FALSE)
#}
#.libPaths(c(librarypath, .libPaths()))
#
#allLibraryPaths <- .libPaths()

library(dplyr)
library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R")
source("CalculateTrophicStructure.R") # Calculator creator.
source("toCheddar.R") # Updated function.
source("TimeSpaceAndTimeSeries-9-Dictionaries.R") # Need to access intervention
# debug(calculateColExtMetrics)

library(parallel)
library(iterators)
library(doParallel)
library(foreach)

# Functions: ##################################################################

# flattenCEs <- function(CE) {
#   id <- strsplit(
#     strsplit(
#       strsplit(basename(CE$Ellipsis$ParentRun), ".", fixed = TRUE)[[1]][1], # Remove .RData.
#       "_", fixed = TRUE)[[1]][-c(1:2)], # Remove TSTS_Type and split seeds off.
#     "-", fixed = TRUE # Separate out the id values.
#   )
#
#
#   if (length(id) < 3) {
#     # I.e., no intervention.
#     id[[3]] <- rep(NA, 4)
#     id[[4]] <- rep(NA, 2)
#   }
#
#   tidytable::data.table(CE$Events) %>% tidytable::rename(
#     EventType = Type.x,
#     SpeciesType = Type.y
#   ) %>% tidytable::mutate(
#     PostIntervention = if("TimeIntervention" %in% names(CE$Ellipsis)) {
#       Times > CE$Ellipsis$TimeIntervention
#     } else {
#       NA
#     },
#     PoolPatch = id[[1]][1],
#     PoolPatchSeed = id[[2]][1],
#     Interactions = id[[1]][2],
#     InteractionsSeed = id[[2]][2],
#     Events = id[[1]][3],
#     EventsSeed = id[[2]][3],
#     InitialConditions = id[[1]][4],
#     InitialConditionsSeed = id[[2]][4],
#     Dispersal = id[[1]][5],
#     NicheDistance = id[[1]][6],
#     Affinity = id[[1]][7],
#     AffinitySeed = id[[2]][5],
#     InterventionPatchType = id[[3]][1],
#     InterventionPatchSeed = id[[4]][1],
#     InterventionTimeType = id[[3]][2],
#     InterventionTimeSeed = id[[4]][2],
#     InterventionDispersal = id[[3]][3],
#     InterventionNicheDistance = id[[3]][4]
#   )
# }

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
Trophic <- foreach::foreach(
  x = iterators::iter(
    dir(datfolders, full.names = TRUE,
        pattern = "(Simulation|Result|Intervention)")
  ), .packages = c("dplyr", "RMTRCode2")
) %op% {
  .libPaths(c(librarypath, .libPaths()))
  library("dplyr")
  library("RMTRCode2")

  x_properties <- strsplit(basename(x), split = splitchar)
  stopifnot(length(x_properties) == 1#,
            #x_properties[[1]][1] == "TSTS",
            #x_properties[[1]][2] == "Simulation"
  )

  filename <- file.path(
    dirname(x),
    if (flag == "TSTS") {
      paste0(c(x_properties[[1]][1],
               "Trophic",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_Trophic_",
             paste0(x_properties[[1]][5:length(x_properties[[1]])],
                    collapse = "_"))
    } else {
      paste0("Trophic_", x)
    }
  )

  if(!overwrite && file.exists(filename)) {
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
      x_mats <- poolmats[[x_poolind]]$InteractionMatrices
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

    # Build the calculator that returns the trophic structures.
    if (x_properties[[1]][2] == "Intervention") {
      # Detect intervention strategy
      interventionType <-
        as.numeric(strsplit(x_properties[[1]][5], split = "-")[[1]][2])
      interventionTimeDictionary <-
        interventionTimeDictionaryOrigin[interventionType,]
      interventionTimeType <- interventionTimeDictionary$InterventionTimespan

      if (interventionTimeType == 0) {
        # We just need two calculators, before and after.
        calculatorBoth <- lapply(
          list(loaded$Ellipsis$Affinity$EffectiveReproductionRateOld,
               loaded$Ellipsis$Affinity$EffectiveReproductionRateIntervention),
          CalculateTrophicStructure,
          Pool = x_pool,
          NumEnvironments = loaded$NumEnvironments,
          InteractionMatrices = x_mats,
          EliminationThreshold = loaded$Parameters$EliminationThreshold
        )
        calculator <- function(t, y) {
          if (t < loaded$Ellipsis$Affinity$TimeIntervention) {
            calculatorBoth[[1]](y)
          } else {
            calculatorBoth[[2]](y)
          }
        }
      } else {
        # We need to create calculators on the fly.
        interpolation <- interpolateMatrices( #TODO: Cannot take function args.
          loaded$Ellipsis$Affinity$EffectiveReproductionRateOld,
          loaded$Ellipsis$Affinity$EffectiveReproductionRateIntervention,
          switchtime = loaded$Ellipsis$Affinity$TimeIntervention,
          timespan = interventionTimeType
        )
        calculator <- function(t, y) {
          engine <- CalculateTrophicStructure(
            Pool = x_pool,
            EffectiveReproductionRate = interpolation(t),
            NumEnvironments = loaded$NumEnvironments,
            InteractionMatrices = x_mats,
            EliminationThreshold = loaded$Parameters$EliminationThreshold
          )
          return(engine(y))
        }
      }
    } else {
      # Just use the "old" ReproductionRate
      calculatorBoth <- lapply(
        list(loaded$Ellipsis$Affinity$EffectiveReproductionRate),
        CalculateTrophicStructure,
        Pool = x_pool,
        NumEnvironments = loaded$NumEnvironments,
        InteractionMatrices = x_mats,
        EliminationThreshold = loaded$Parameters$EliminationThreshold
      )
      calculator <- function(t = NULL, y) {
        calculatorBoth[[1]](y)
      }
    }

    # Apply the calculator to each time step and aggregate the results.
    # This is likely slow and might be space intensive.
    Trophic <- apply(loaded$Abundance, MARGIN = 1,
                     function(ro) {calculator(ro[1], ro[-1])},
                     simplify = FALSE)

    Trophic <- list(TrophicStructure = Trophic, Ellipsis = list())
    if ("ParentRun" %in% names(loaded$Ellipsis))
      Trophic$Ellipsis$GrandparentRun <- loaded$Ellipsis$ParentRun
    Trophic$Ellipsis$ParentRun <- x

    # So now Trophic contains both neutral and dynamic observed events through
    # time, the appropriate bins to ease plotting, and simulation metadata.
    save(Trophic, file = filename)

    if (alsoload) {
      Trophic # return the object to the foreach loop.
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
  CE$Ellipsis$TimeStart <-
    loaded$Abundance[1, 1] / loaded$ReactionTime

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
      Times < CE$Ellipsis$TimeStart
    ),
    CE$Events
  )
  CE
}

# Recombine into a single set.
ColExt <- c(ColExtBase, ColExtIntervention)

# Save this almost processed object so we don't miss out.
save(ColExt, file = "ColExt9a9_full.RData")

# Flatten the object to facilitate plotting. Fairly fast believe it or not.
ColExt <- tidytable::bind_rows(lapply(ColExt, flattenCEs))

# Save the flat object for combination with the flattened diversities.
save(ColExt, file = "ColExt9a9_flat.RData")

if (cores > 1)
  parallel::stopCluster(clust)
