# Introduction: ###############################################################
# Sequel to TimeSpaceAndTimeSeries-4c-Diversity.R.
# We're customising the previous version to the new spec of 6a.
# E.g. the preference for "affinity" instead of "niche".
#
# As always, please see the previous files for some design choices,
# although we aim to improve design at each stage.

# Parameters: #################################################################
# datfolders <- c(
#   # "TSTS_Simulations_1-1_1-1_2024-02-13",
#   # "TSTS_Simulations_1-1_2-2_2024-02-14",
#   # "TSTS_Simulations_2-1_2-2_2024-02-14",
#   # "TSTS_Simulations_10-1_2-2_2024-02-15",
#   # "TSTS_Simulations_6-1_2-2_2024-02-15"
#   # "TSTS_Simulations_11-1_3-3_2024-02-23",
#   # "TSTS_Simulations_11-1_4-4_2024-02-23"
#   # "TSTS_Simulations_17-1_5-5_2024-03-08",
#   # "TSTS_Simulations_18-1_6-6_2024-03-08"
#   # "TSTS_Simulations_19-1_7-7_2024-03-08",
#   # "TSTS_Simulations_20-1_8-8_2024-03-08"
#   #"TSTS_Simulations_18-1_6-6_2024-05-23"
#   # "TSTS_Simulations_138-1_10-10_2024-08-22"
#   "TSTS_Simulations_158-1_9-9_2024-08-22"
# )
datfolders <- dir(pattern = "TSTS_Simulations_")

cores <- 4

rows_per_event <- 1

# Libraries: ##################################################################
library(dplyr)
library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R") # Abundance metrics.

library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

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

    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects

    if (!"ReactionTime" %in% names(loaded$Ellipsis)) {
      loaded$Ellipsis$ReactionTime <- loaded$ReactionTime
    }

    if (loaded$Ellipsis$Timescale == "Simulation") {
      loaded$Events$Times <- loaded$Events$Times / loaded$Ellipsis$ReactionTime
      loaded$Abundance[, 1] <- loaded$Abundance[, 1] / loaded$Ellipsis$ReactionTime
      loaded$Ellipsis$Timescale <- "Characteristic"
    }

    if (!is.null(x_pool)) {
      numberOfSpecies <- nrow(x_pool)
    } else {
      numberOfSpecies <- (ncol(loaded$Abundance) - 1) / loaded$NumEnvironments
    }

    if (!is.null(x_pool) &&
        any(grepl(pattern = "(Niche|Affinity)", x = names(x_pool)))) {
      # Proceed for each niche, then grab the overall.
      Niches <- grep(pattern = "(Niche|Affinity)", x = names(x_pool),
                     value = TRUE)

      # This is the single niche implementation!!!
      # We may want to, in the future, group_by all niches simultaneously.
      # Alternatively, we may want to discretize the niches first.

      x_pool <- x_pool %>% dplyr::mutate(
        dplyr::across(
          dplyr::all_of(Niches),
          .fns = ~ if(length(unique(.x)) >
                      ceiling((loaded$NumEnvironments + 1)/2)) {
            # Where'd the formula come from? We have 10 patches, arranged in
            # a ring. Assume that there are two patches that are maximally
            # distant in niche space and the two paths are monotonic.
            # These paths are of length n / 2 + 1 = 6. If we had 11, then the
            # extra patch would be on adjacent to one of the extremes (w.l.o.g.)
            # but then either one could be the extreme. So 11 also => 6.
            cut(.x, breaks = max(ceiling((loaded$NumEnvironments + 1)/2), 2))
          } else {
            .x
          }
        )
        )

      Diversity <- lapply(
        Niches, function(colname) {
          x_pool %>% dplyr::group_by(
            dplyr::across(dplyr::all_of(colname))
          ) %>% dplyr::group_map(
            .f = function(.x, .y) {
              # .x$ID is numeric!
              idcolumns <- c(
                1,
                rep(.x$ID, loaded$NumEnvironments) +
                  1 +  # Time Column
                  numberOfSpecies * rep((1:loaded$NumEnvironments) - 1,
                                        each = length(.x$ID)) # Env Indexing
              )
              loaded_subset <- loaded
              loaded_subset$Abundance <- loaded_subset$Abundance[, idcolumns]

              divs <- RMTRCode2::thinAndCalculateDiversities(#Calculate_Diversity(
                loaded_subset,
                nspecies = c(Basal = sum(.x$Type == "Basal"),
                             Consumer = sum(.x$Type == "Consumer")),
                # My standard approach for nspecies.
                preferred_rows_per_event = rows_per_event,
                divide_time_by = 1
              )

              # thinAndCalculateDiversities yields lists
              divs$NicheValues <- .y
              return(divs)
            }
          )
        }
      )

      names(Diversity) <- Niches

      Diversity$Diversities <-
        RMTRCode2::thinAndCalculateDiversities(#Calculate_Diversity(
          loaded,
          nspecies = table(x_pool$Type),
          # My standard approach for nspecies.
          preferred_rows_per_event = rows_per_event,
          divide_time_by = 1
        )

      Diversity$Ellipsis <- loaded$Ellipsis
    } else {
      Diversity <- list(
        Diversities = RMTRCode2::thinAndCalculateDiversities(#Calculate_Diversity(
          loaded,
          nspecies = c(Basal = 34, Consumer = 66) * numberOfSpecies / 100,
          # My standard approach for nspecies.
          preferred_rows_per_event = rows_per_event,
          divide_time_by = 1
        ),
        Ellipsis = loaded$Ellipsis
      )
    }

    save(Diversity, file = filename)
  }

  return(Diversity)
}

DivAbund <- foreach::foreach(
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
               "DivAbund",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_DivAbund_",
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

    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects

    if (!"ReactionTime" %in% names(loaded$Ellipsis)) {
      loaded$Ellipsis$ReactionTime <- loaded$ReactionTime
    }

    if (loaded$Ellipsis$Timescale == "Simulation") {
      loaded$Events$Times <- loaded$Events$Times / loaded$Ellipsis$ReactionTime
      loaded$Abundance[, 1] <- loaded$Abundance[, 1] / loaded$Ellipsis$ReactionTime
      loaded$Ellipsis$Timescale <- "Characteristic"
    }

    if (!is.null(x_pool)) {
      numberOfSpecies <- nrow(x_pool)
    } else {
      numberOfSpecies <- (ncol(loaded$Abundance) - 1) / loaded$NumEnvironments
    }

    if (!is.null(x_pool) &&
        any(grepl(pattern = "(Niche|Affinity)", x = names(x_pool)))) {
      # Proceed for each niche, then grab the overall.
      Niches <- grep(pattern = "(Niche|Affinity)", x = names(x_pool),
                     value = TRUE)

      # This is the single niche implementation!!!
      # We may want to, in the future, group_by all niches simultaneously.
      # Alternatively, we may want to discretize the niches first.

      x_pool <- x_pool %>% dplyr::mutate(
        dplyr::across(
          dplyr::all_of(Niches),
          .fns = ~ if(length(unique(.x)) >
                      ceiling((loaded$NumEnvironments + 1)/2)) {
            # Where'd the formula come from? We have 10 patches, arranged in
            # a ring. Assume that there are two patches that are maximally
            # distant in niche space and the two paths are monotonic.
            # These paths are of length n / 2 + 1 = 6. If we had 11, then the
            # extra patch would be on adjacent to one of the extremes (w.l.o.g.)
            # but then either one could be the extreme. So 11 also => 6.
            cut(.x, breaks = max(ceiling((loaded$NumEnvironments + 1)/2), 2))
          } else {
            .x
          }
        )
      )

      DivAbund <- lapply(
        Niches, function(colname) {
          x_pool %>% dplyr::group_by(
            dplyr::across(dplyr::all_of(colname))
          ) %>% dplyr::group_map(
            .f = function(.x, .y) {
              # .x$ID is numeric!
              idcolumns <- c(
                1,
                rep(.x$ID, loaded$NumEnvironments) +
                  1 +  # Time Column
                  numberOfSpecies * rep((1:loaded$NumEnvironments) - 1,
                                        each = length(.x$ID)) # Env Indexing
              )
              loaded_subset <- loaded
              loaded_subset$Abundance <- loaded_subset$Abundance[, idcolumns]

              divs <- thinAbundance(
                loaded_subset$Abundance, loaded_subset$Events,
                loaded_subset$Parameters$EliminationThreshold,
                rows_per_event
              ) %>% calculateAbundanceMetrics(
                c(Basal = sum(.x$Type == "Basal"),
                  Consumer = sum(.x$Type == "Consumer")),
                loaded_subset$NumEnvironments
              )

              # thinAndCalculateDiversities yields lists
              divs$NicheValues <- .y
              return(divs)
            }
          )
        }
      )

      names(DivAbund) <- Niches

      DivAbund$DivAbund <- thinAbundance(
        loaded$Abundance, loaded$Events,
        loaded$Parameters$EliminationThreshold,
        rows_per_event
      ) %>% calculateAbundanceMetrics(
        table(x_pool$Type),
        loaded$NumEnvironments
      )

      DivAbund$Ellipsis <- loaded$Ellipsis
    } else {
      DivAbund <- list(
        DivAbund = thinAbundance(
          loaded$Abundance, loaded$Events,
          loaded$Parameters$EliminationThreshold,
          rows_per_event
        ) %>% calculateAbundanceMetrics(
          c(Basal = 34, Consumer = 66) * numberOfSpecies / 100,
          loaded$NumEnvironments
        ),
        Ellipsis = loaded$Ellipsis
      )
    }

    save(DivAbund, file = filename)
  }

  return(DivAbund)
}

SpeciesPresence <- foreach::foreach(
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
               "Presence",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_Presence_",
             paste0(x_properties[[1]][5:length(x_properties[[1]])],
                    collapse = "_"))
    }
  )

  if(file.exists(filename)) {
    load(filename)
  } else {
    x_dir <- dirname(x)
    x_poolind <- which(unlist(lapply(poolmats, function(y) y$Dir == x_dir)))
    if(any(x_poolind)) {
      x_pool <- poolmats[[x_poolind]]$Pool
    } else {
      x_pool <- NULL
    }

    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects

    if (!"ReactionTime" %in% names(loaded$Ellipsis)) {
      loaded$Ellipsis$ReactionTime <- loaded$ReactionTime
    }

    if (loaded$Ellipsis$Timescale == "Simulation") {
      loaded$Events$Times <- loaded$Events$Times / loaded$Ellipsis$ReactionTime
      loaded$Abundance[, 1] <- loaded$Abundance[, 1] / loaded$Ellipsis$ReactionTime
      loaded$Ellipsis$Timescale <- "Characteristic"
    }

    SpeciesPresence <- list(
      SpeciesPresences = RMTRCode2::Calculate_Species(
        loaded, bintimes = TRUE # Not well implemented: only two states.
      ) %>% dplyr::filter(Abundance > loaded$Parameters$EliminationThreshold),
      Ellipsis = loaded$Ellipsis
    )

    if ("TimeFloor" %in% names(SpeciesPresence$SpeciesPresences)) {
      SpeciesPresence$SpeciesPresences <-
        SpeciesPresence$SpeciesPresences %>% dplyr::rename(
          Time = TimeFloor
        )
    }

    if (!is.null(x_pool)) {
      # Need to assign patch and species affinities in order to know how
      # well they are pairing with each other.
      SpeciesPresence$SpeciesPresences <-
        SpeciesPresence$SpeciesPresences %>% dplyr::left_join(
          x_pool %>% dplyr::select(
            ID, Size, Type, dplyr::starts_with("Affinity")
          ),
          by = c("Species" = "ID")
        )
    }


    affinityNames <-
      c("PatchAffinities", "PatchAffinitiesOld", "PatchAffinitiesIntervention")
    interventionTime <-
      if ("TimeIntervention" %in% names(loaded$Ellipsis$Affinity)) {
        loaded$Ellipsis$Affinity$TimeIntervention
      } else {
        # We'll take a proxy of the first event...
        min(loaded$Events$Times)
      }

    if (affinityNames[1] %in% names(loaded$Ellipsis$Affinity)) {
      SpeciesPresence$SpeciesPresences$EnvAffinity <-
        loaded$Ellipsis$Affinity[[affinityNames[1]]][
          SpeciesPresence$SpeciesPresences$Environment, 1
          ]
    } else if (
      affinityNames[2] %in% names(loaded$Ellipsis$Affinity) &&
      affinityNames[3] %in% names(loaded$Ellipsis$Affinity)
    ) {
      SpeciesPresence$SpeciesPresences <- SpeciesPresence$SpeciesPresences %>%
        dplyr::mutate(
          EnvAffinity = dplyr::case_when(
            Time < interventionTime |
              (!Environment %in% loaded$Ellipsis$Affinity$PatchInterventions) ~
              loaded$Ellipsis$Affinity[[affinityNames[2]]][
                Environment, 1
                ],
            Time >= interventionTime &
              Environment %in% loaded$Ellipsis$Affinity$PatchInterventions ~
              loaded$Ellipsis$Affinity[[affinityNames[3]]][
                Environment, 1
                ],
            TRUE ~ NA_real_
          )
        )
    } else if (!is.null(x_pool)) {
      # Only stores the base patch "affinities".
      SpeciesPresence$SpeciesPresences$EnvAffinity <-
        poolmats[[x_poolind]]$PatchAffinities[
          SpeciesPresence$SpeciesPresences$Environment,
          ]
    }

    # Unsure where grouping was coming from, but this appears to be inflating
    # file sizes.
    SpeciesPresence$SpeciesPresences <-
      SpeciesPresence$SpeciesPresences %>% dplyr::ungroup()

    save(SpeciesPresence, file = filename)
  }

  return(SpeciesPresence)
}

# Cleanup: ####################################################################
if (exists("clust")) {
  parallel::stopCluster(clust)
}
