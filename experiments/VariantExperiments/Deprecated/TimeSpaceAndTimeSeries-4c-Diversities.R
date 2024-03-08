datfolders <- c(
  # "TSTS_Simulations_1-1-1_2024-01-16",
  # "TSTS_Simulations_1-2-1_2024-01-10",
  # "TSTS_Simulations_1-2-2_2024-01-10",
  # "TSTS_Simulations_1-2-3_2024-01-10",
  # "TSTS_Simulations_1-3-1_2024-01-12",
  # "TSTS_Simulations_1-3-2_2024-01-12",
  # "TSTS_Simulations_1-3-3_2024-01-15"#,
  # "TSTS_Simulations_2-1-6_2024-01-19",
  # "TSTS_Simulations_2-2-6_2024-01-19",
  # "TSTS_Simulations_2-3-6_2024-01-19"#,
  # "TSTS_Simulations_3-1-7_2024-01-22",
  # "TSTS_Simulations_3-2-7_2024-01-22",
  # "TSTS_Simulations_3-3-7_2024-01-22"#,
  "Data_2024-01-29"
)

cores <- 3


# Libraries: ##################################################################
library(dplyr)
library(RMTRCode2)

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
          all(x[1] == x)) ){
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
  dir(datfolders, full.names = TRUE, pattern = "PoolMats"), function(x) {
    names <- load(x)
    return(c(mget(names), "Dir" = x))
  })

Diversity <- foreach::foreach(
  x = iterators::iter(
    dir(datfolders, full.names = TRUE, pattern = "(Simulation|Result)")
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
      paste0(x_properties[[1]][1],
             "_Diversity_",
             x_properties[[1]][3])
    } else if (flag == "Data") {
      paste0("TSTS_Diversity_",
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

    if (!is.null(x_pool)) {
      numberOfSpecies <- nrow(x_pool)
    } else {
      numberOfSpecies <- (ncol(loaded$Abundance) - 1) / loaded$NumEnvironments
    }

    if (!is.null(x_pool) && any(grepl(pattern = "Niche", x = names(x_pool)))) {
      # Proceed for each niche, then grab the overall.
      Niches <- grep(pattern = "Niche", x = names(x_pool), value = TRUE)

      Diversity <- lapply(
        Niches, function(colname) {
          x_pool %>% dplyr::group_by(
            dplyr::across(dplyr::all_of(colname))
          ) %>% dplyr::group_modify(
            .f = function(.x, .y) {
              # .x$ID is numeric!
              idcolumns <- c(1, rep(.x$ID, loaded$NumEnvironments) +
                1 +  # Time Column
                numberOfSpecies * rep(1:loaded$NumEnvironments,
                                      each = length(.x$ID)) # Env Indexing
              )
              loaded_subset <- loaded
              loaded_subset$Abundance <- loaded_subset$Abundance[, idcolumns]

              RMTRCode2::Calculate_Diversity(
                loaded,
                nspecies = c(Basal = sum(.x$Type == "Basal"),
                             Consumer = sum(.x$Type == "Consumer"))
                # My standard approach for nspecies.
              )
            }
          )
        }
      )

      names(Diversity) <- Niches

      Diversity$Diversities <- RMTRCode2::Calculate_Diversity(
        loaded,
        nspecies = c(Basal = 34, Consumer = 66) * numberOfSpecies / 100
        # My standard approach for nspecies.
      )

      Diversity$Ellipsis <- loaded$Ellipsis
    } else {
      Diversity <- list(
        Diversities = RMTRCode2::Calculate_Diversity(
          loaded,
          nspecies = c(Basal = 34, Consumer = 66) * numberOfSpecies / 100
          # My standard approach for nspecies.
        ),
        Ellipsis = loaded$Ellipsis
      )
    }

    save(Diversity, file = filename)
  }

  return(Diversity)
}


SpeciesPresence <- foreach::foreach(
  x = iterators::iter(
    dir(datfolders, full.names = TRUE, pattern = "(Simulation|Result)")
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
      paste0(x_properties[[1]][1],
             "_Presence_",
             x_properties[[1]][3])
    } else if (flag == "Data") {
      paste0("TSTS_Presence_",
             paste0(x_properties[[1]][5:length(x_properties[[1]])],
                    collapse = "_"))
    }
  )

  if(file.exists(filename)) {
    load(filename)
  } else {
    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects

    if (!"ReactionTime" %in% names(loaded$Ellipsis)) {
      loaded$Ellipsis$ReactionTime <- loaded$ReactionTime
    }

    SpeciesPresence <- list(
      SpeciesPresences = RMTRCode2::Calculate_Species(
        loaded
      ),
      Ellipsis = loaded$Ellipsis
    )

    save(SpeciesPresence, file = filename)
  }

  return(SpeciesPresence)
}

# Cleanup: ####################################################################
if (exists("clust")) {
  parallel::stopCluster(clust)
}
