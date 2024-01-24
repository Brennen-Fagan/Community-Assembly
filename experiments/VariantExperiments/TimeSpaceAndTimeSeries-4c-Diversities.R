datfolders <- c(
  "TSTS_Simulations_1-1-1_2024-01-16",
  "TSTS_Simulations_1-2-1_2024-01-10",
  "TSTS_Simulations_1-2-2_2024-01-10",
  "TSTS_Simulations_1-2-3_2024-01-10",
  "TSTS_Simulations_1-3-1_2024-01-12",
  "TSTS_Simulations_1-3-2_2024-01-12",
  "TSTS_Simulations_1-3-3_2024-01-15"#,
  # "TSTS_Simulations_2-1-6_2024-01-19",
  # "TSTS_Simulations_2-2-6_2024-01-19",
  # "TSTS_Simulations_2-3-6_2024-01-19"#,
  # "TSTS_Simulations_3-1-7_2024-01-22",
  # "TSTS_Simulations_3-2-7_2024-01-22",
  # "TSTS_Simulations_3-3-7_2024-01-22"
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

Diversity <- foreach::foreach(
  x = iterators::iter(
    dir(datfolders, full.names = TRUE, pattern = "Simulation")
  ), .packages = c("dplyr", "RMTRCode2")
) %op% {
  x_properties <- strsplit(basename(x), split = "_")
  stopifnot(length(x_properties) == 1,
            x_properties[[1]][1] == "TSTS",
            x_properties[[1]][2] == "Simulation")
  filename <- file.path(
    dirname(x),
    paste0(x_properties[[1]][1],
           "_Diversity_",
           x_properties[[1]][3])
  )
  
  if(file.exists(filename)) {
    load(filename)
  } else {
    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects
    
    numberOfSpecies <- (ncol(loaded$Abundance) - 1)/loaded$NumEnvironments
    Diversity <- list(
      Diversities = RMTRCode2::Calculate_Diversity(
        loaded,
        nspecies = c(Basal = 34, Consumer = 66) * numberOfSpecies / 100
        # My standard approach for nspecies.
      ),
      Ellipsis = loaded$Ellipsis
    )
    
    save(Diversity, file = filename)
  }
  
  return(Diversity)
}


SpeciesPresence <- foreach::foreach(
  x = iterators::iter(
    dir(datfolders, full.names = TRUE, pattern = "Simulation")
  ), .packages = c("dplyr", "RMTRCode2")
) %op% {
  x_properties <- strsplit(basename(x), split = "_")
  stopifnot(length(x_properties) == 1,
            x_properties[[1]][1] == "TSTS",
            x_properties[[1]][2] == "Simulation")
  filename <- file.path(
    dirname(x),
    paste0(x_properties[[1]][1],
           "_Presence_",
           x_properties[[1]][3])
  )
  
  if(file.exists(filename)) {
    load(filename)
  } else {
    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects
    
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