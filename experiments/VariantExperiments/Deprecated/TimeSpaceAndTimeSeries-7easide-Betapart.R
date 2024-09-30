# Introduction: ###############################################################
# Sequel to TimeSpaceAndTimeSeries-6c-Diversity.R.
# We're customising the previous version to the new spec of 7b and 7d and to
# now use the betapart package.

# Parameters: #################################################################
datfolders <- dir(pattern = "TSTS_Simulations_")
rows_per_event <- 1

cores <- 1

# Libraries: ##################################################################
library(dplyr)
library(RMTRCode2)
library(betapart)

source("TimeSpaceAndTimeSeries-0-Functions.R")

library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

runBetapart <- function(loaded, numberOfSpecies) {
  # out: Dataframe with 6 columns: Time1, Time2, Env, Measurement, Value
  list(Betapart = dplyr::bind_rows(
    # BraySpace = 
      if(loaded$NumEnvironments > 1) dplyr::bind_rows(lapply(
      1:nrow(loaded$Abundance), function(i, a) {
        outs <- betapart::beta.multi.abund(
          matrix(a[i, -1], byrow = TRUE, ncol = numberOfSpecies)
        )
        data.frame(
          Time1 = a[i, 1][[1]],
          Time2 = a[i, 1][[1]],
          Env = "Multi", 
          Measurement = names(outs),
          Value = unlist(outs, use.names = FALSE),
          stringsAsFactors = FALSE
        )
      }, a = loaded$Abundance)),
    # BrayTimePatch = 
    dplyr::bind_rows(lapply(
      2:nrow(loaded$Abundance), function(i, a) {
        dplyr::bind_rows(
          lapply(1:loaded$NumEnvironments, function(j) {
            outs <- betapart::beta.multi.abund(
              a[c(-1,0)+i, 
                1+numberOfSpecies*(j - 1)+1:numberOfSpecies]
            )
            data.frame(
              Time1 = a[i - 1, 1][[1]],
              Time2 = a[i, 1][[1]],
              Env = toString(j), 
              Measurement = names(outs),
              Value = unlist(outs, use.names = FALSE),
              stringsAsFactors = FALSE
            )
          })
        )
      }, a = loaded$Abundance)),
    # JaccardSpace = 
      if(loaded$NumEnvironments > 1) dplyr::bind_rows(lapply(
      1:nrow(loaded$Abundance), function(i, a) {
        outs <- betapart::beta.multi(
          # as.numeric but without ruining matrix struct.
          1*(matrix(a[i, -1] > 0, byrow = TRUE, ncol = numberOfSpecies)),
          index.family = "jaccard"
        )
        data.frame(
          Time1 = a[i, 1][[1]],
          Time2 = a[i, 1][[1]],
          Env = "Multi", 
          Measurement = names(outs),
          Value = unlist(outs, use.names = FALSE),
          stringsAsFactors = FALSE
        )
      }, a = loaded$Abundance)),
    # JaccardTimePatch = 
    dplyr::bind_rows(lapply(
      2:nrow(loaded$Abundance), function(i, a) {
        dplyr::bind_rows(
          lapply(1:loaded$NumEnvironments, function(j) {
            outs <- betapart::beta.multi(
              1*(a[c(-1,0)+i, # as.numeric but without ruining matrix struct.
                1+numberOfSpecies*(j - 1)+1:numberOfSpecies] > 0),
              index.family = "jaccard"
            )
            data.frame(
              Time1 = a[i - 1, 1][[1]],
              Time2 = a[i, 1][[1]],
              Env = toString(j), 
              Measurement = names(outs),
              Value = unlist(outs, use.names = FALSE),
              stringsAsFactors = FALSE
            )
          })
        )
      }, a = loaded$Abundance)),
    # JaccardTimeRegion = 
    dplyr::bind_rows(lapply(
      2:nrow(loaded$Abundance), function(i, a) {
        outs <- betapart::beta.temp(
          # as.numeric but without ruining matrix struct.
          1*(matrix(a[i - 1, -1] > 0, byrow = TRUE, ncol = numberOfSpecies)),
          1*(matrix(a[i, -1] > 0, byrow = TRUE, ncol = numberOfSpecies)),
          index.family = "jaccard"
        )
        data.frame(
          Time1 = a[i - 1, 1][[1]],
          Time2 = a[i, 1][[1]],
          Env = "Multi", 
          Measurement = names(outs),
          Value = unlist(outs, use.names = FALSE),
          stringsAsFactors = FALSE
        )
      }, a = loaded$Abundance))
  ))
}

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

# Prep Load: ##################################################################
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

# Run Calculations: ###########################################################
betapartition <- foreach::foreach(
  x = iterators::iter(
    dir(datfolders, full.names = TRUE,
        pattern = "(Simulation|Result|Intervention)")
  ), .packages = c("dplyr", "RMTRCode2", "betapart")
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
               "Betapart",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_Betapart_",
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
    
    # Make sure spurious zeros from the calculation aren't preserved.
    toEliminate <- 
      loaded$Abundance[, -1] < loaded$Parameters$EliminationThreshold
    loaded$Abundance[, -1][toEliminate] <- 0
    
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
              
              # divs <- RMTRCode2::thinAndCalculateDiversities(#Calculate_Diversity(
              #   loaded_subset,
              #   nspecies = c(Basal = sum(.x$Type == "Basal"),
              #                Consumer = sum(.x$Type == "Consumer")),
              #   # My standard approach for nspecies.
              #   preferred_rows_per_event = rows_per_event,
              #   divide_time_by = 1
              # )
              divs <- runBetapart(loaded_subset, nrow(.x))
              
              # thinAndCalculateDiversities yields lists
              divs$NicheValues <- .y
              return(divs)
            }
          )
        }
      )
      
      names(Diversity) <- Niches
      
      Diversity$Diversities <- runBetapart(loaded, numberOfSpecies)
      
      Diversity$Ellipsis <- loaded$Ellipsis
    } else {
      Diversity <- list(
        Diversities = runBetapart(loaded, numberOfSpecies),
        Ellipsis = loaded$Ellipsis
      )
    }
    
    save(Diversity, file = filename)
  }
  
  return(Diversity)
}


# Cleanup: ####################################################################
if (exists("clust")) {
  parallel::stopCluster(clust)
}
