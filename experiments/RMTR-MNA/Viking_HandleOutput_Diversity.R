# Script Goals: ###############################################################
#   0. In Parallel,
#   1. Load the Data,
#   2. Extract Alpha, Beta, Gamma diversities over time,
#   3. Record Data as Present,
#   4. Determine if we can compress the Data further.
#   5. Save the Data in a more compressed format if so,
#   6. Save the biodiversity values separately, but alongside its parameters.
#   7. Determine Data that is missing,
#   8. Save list of missing data to a file.

# Parameters: #################################################################
divide_time_by <- 1E4
preferred_rows_per_event <- 1.5
directory <- "."

outputLocation <- file.path(directory, paste0("SaveDiversity_", Sys.Date()))
if (!dir.exists(outputLocation)) {
  dir.create(outputLocation, showWarnings = FALSE)
}

# Libraries: ##################################################################
print("Loading libraries")
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

packages <- c(
  "dplyr",      # Data Manipulation
  "Matrix",     # Common Format
  "readr",      # CSV manipulation
  "vegan",      # Biodiversity measurements
  "parallel",   # Base parallel dependency
  "doParallel", # For compatible cluster
  "foreach",    # For parallel evaluation
  "iterators"   # For parallel evaluation
)

# Does not work because it tries to use the system-wide libraries... Oops.
# update.packages(repos = 'https://cloud.r-project.org',
#                 oldPkgs = packages, ask = FALSE)

for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, lib = librarypath,
                     repos = 'https://cloud.r-project.org',
                     dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

if (!require("RMTRCode2", character.only = TRUE)) {
  install.packages(
    "RMTRCode2_1.0.tar.gz", lib = librarypath,
    repos = NULL, type = "source"
  )
}
library(RMTRCode2)#, lib.loc = librarypath) # lib.loc shouldn't be necessary.

# Receive the Selection: ######################################################
print("Receiving number of cores:")

cargs <- as.numeric(commandArgs(trailingOnly = TRUE))

print(cargs)

clust <- parallel::makeCluster(cargs[1])
doParallel::registerDoParallel(clust)

# Files: ######################################################################
print("Identifying files.")

directories <- dir(path = ".", pattern = "SaveOutput[_]",
                   full.names = TRUE, include.dirs = TRUE)
files <- dir(path = directories,
             pattern = "^MNA[-].+[.]RData$",
             full.names = TRUE, recursive = TRUE,
             include.dirs = TRUE, no.. = TRUE)

# Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
filesParameters <- dir(path = ".", pattern = "^MNA[-].+Cases[.]csv$",
                       full.names = TRUE, recursive = TRUE,
                       include.dirs = TRUE, no.. = TRUE)
parametersMaster <- load("MNA-MasterParameters.RData")
parametersCSVs <- lapply(filesParameters, readr::read_csv)

# 0 In Parallel: ##############################################################
results <- foreach::foreach(
  file = iterators::iter(files), .packages = c("dplyr")
) %dopar% {
  # Retrieve the trailing two id numbers from before the file extension.
  idNums <- as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                                     split = "-", fixed = TRUE)[[1]],
                            n = 2))
  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    fileContents <- get(fileContents)

    # 2 Extract Alpha, Beta, Gamma diversities over time: ######################
    # Recycling from SecondAttempt-Doc-Analysis2-Gallery.Rmd.
    Diversity <- lapply(
      fileContents,
      FUN = function(loaded) {
        # We can't handle all of the data that we are going to be looking at;
        # a small sample had ~180k rows for ~5.3k events = ~34 rows per event.
        # To reduce it, we will divide time up so that there are about the
        # preferred number of rows per event.
        bythin <- floor((nrow(loaded$Abundance)
                         / nrow(loaded$Events))
                        / preferred_rows_per_event)

        loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                                 to = nrow(loaded$Abundance),
                                                 by = bythin), ]

        # Remove illegal values (that the numerical engine uses as inbetweens).
        toEliminate <- loaded$Abundance[, -1] <
          loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
        loaded$Abundance[, -1][toEliminate] <- 0
        loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divide_time_by

        # Convert to Binary, since we are only interested in richness.
        loaded$Abundance[, -1] <- loaded$Abundance[, -1] > 0

        # loaded$Abundance is now prepared.
        # Now, we calculate richnesses.

        ### Alpha Diversity: ##################################################
        diversity_alpha <- lapply(
          1:loaded$NumEnvironments,
          function(i, abund, numSpecies) {
            time <- abund[, 1]
            env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
            richness <- rowSums(env)
            species <- apply(
              env, MARGIN = 1,
              FUN = function(x) {
                toString(which(x > 0))
              }
            )
            data.frame(Time = time,
                       Richness = richness,
                       Species = species,
                       Environment = i,
                       stringsAsFactors = FALSE)
          },
          abund = loaded$Abundance,
          numSpecies = (ncol(loaded$Abundance) - 1) / loaded$NumEnvironments
        )

        diversity_alpha <- dplyr::bind_rows(diversity_alpha)

        ### Gamma Diversity: ##################################################
        diversity_gamma <- diversity_alpha %>% dplyr::group_by(
          Time
        ) %>% dplyr::summarise(
          Mean = mean(Richness),
          SpeciesTotal = toString(sort(unique(unlist(strsplit(paste(
            Species, collapse = ", "), split = ", ", fixed = TRUE))))),
          Gamma = unlist(lapply(
            strsplit(
              SpeciesTotal, split = ", ", fixed = TRUE
            ),
            function(x) length(x[x!=""])
          ))
        ) %>% tidyr::pivot_longer(
          cols = c(Mean, Gamma),
          names_to = "Aggregation",
          values_to = "Richness"
        ) %>% dplyr::select(-SpeciesTotal)

        ### Beta Diversity (Jaccard, Space): ##################################
        diversity_beta <- apply(
          loaded$Abundance,
          MARGIN = 1, # Rows
          function(row, envs) {
            time <- row[1]
            # Vegan complains about rows with all 0's.
            # The warning is generic, so we cannot silence it specifically.
            dists <- suppressWarnings(vegan::vegdist(
              method = "jaccard",
              x = matrix(row[-1] > 0, nrow = envs, byrow = TRUE)
            ))

            dataf <- expand.grid(
              Env1 = 1:envs,
              Env2 = 1:envs
            ) %>% dplyr::filter(
              Env1 < Env2
            ) %>% dplyr::mutate(
              Time = time,
              Jaccard = dists
            )

            return(dataf)
          },
          envs = loaded$NumEnvironments
        )

        ### Return Diversities: ###############################################
        return(list(
          alpha = diversity_alpha,
          beta = diversity_beta,
          gamma = diversity_gamma
        ))
      }
    )
    # 6. Save the biodiversity values...: #####################################
    # separately, but alongside its parameters.
    # idNums indicates which variation we are dealing with, as well as the
    # associated parameters.
    # This is split between MNA-MasterParameters.RData and the .csv files.
    theParameters <- parametersCSVs[[idNums[1]]][(idNums[2] %/% 10) + 1,]
    Diversity$PoolMod <- systemMods$PoolBodySizes[theParameters$Framework,]
    Diversity$NoiseMod <- systemMods$InteractionsNoiseMultiplier[theParameters$Framework]
    Diversity$NeutralMod<- systemMods$NeutralRateMultipliers[theParameters$Neutral,]
    Diversity$SpaceMod <- systemMods$SpaceDistanceMultiplier[theParameters$Space]

    save(Diversity,
         file = file.path(
           outputLocation,
           paste0("MNA-Master", LETTERS[idNums[1]], "-Cases-Diversity-",
                  idNums[1], "-", idNums[2], ".RData")
         )
    )

    # 4. Determine if we can compress the Data further. #######################
    # 5. Save the Data in a more compressed format if so, #####################
    # So I have an attempt at this running on Research 1, but the idea, below,
    # seems to run too slowly to be worthwhile.
    # Furthermore, in the first ~21k (/ ~179k), only ~6k were removed => 71%.
    # temp <- results[[1]]$Abundance
    # temp <- cbind(starttime = temp[, 1], temp)
    # current <- 1
    # for (r in 2:nrow(results[[1]]$Abundance)) {
    #   if (isTRUE(all.equal(temp[current, c(-1, -2)], results[[1]]$Abundance[r, -1]))) {
    #     temp[current, 2] <- results[[1]]$Abundance[r, 1]
    #     temp <- temp[-(current + 1), ]
    #   } else {
    #     current <- current + 1
    #   }
    # }

    ###   3. Record Data as Present: ##########################################
    # This is already in the file list, but we can make it easier by taking the
    # id numbers. If there is an error at any point, record the relevant id
    # numbers and the error to know that something went wrong.
    return(idNums)
  }, error = function(e, id = idNums) {
    return(list(ID = id,
                error = e))
  })
}

# Cleanup: ####################################################################
parallel::stopCluster(clust)

# 7. Determine Data that is missing: ##########################################
results_iflist <- unlist(lapply(results, is.list))
results_success <- do.call(rbind, results[!results_iflist])
results_error <- do.call(rbind,
                         lapply(results[results_iflist], function(x) x$ID))
results_present <- dplyr::bind_rows(
  data.frame(results_success),
  data.frame(results_error)
)
names(results_present) <- c("Case", "Simulation")
results_present <- results_present %>% dplyr::arrange(Case, Simulation)
results_missing <- data.frame(
  Case = rep(1:3, unlist(lapply(parametersCSVs, nrow)) * 10),
  Simulation = unlist(lapply(
    unlist(lapply(parametersCSVs, nrow)) * 10,
    seq, from = 1, by = 1
  ))
)
results_missing <- dplyr::anti_join(results_missing, results_present)

# 8. Save list of missing data to a file: #####################################
write.csv(results_missing,
          file = file.path(outputLocation,
                           paste0("MNA-MasterMissingRuns.csv"))
)
