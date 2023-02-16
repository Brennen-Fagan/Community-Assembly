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

outputLocation <- file.path(".", paste0("SaveDiversity_",
                                        "2022-06-17" #Sys.Date()
                                        ))
if (!dir.exists(outputLocation)) {
  dir.create(outputLocation, showWarnings = FALSE)
}

# Libraries: ##################################################################
print("Loading libraries")
librarypath <- file.path(".", "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

packages <- c(
  "dplyr",      # Data Manipulation
  "Matrix",     # Common Format
  "readr",
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

# # Receive the Selection: ######################################################
print("Receiving number of cores:")

cargs <- 4
# cargs <- as.numeric(commandArgs(trailingOnly = TRUE))

clust <- parallel::makeCluster(cargs[1], outfile = "")
doParallel::registerDoParallel(clust)

# Functions: ##################################################################

# Recycling from SecondAttempt-Doc-Analysis2-Gallery.Rmd.
thinAndCalculateDiversities <- function(loaded, nspecies) {
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
  print("prep")

  ### Alpha Diversity: ##################################################
  diversity_alpha <- lapply(
    1:loaded$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      env_basal <- env[, 1:nspecies[1]]
      env_consumer <- env[, nspecies[1] + 1:nspecies[2]]
      richness <- rowSums(env)
      richness_basal <- rowSums(env_basal)
      richness_consumer <- rowSums(env_consumer)
      species <- apply(
        env, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_basal <- apply(
        env_basal, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_consumer <- apply(
        env_consumer, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0) + nspecies[1]) # + basals doesn't change it, but
          # results in a more intuitive numbering system of species.
        }
      )
      data.frame(Time = time,
                 Richness = richness,
                 Richness_Basal = richness_basal,
                 Richness_Consumer = richness_consumer,
                 Species = species,
                 Species_Basal = species_basal,
                 Species_Consumer = species_consumer,
                 Environment = i,
                 stringsAsFactors = FALSE)
    },
    abund = loaded$Abundance,
    numSpecies = sum(nspecies)
  )

  diversity_alpha <- dplyr::bind_rows(diversity_alpha)

  print("alpha")
  ### Gamma Diversity: ##################################################
  diversity_gamma <- diversity_alpha %>% dplyr::group_by(
    Time
  ) %>% dplyr::summarise(
    Mean = mean(Richness),
    Mean_Basal = mean(Richness_Basal),
    Mean_Consumer = mean(Richness_Consumer),
    Var = var(Richness),
    Var_Basal = var(Richness_Basal),
    Var_Consumer = var(Richness_Consumer),

    SpeciesTotal = toString(sort(unique(unlist(strsplit(paste(
      Species, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Basal = toString(sort(unique(unlist(strsplit(paste(
      Species_Basal, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Consumer = toString(sort(unique(unlist(strsplit(paste(
      Species_Consumer, collapse = ", "), split = ", ", fixed = TRUE))))),

    Gamma = unlist(lapply(
      strsplit(
        SpeciesTotal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Gamma_Basal = unlist(lapply(
      strsplit(
        SpeciesTotal_Basal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Gamma_Consumer = unlist(lapply(
      strsplit(
        SpeciesTotal_Consumer, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    ))
  ) %>% tidyr::pivot_longer(
    cols = c(Mean, Var, Gamma),
    names_to = "Aggregation",
    values_to = "Richness"
  ) %>% dplyr::mutate(
    Basals = case_when(
      Aggregation == "Mean" ~ Mean_Basal,
      Aggregation == "Var" ~ Var_Basal,
      Aggregation == "Gamma" ~ as.numeric(Gamma_Basal),
      TRUE ~ -1
    ),
    Consumers = case_when(
      Aggregation == "Mean" ~ Mean_Consumer,
      Aggregation == "Var" ~ Var_Consumer,
      Aggregation == "Gamma" ~ as.numeric(Gamma_Consumer),
      TRUE ~ -1
    )
  ) %>% dplyr::select(
    -Mean_Basal, -Var_Basal, -Gamma_Basal,
    -Mean_Consumer, -Var_Consumer, -Gamma_Consumer,
    -SpeciesTotal, -SpeciesTotal_Basal, -SpeciesTotal_Consumer
  )

  print("gamma")
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
  print("beta")

  ### Return Diversities: ###############################################
  return(list(
    alpha = diversity_alpha %>% dplyr::select(
      -Species_Basal, -Species_Consumer
    ),
    beta = diversity_beta,
    gamma = diversity_gamma
  ))
}

# Files: ######################################################################
print("Identifying files.")

directories <- dir(path = ".", pattern = "SaveOutput[_]",
                   full.names = TRUE, include.dirs = TRUE)
files <- dir(path = directories,
             pattern = "^MNA[-]Master.+[.]RData$",
             full.names = TRUE, recursive = TRUE,
             include.dirs = TRUE, no.. = TRUE)

# Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
filesParameters <- dir(path = "Prepared_2021-12-20", pattern = "^MNA[-].+Cases[.]csv$",
                       full.names = TRUE, recursive = TRUE,
                       include.dirs = TRUE, no.. = TRUE)
parametersMaster <- load(file.path("Prepared_2021-12-20", "MNA-MasterParameters.RData"))
parametersCSVs <- lapply(filesParameters, utils::read.csv)

# Note the problem observed 31/01/2022, 1400 in the calendar.
# To accomodate, reassign inf to the *end* of the Space list.
systemMods$SpaceDistanceMultiplier <- systemMods$SpaceDistanceMultiplier[c(2:9, 1)]

# 0 In Parallel: ##############################################################
results <- foreach::foreach(
  file = iterators::iter(files), .packages = "dplyr"
) %dopar% {
  print(file)
  # Retrieve the trailing id numbers from before the file extension.
  # 1st: Set, 2nd: CaseNumber, 3rd: History, 4th: Part
  # Note, if last two are not present, all histories are bundled together.
  # If last two are present, each file has a single part of a single history.
  idNums <- suppressWarnings(na.omit(
    as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                             split = "-", fixed = TRUE)[[1]],
                    n = 4))),
    classes = "simpleWarning")

  print(idNums)

  fileName <- file.path(
    outputLocation,
    paste0("MNA-Master", LETTERS[idNums[1]], "-Cases-DiversityBC-",
           paste(idNums, collapse = "-"), ".RData")
  )

  if (file.exists(fileName)) {
    print("Already done.")

    if (length(idNums) == 2){
      return(c(idNums, NA, NA))
    } else {
      return(idNums)
    }
  }

  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    fileContents <- get(fileContents)

    print("loaded")

    # 2 Extract Alpha, Beta, Gamma diversities over time: ######################
    if (length(idNums) == 2) { # fileContents is a list of results
      Diversity <- list(
        Diversities = lapply(
          fileContents,
          FUN = thinAndCalculateDiversities,
          nspecies = systemBase$species
        ))
    } else if (length(idNums) == 4) {# fileContents is a result
      Diversity <- list(
        Diversities = list(thinAndCalculateDiversities(
          fileContents,
          nspecies = systemBase$species
        ))
      )
    } else stop("idNums is not 2 or 4.")
    print("diversity")
    # 6. Save the biodiversity values...: #####################################
    # separately, but alongside its parameters.
    # idNums indicates which variation we are dealing with, as well as the
    # associated parameters.
    # This is split between MNA-MasterParameters.RData and the .csv files.
    theParameters <- parametersCSVs[[idNums[1]]][((idNums[2] - 1) %/% 10) + 1,]
    Diversity$PoolMod <- systemMods$PoolBodySizes[theParameters$Framework, ]
    Diversity$NoiseMod <- systemMods$InteractionsNoiseMultiplier[theParameters$Framework]
    Diversity$NeutralMod<- systemMods$NeutralRateMultipliers[theParameters$Neutral,]
    Diversity$SpaceMod <- systemMods$SpaceDistanceMultiplier[theParameters$Space]

    save(Diversity,
         file = fileName
    )
    print("save")

    # Ran into memory issue. We'll see if this fixes it.
    rm(Diversity)
    rm(fileContents)
    gc()

    ###   3. Record Data as Present: ##########################################
    # This is already in the file list, but we can make it easier by taking the
    # id numbers. If there is an error at any point, record the relevant id
    # numbers and the error to know that something went wrong.
    if (length(idNums) == 2){
      return(c(idNums, NA, NA))
    } else {
      return(idNums)
    }
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
names(results_present) <- c("Set", "CaseNumber", "History", "Part")
results_present <- results_present %>% dplyr::arrange(Set, CaseNumber)
results_missing <- data.frame(
  Set = rep(1:3, unlist(lapply(parametersCSVs[1:3], nrow)) * 10),
  CaseNumber = unlist(lapply(
    unlist(lapply(parametersCSVs[1:3], nrow)) * 10,
    seq, from = 1, by = 1
  ))
)
results_missing <- dplyr::anti_join(
  results_missing, results_present, by = c("Set", "CaseNumber")
)

# 8. Save list of missing data to a file: #####################################
write.csv(results_missing,
          file = file.path(outputLocation,
                           paste0("MNA-MasterMissingRuns.csv"))
)
write.csv(results_error,
          file = file.path(outputLocation,
                           paste0("MNA-MasterErrorRuns.csv"))
)
