library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

directory <- '.'
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Interventions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9b-SimulationFunction.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9d-InterventionFunction.R"))

dirTag <- "TSTS_Simulations"
dirDate <- "2025-05-12"
baseTag <-  "TSTS_Simulation" # Note the distinction ("s").

cores <- 3

seeds <- withRandom({runif(10)*1e8}, seed = 94311482)

parameterChoices <- data.frame(
  poolpatchDictionaryChoice = rep(142486, 3), # regular pools
  poolpatchSeedChoice = seeds[1], # shared pool
  dynamicsDictionaryChoice = rep(4929, 3), # regular dynamics
  dynamicsSeedChoice = seeds[2], # shared dynamics
  eventsDictionaryChoice = rep(28, 3), # regular events
  eventsSeedChoice = seeds[3], # shared events
  initialConditionsDictionaryChoice = rep(1, 3), # empty habitat
  initialConditionsSeedChoice = seeds[4], # shared ICs
  dispersalDictionaryChoice = rep(NA, 3), # no dispersal
  distanceDictionaryChoice = rep(5, 3), # Start with 5^(1 - 2*dist).
  affinityDictionaryChoice = c(1, 15, 29), # (0), (0), (0.5), (1) with 100% 0.
  affinitySeedChoice = seeds[5:7], # different affinities
  interventionPatchDictionaryChoice = rep(111, 3), # ->0
  interventionPatchSeedChoice = seeds[8], # consistency
  interventionTimeDictionaryChoice = rep(3, 3), # 5% in.
  interventionTimeSeedChoice = seeds[9], # consistency
  interventionDispersalDictionaryChoice =  rep(NA, 3), # no dispersal,
  interventionDistanceDictionaryChoice = rep(5, 3) # 5 fold change both directions
)

eventsDictionaryOrigin[28,]$EventsNumberMultiplier <- 200

if (cores > 1) {
  clust <- parallel::makeCluster(min(nrow(parameterChoices), cores))
  doParallel::registerDoParallel(clust)
  `%op%` <- `%dopar%`
} else {
  `%op%` <- `%do%`
}

toExport <- unlist(lapply(
  1:9, # Non-seed, non-previous columns of parameterChoices
  function(id, pcs, dicts) {
    indices <- unique(unlist(pcs[, id]))
    dictChoices <- dicts[[id]][indices, ]
    stringCols <- unlist(lapply(dictChoices,
                                function(x) all(is.character(x))))
    dictChoices[, stringCols] %>% unlist()
  },
  pcs = parameterChoices %>% dplyr::select(-dplyr::contains("Seed")),
  dicts = list(poolpatchDictionaryOrigin, dynamicsDictionaryOrigin,
               eventsDictionaryOrigin, initialConditionsDictionaryOrigin,
               dispersalDictionaryOrigin,
               distanceDictionaryOrigin, affinityDictionaryOrigin,
               interventionPatchDictionaryOrigin,
               interventionTimeDictionaryOrigin)
)) %>% unique()
toExport <- toExport[!grepl("=", toExport, fixed = TRUE) &
                       !grepl("::", toExport, fixed = TRUE) &
                       !is.na(toExport) ]
toExport <- toExport[toExport %in% ls()]

success <- foreach::foreach(
  pc = iterators::iter(parameterChoices, by = "row"),
  .packages = c("RMTRCode2", "dplyr"), .export = toExport
) %op% {
  # ) %do% {
  pc <- as.list(pc) # untibble so we are passing numerics and strings.
  fileID <- file.path(
    paste0(
      dirTag, "_", pc[[1]], "-", pc[[3]],
      "_", pc[[2]], "-", pc[[4]], "_", dirDate
    ),
    paste0(
      baseTag,
      "_", pc[[1]], "-", pc[[3]], "-", pc[[5]],
      "-", pc[[7]], "-", pc[[9]], "-", pc[[10]], "-", pc[[11]],
      "_", pc[[2]], "-", pc[[4]], "-", pc[[6]], "-", pc[[8]], "-", pc[[12]],
      ".RData"
    )
  )
  if (file.exists(fileID))
    interventionWrapper(
      # ID = list(
      #   Tag = fileTags,
      #   poolpatchDictionaryChoice = pc$pp, [1]
      #   poolpatchSeedChoice = pc$ppSeed, [2]
      #   dynamicsDictionaryChoice = pc$dyn, [3]
      #   dynamicsSeedChoice = pc$dynSeed, [4]
      #   eventsDictionaryChoice = pc$events, [5]
      #   eventsSeedChoice = pc$eventsSeed, [6]
      #   initialConditionsDictionaryChoice = pc$initconds, [7]
      #   initialConditionsSeedChoice = pc$initcondsSeed, [8]
      #   dispersalDictionaryChoice = pc$dispersal, [9]
      #   distanceDictionaryChoice = pc$distance, [10]
      #   affinityDictionaryChoice = pc$affinity, [11]
      #   affinitySeedChoice = pc$affinitySeed, [12]
      #   Date = fileDates
      # ),
      ID = fileID,
      interventionPatchDictionaryChoice = pc[[13]],
      interventionPatchSeedChoice = pc[[14]],
      interventionTimeDictionaryChoice = pc[[15]],
      interventionTimeSeedChoice = pc[[16]],
      interventionDispersalDictionaryChoice = pc[[17]],
      interventionDistanceDictionaryChoice = pc[[18]],
      returnResults = FALSE,
      saveResults = TRUE,
      skipIfSaveExists = TRUE
    )
}

if (cores > 1)
  parallel::stopCluster(clust)
