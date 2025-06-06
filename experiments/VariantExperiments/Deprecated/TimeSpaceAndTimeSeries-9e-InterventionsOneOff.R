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
dirDate <- "2025-04-15"
baseTag <-  "TSTS_Simulation" # Note the distinction ("s").

seeds <- withRandom({runif(8)*1e8}, seed = 31300093)

parameterChoices <- data.frame(
  poolpatchDictionaryChoice = c(142486, 142486), # regular pools
  poolpatchSeedChoice = seeds[1], # shared pool
  dynamicsDictionaryChoice = c(4929, 4929), # regular dynamics
  dynamicsSeedChoice = seeds[2], # shared dynamics
  eventsDictionaryChoice = c(28, 28), # regular events
  eventsSeedChoice = seeds[3], # shared events
  initialConditionsDictionaryChoice = c(1, 1), # empty habitat
  initialConditionsSeedChoice = seeds[4], # shared ICs
  dispersalDictionaryChoice = c(NA, NA), # no dispersal
  distanceDictionaryChoice = c(1, 1), # Start with rho.noop
  affinityDictionaryChoice = c(20, 21), # runif then evensplit, start at 0.5s
  affinitySeedChoice = seeds[5:6], # different affinities
  interventionPatchDictionaryChoice = c(112, 114), # ->0.25, ->0.75
  interventionPatchSeedChoice = seeds[7], # consistency
  interventionTimeDictionaryChoice = c(1, 1), # mean method, no randomness
  interventionTimeSeedChoice = seeds[8], # consistency
  interventionDispersalDictionaryChoice =  c(NA, NA), # no dispersal,
  interventionDistanceDictionaryChoice = c(5, 5) # 5 fold change both directions
)

toExport <- unlist(lapply(
  1:9, # Non-seed, non-previous columns of parameterChoices
  function(id, pcs, dicts) {
    indices <- unique(unlist(pcs[, id]))
    dictChoices <- dicts[[id]][indices, ]
    stringCols <- unlist(lapply(dictChoices,
                                function(x) all(is.character(x))))
    dictChoices[, stringCols] %>% unlist()
  },
  pcs = parameterChoices %>% dplyr::select(-dplyr::ends_with("Seed")),
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
) %dopar% {
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
