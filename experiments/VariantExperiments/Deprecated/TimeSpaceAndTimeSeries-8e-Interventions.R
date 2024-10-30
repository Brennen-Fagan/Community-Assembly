library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

directory <- '.'
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Interventions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-8b-SimulationFunction.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-8d-InterventionFunction.R"))

dirTag <- "TSTS_Simulations"
dirDate <- "2024-10-22"
baseTag <-  "TSTS_Simulation" # Note the distinction

runSimulations <- FALSE
source(file.path(directory, "TimeSpaceAndTimeSeries-8c-Simulations.R"))
# For Side Effects, specifically parameterChoices.

stopifnot(exists("parameterChoices"))

# All combinations!!!
parameterChoices <- parameterChoices %>% dplyr::full_join(
  with(experiments, expand.grid(
    ip = iPDO$ID,
    ipSeed = 1, # No random components
    it = iTDO$ID,
    itSeed = 1, # No random components
    itDisp = iDispChoice,
    itDist = iDistChoice,
    stringsAsFactors = FALSE
  )),
  by = character() # cross-join
)

# Run across each row of parameterChoices: ####################################
clust <- parallel::makeCluster(12)
doParallel::registerDoParallel(clust)

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
    ID = file.path(
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
    ),
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

parallel::stopCluster(clust)
