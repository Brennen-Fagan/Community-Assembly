# As in 9c-Simulations, but running "4" total specific example simulations.
# These are meant to look for if systems ending in land-use (0) end in the
# same species community after a "long enough" time period.
# So we're running (0), (0) -> (0), (0.5) -> (0), and (1) -> (0).
# After previous experience, we expect 1) and 2) to be same throughout, while
# 3) and 4) are different even after coming to the same land-use.
# NOTE THE CHANGE TO EVENTSDICTIONARY TO INCREASE THE LENGTH OF THE SIMULATION.
logisticCarryingCapacity <- NULL

directory <- '.'
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9-Dictionaries.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9b-SimulationFunction.R"))

library(parallel)
library(doParallel)
library(foreach)
library(iterators)

cores <- 3

# > runif(1)*1e8
# [1] 94311482
seeds <- withRandom({runif(7)*1e8}, seed = 94311482)

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
  affinitySeedChoice = seeds[5:7] # different affinities
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
  1:7, # Non-seed columns of parameterChoices
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
               distanceDictionaryOrigin, affinityDictionaryOrigin)
)) %>% unique()
toExport <- toExport[!grepl("=", toExport, fixed = TRUE) &
                       !grepl("::", toExport, fixed = TRUE) &
                       !is.na(toExport) ]
toExport <- toExport[toExport %in% ls()]

success <- foreach::foreach(
  pc = iterators::iter(parameterChoices, by = "row"),
  .packages = c("RMTRCode2", "dplyr"), .export = toExport
) %op% {
  pc <- unlist(pc) # untibble so we are passing numerics.
  simulationWrapper(
    poolpatchDictionaryChoice = pc[1],
    poolpatchSeedChoice = pc[2],
    dynamicsDictionaryChoice = pc[3],
    dynamicsSeedChoice = pc[4],
    eventsDictionaryChoice = pc[5],
    eventsSeedChoice = pc[6],
    initialConditionsDictionaryChoice = pc[7],
    initialConditionsSeedChoice = pc[8],
    dispersalDictionaryChoice = pc[9],
    distanceDictionaryChoice = pc[10],
    affinityDictionaryChoice = pc[11],
    affinitySeedChoice = pc[12],
    logisticCarryingCapacity = logisticCarryingCapacity,
    returnResults = FALSE,
    saveResults = TRUE,
    skipIfSaveExists = TRUE
  )
}
