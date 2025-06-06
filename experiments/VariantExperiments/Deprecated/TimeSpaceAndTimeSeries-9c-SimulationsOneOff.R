# As in 9c-Simulations, but we're only running two specific example simulations.
# These are based on 9a9, but are meant to be archetypal.
logisticCarryingCapacity <- NULL

directory <- '.'
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9-Dictionaries.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9b-SimulationFunction.R"))

library(parallel)
library(doParallel)
library(foreach)
library(iterators)

cores <- 1

# > runif(1)*1e8
# [1] 31300093
seeds <- withRandom({runif(6)*1e8}, seed = 31300093)

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
  affinitySeedChoice = seeds[5:6] # different affinities
)

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
  pcs = parameterChoices %>% dplyr::select(-dplyr::ends_with("Seed")),
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
