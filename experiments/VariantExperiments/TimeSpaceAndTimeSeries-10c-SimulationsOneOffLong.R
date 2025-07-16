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
source(file.path(directory, "TimeSpaceAndTimeSeries-10-Dictionaries.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-10b-SimulationFunction.R"))

library(parallel)
library(doParallel)
library(foreach)
library(iterators)

cores <- 3

# > runif(1)*1e8
# [1] 8258686
# truncate so we don't have bad file names.
seeds <- floor(withRandom({runif(8)*1e8}, seed = 8258686))

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
  speciesAffinityDictionaryChoice = rep(1, 3), # 100% 0 species
  speciesAffinitySeedChoice = seeds[5], # shared species' preference seed
  patchAffinityDictionaryChoice = c(1, 3, 5), # (0), (0.5), (1) starting LUs.
  patchAffinitySeedChoice = seeds[6:8] # different patch land-uses
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
  1:8, # Non-seed columns of parameterChoices
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
               dispersalDictionaryOrigin, distanceDictionaryOrigin,
               speciesAffinityDictionaryOrigin, patchAffinityDictionaryOrigin)
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
    speciesAffinityDictionaryChoice = pc[11],
    speciesAffinitySeedChoice = pc[12],
    patchAffinityDictionaryChoice = pc[13],
    patchAffinitySeedChoice = pc[14],
    logisticCarryingCapacity = logisticCarryingCapacity,
    returnResults = FALSE,
    saveResults = TRUE,
    skipIfSaveExists = TRUE
  )
}
