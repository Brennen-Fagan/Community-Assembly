# runSimulations <- TRUE
runSimulationsFlag <- # Default to TRUE if runSimulations does not exist.
  (!exists("runSimulations") || get("runSimulations") == TRUE)

simulationsTarget <- switch(
  EXPR = {if(exists("simulationsTargetIndex")) simulationsTargetIndex else "NA"},
  "1" = "TimeSpaceAndTimeSeries-8a-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-8a-DictionaryIDs.R" # Default
)
directory <- '.'
source(file.path(directory, simulationsTarget))

if (runSimulationsFlag) {
  source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
  source(file.path(directory, "TimeSpaceAndTimeSeries-0-Dictionaries.R"))
  source(file.path(directory, "TimeSpaceAndTimeSeries-8b-SimulationFunction.R"))

  library(parallel)
  library(doParallel)
  library(foreach)
  library(iterators)
}

# Setup Notes: ################################################################
# Seeds already used through 7:
#   1:29 (poolpatch), 1:29 (dynamics), 1:52 (events)

# Parameters: #################################################################
numberPools <- 4
numberHistories <- 1
numberAffinities <- 1
poolSeedOffset <- 30
historiesSeedOffset <- 53
affinitiesSeedOffset <- 1

parameterChoices <- with(
  experiments,
  expand.grid( # All combinations!
    pp = ppDO$ID,
    dyn = dynDO$ID,
    events = eDO$ID,
    dispersal = dispDO$ID,
    affinity = aDO$ID,
    distance = distDO$ID
  )
)

# Add seeds
# 4 unique pool-patch replicates, each pool assigned (one set of) one patch.
# Additionally, each pool-patch combination is assigned one history.
# Similarly, each pool-patch-affinity combination is assigned one affinity seed.

parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp
) %>% dplyr::mutate(
  ppSeed = (dplyr::cur_group_id() - 1) * numberPools + poolSeedOffset
) %>% dplyr::ungroup(
)

parameterChoices <- dplyr::bind_rows(
  parameterChoices,
  if (numberPools > 1) {
    dplyr::bind_rows(lapply(1:(numberPools - 1), function(i) {
      parameterChoices %>% dplyr::mutate(
        ppSeed = ppSeed + i
      )
    }))
  }
) %>% dplyr::mutate(
  dynSeed = ppSeed
)

parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp, ppSeed # Maybe need to add events here?
) %>% dplyr::mutate(
  # Repeats across distances and dispersals
  eventsSeed =
    (dplyr::cur_group_id() - 1) * numberHistories + historiesSeedOffset
) %>% dplyr::ungroup(
)

parameterChoices <- dplyr::bind_rows(
  parameterChoices,
  if (numberHistories > 1) {
    dplyr::bind_rows(lapply(1:(numberHistories - 1), function(i) {
      parameterChoices %>% dplyr::mutate(
        eventsSeed = eventsSeed + i
      )
    }))
  }
)

parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp, ppSeed, affinity # Maybe need to add events here?
) %>% dplyr::mutate(
  # Repeats across distances and dispersals
  affinitySeed =
    (dplyr::cur_group_id() - 1) * numberAffinities + affinitiesSeedOffset
) %>% dplyr::ungroup(
)

parameterChoices <- dplyr::bind_rows(
  parameterChoices,
  if (numberAffinities > 1) {
    dplyr::bind_rows(lapply(1:(numberAffinities - 1), function(i) {
      parameterChoices %>% dplyr::mutate(
        affinitySeed = affinitySeed + i
      )
    }))
  }
)

# Order matters due to some transforms for the wrapper.
# poolpatchDictionaryChoice = pc$pp, [1]
# poolpatchSeedChoice = pc$ppSeed, [2]
# dynamicsDictionaryChoice = pc$dyn, [3]
# dynamicsSeedChoice = pc$dynSeed, [4]
# eventsDictionaryChoice = pc$events, [5]
# eventsSeedChoice = pc$eventsSeed, [6]
# dispersalDictionaryChoice = pc$dispersal, [7]
# distanceDictionaryChoice = pc$distance, [8]

parameterChoices <- parameterChoices %>% dplyr::select(
  dplyr::starts_with("pp"),
  dplyr::starts_with("dyn"),
  dplyr::starts_with("events"),
  dplyr::starts_with("dispersal"),
  dplyr::starts_with("distance"),
  dplyr::starts_with("affinity"),
  dplyr::everything()
) %>% dplyr::arrange(
  dplyr::across(dplyr::everything())
)

# Run across each row of parameterChoices: ####################################
if (runSimulationsFlag) {
  clust <- parallel::makeCluster(12)
  doParallel::registerDoParallel(clust)
  toExport <- unlist(lapply(
    1:6, # Non-seed columns of parameterChoices
    function(id, pcs, dicts) {
      indices <- unique(unlist(pcs[, id]))
      dictChoices <- dicts[[id]][indices, ]
      stringCols <- unlist(lapply(dictChoices,
                                  function(x) all(is.character(x))))
      dictChoices[, stringCols] %>% unlist()
    },
    pcs = parameterChoices %>% dplyr::select(-dplyr::ends_with("Seed")),
    dicts = list(poolpatchDictionaryOrigin, dynamicsDictionaryOrigin,
                 eventsDictionaryOrigin, dispersalDictionaryOrigin,
                 distanceDictionaryOrigin, affinityDictionaryOrigin)
  )) %>% unique()
  toExport <- toExport[!grepl("=", toExport, fixed = TRUE) &
                         !grepl("::", toExport, fixed = TRUE) &
                         !is.na(toExport) ]
  toExport <- toExport[toExport %in% ls()]

  success <- foreach::foreach(
    pc = iterators::iter(parameterChoices, by = "row"),
    .packages = c("RMTRCode2", "dplyr"), .export = toExport
  ) %dopar% {
    pc <- unlist(pc) # untibble so we are passing numerics.
    simulationWrapper(
      poolpatchDictionaryChoice = pc[1],
      poolpatchSeedChoice = pc[2],
      dynamicsDictionaryChoice = pc[3],
      dynamicsSeedChoice = pc[4],
      eventsDictionaryChoice = pc[5],
      eventsSeedChoice = pc[6],
      dispersalDictionaryChoice = pc[7],
      distanceDictionaryChoice = pc[8],
      affinityDictionaryChoice = pc[9],
      affinitySeedChoice = pc[10],
      returnResults = FALSE,
      saveResults = TRUE,
      skipIfSaveExists = TRUE
    )
  }

  parallel::stopCluster(clust)
}
