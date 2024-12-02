# runSimulations <- TRUE
runSimulationsFlag <- # Default to TRUE if runSimulations does not exist.
  (!exists("runSimulations") || get("runSimulations") == TRUE)

simulationsTargets <- c(
  "TimeSpaceAndTimeSeries-9a-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-9a2-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-9a3-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-9a4-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-9a5-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-9a6-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-9a7-DictionaryIDs.R",
  "TimeSpaceAndTimeSeries-9a8-DictionaryIDs.R"#,
)

simulationsTarget <- switch(
  EXPR = {
    if(exists("simulationsTargetIndex")) simulationsTargetIndex else "NA"
  },
  "1" = simulationsTargets[1],
  "2" = simulationsTargets[2],
  "3" = simulationsTargets[3],
  "4" = simulationsTargets[4],
  "5" = simulationsTargets[5],
  "6" = simulationsTargets[6],
  "7" = simulationsTargets[7],
  "8" = simulationsTargets[8],
  simulationsTargets[length(simulationsTargets)] # Default
)
directory <- '.'
source(file.path(directory, simulationsTarget))

if (runSimulationsFlag) {
  source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
  source(file.path(directory, "TimeSpaceAndTimeSeries-9-Dictionaries.R"))
  source(file.path(directory, "TimeSpaceAndTimeSeries-9b-SimulationFunction.R"))

  library(parallel)
  library(doParallel)
  library(foreach)
  library(iterators)
}

cores <- 8 # Normally happy to put higher, but there's not the redundancy of
# pools in this go round, so each pool needs to be made, and that takes a lot
# of cores!

# Setup Notes: ################################################################
# Seeds already used through 8:
#   1:34 (poolpatch), 1:34 (dynamics), 1:57 (events)
# Lost track somewhere of initialConditionsSeedOffset and affinitiesSeedOffset,
# so setting to max of these two values: 81 (+ 1).
# through 9a: 84-84-107-131-131
# Parameters: #################################################################
if (simulationsTarget == simulationsTargets[1]) {
  numberPools <- 1
  numberHistories <- 1
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 35
  historiesSeedOffset <- 58
  initialConditionsSeedOffset <- 82
  affinitiesSeedOffset <- 82
} else if (simulationsTarget == simulationsTargets[2]) {
  numberPools <- 4
  numberHistories <- 4
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 85
  historiesSeedOffset <- 108
  initialConditionsSeedOffset <- 132
  affinitiesSeedOffset <- 132
} else if (simulationsTarget == simulationsTargets[3]) {
  numberPools <- 1
  numberHistories <- 1
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 89
  historiesSeedOffset <- 123
  initialConditionsSeedOffset <- 135
  affinitiesSeedOffset <- 135
} else if (simulationsTarget == simulationsTargets[4]) {
  numberPools <- 3
  numberHistories <- 4
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 188
  historiesSeedOffset <- 222
  initialConditionsSeedOffset <- 234
  affinitiesSeedOffset <- 234
} else if (simulationsTarget == simulationsTargets[5]) {
  numberPools <- 1
  numberHistories <- 1
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 191
  historiesSeedOffset <- 234
  initialConditionsSeedOffset <- 237
  affinitiesSeedOffset <- 237
} else if (simulationsTarget == simulationsTargets[6]) {
  numberPools <- 1
  numberHistories <- 1
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 192
  historiesSeedOffset <- 235
  initialConditionsSeedOffset <- 238
  affinitiesSeedOffset <- 242
} else if (simulationsTarget == simulationsTargets[7]) {
  numberPools <- 1
  numberHistories <- 1
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 242
  historiesSeedOffset <- 285
  initialConditionsSeedOffset <- 288
  affinitiesSeedOffset <- 292
} else if (simulationsTarget == simulationsTargets[8]) {
  numberPools <- 2
  numberHistories <- 1
  numberAffinities <- 1
  numberInitConds <- 1
  poolSeedOffset <- 341
  historiesSeedOffset <- 384
  initialConditionsSeedOffset <- 387
  affinitiesSeedOffset <- 391
} else {
  stop("Settings not detected properly.")
}


parameterChoices <- dplyr::bind_rows(lapply(
  experiments, function(experiment) {
    with(
      experiment,
      expand.grid( # All combinations!
        pp = ppDO$ID,
        dyn = dynDO$ID,
        events = eDO$ID,
        initconds = icDO$ID,
        dispersal = dispDO$ID,
        affinity = aDO$ID,
        distance = distDO$ID
      )
    )})) %>% dplyr::arrange(
      pp, dyn, events, initconds, dispersal, affinity, distance
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

parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp, ppSeed, initconds # Maybe need to add events here?
) %>% dplyr::mutate(
  # Repeats across distances and dispersals
  initcondsSeed =
    (dplyr::cur_group_id() - 1) * numberInitConds + initialConditionsSeedOffset
) %>% dplyr::ungroup(
)

parameterChoices <- dplyr::bind_rows(
  parameterChoices,
  if (numberInitConds > 1) {
    dplyr::bind_rows(lapply(1:(numberInitConds - 1), function(i) {
      parameterChoices %>% dplyr::mutate(
        initcondsSeed = initcondsSeed + i
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
# initialConditionsDictionaryChoice = pc$initconds, [7]
# initialConditionsSeedChoice = pc$initcondsSeed, [8]
# dispersalDictionaryChoice = pc$dispersal, [9]
# distanceDictionaryChoice = pc$distance, [10]
# affinityDictionaryChoice = pc$affinity, [11]
# affinitySeedChoice = pc$affinitySeed, [12]

parameterChoices <- parameterChoices %>% dplyr::select(
  dplyr::starts_with("pp"),
  dplyr::starts_with("dyn"),
  dplyr::starts_with("events"),
  dplyr::starts_with("initconds"),
  dplyr::starts_with("dispersal"),
  dplyr::starts_with("distance"),
  dplyr::starts_with("affinity"),
  dplyr::everything()
) %>% dplyr::arrange(
  dplyr::across(dplyr::everything())
)

# Run across each row of parameterChoices: ####################################
if (runSimulationsFlag) {
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
      returnResults = FALSE,
      saveResults = TRUE,
      skipIfSaveExists = TRUE
    )
  }

  if (cores > 1)
    parallel::stopCluster(clust)
}
