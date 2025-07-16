# runSimulations <- TRUE
runSimulationsFlag <- # Default to TRUE if runSimulations does not exist.
  (!exists("runSimulations") || get("runSimulations") == TRUE)

logisticCarryingCapacity <-
  NULL # OR
# list(Basal = 1000) # 300 or 3000 also have arguments.

simulationsTargets <- c(
  "TimeSpaceAndTimeSeries-10a1-DictionaryIDs.R"#,
)

simulationsTarget <- switch(
  EXPR = {
    if(exists("simulationsTargetIndex")) simulationsTargetIndex else "NA"
  },
  "1" = simulationsTargets[1],
  # "2" = simulationsTargets[2],
  # "3" = simulationsTargets[3],
  # "4" = simulationsTargets[4],
  # "5" = simulationsTargets[5],
  # "6" = simulationsTargets[6],
  # "7" = simulationsTargets[7],
  # "8" = simulationsTargets[8],
  # "9" = simulationsTargets[9],
  simulationsTargets[length(simulationsTargets)] # Default
)
directory <- '.'
source(file.path(directory, simulationsTarget))

if (runSimulationsFlag) {
  source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
  source(file.path(directory, "TimeSpaceAndTimeSeries-10-Dictionaries.R"))
  source(file.path(directory, "TimeSpaceAndTimeSeries-10b-SimulationFunction.R"))

  library(parallel)
  library(doParallel)
  library(foreach)
  library(iterators)
}

if (!exists("cores"))
  cores <- 1 # Normally happy to put higher, but there's not the redundancy of
# pools in this go round, so each pool needs to be made, and that takes a lot
# of cores!

# Setup Notes: ################################################################
# Parameters: #################################################################
if (simulationsTarget == simulationsTargets[1]) {
    numberPools <- 44 # Power Analysis: detect 95% non-zero if 2 richdiff at 8sd.
    numberHistories <- 1 # Could rerun pools w/diff hists if needed, use offset.
    numberSAffinities <- 1 # Mostly for runif case, but should just do as above.
    numberPAffinities <- 1 # Essentially unimplemented at the moment.
    numberInitConds <- 1 # Starting from none anyways.
    # New seeds, so we can restart
    poolSeedOffset <- 1
    historiesSeedOffset <- 1
    dynamicsSeedOffset <- 1
    initialConditionsSeedOffset <- 1
    speciesAffinitiesSeedOffset <- 1
    patchAffinitiesSeedOffset <- 1
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
        saff = sADO$ID,
        paff = pADO$ID,
        distance = distDO$ID
      )
    )})) %>% dplyr::arrange(
      pp, dyn, events, initconds, dispersal, saff, paff, distance
    )

# Add seeds
# 4 unique pool-patch replicates, each pool assigned (one set of) one patch.
# Additionally, each pool-patch combination is assigned one history.
# Similarly, each pool-patch-affinity combination is assigned one affinity seed.

# Set up base seeds. pp and dyn are closely related.
parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp, dyn
) %>% dplyr::mutate(
  ppSeed = (dplyr::cur_group_id() - 1) * numberPools + poolSeedOffset,
  dynSeed = (dplyr::cur_group_id() - 1) * numberPools + dynamicsSeedOffset
) %>% dplyr::ungroup(
)

# Adjust for the number of replicates.
parameterChoices <- dplyr::bind_rows(
  parameterChoices,
  if (numberPools > 1) {
    dplyr::bind_rows(lapply(1:(numberPools - 1), function(i) {
      parameterChoices %>% dplyr::mutate(
        ppSeed = ppSeed + i,
        dynSeed = dynSeed + i
      )
    }))
  }
)

# For [1], eventsSeed will match ppSeed and dynSeed because we are only running
# one history per pool. This quickly changes if numberHistories ever changes.
parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp, ppSeed, dyn, dynSeed
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

# Next, need to consider the affinity structures decomposed.
# Each pool needs a seed assigned to it for each species affinity structure.
parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp, ppSeed, saff # Maybe need to add events here?
) %>% dplyr::mutate(
  # Repeats across distances and dispersals
  saffSeed =
    (dplyr::cur_group_id() - 1) * numberSAffinities + speciesAffinitiesSeedOffset
) %>% dplyr::ungroup(
)

parameterChoices <- dplyr::bind_rows(
  parameterChoices,
  if (numberSAffinities > 1) {
    dplyr::bind_rows(lapply(1:(numberSAffinities - 1), function(i) {
      parameterChoices %>% dplyr::mutate(
        saffSeed = saffSeed + i
      )
    }))
  }
)

# It is more ambiguous when to switch patch structures.
# Patch affinities could be treated as independent of pool, dynamics, history,
# etc. In practice though, we'd probably want to look for correlation structures
# between pool randomisations and patch randomisations.
# I'm going to opt for keeping it tied to pool as a result.
parameterChoices <- parameterChoices %>% dplyr::group_by(
  pp, ppSeed, paff # Maybe need to add events here?
) %>% dplyr::mutate(
  # Repeats across distances and dispersals
  paffSeed =
    (dplyr::cur_group_id() - 1) * numberPAffinities + patchAffinitiesSeedOffset
) %>% dplyr::ungroup(
)

parameterChoices <- dplyr::bind_rows(
  parameterChoices,
  if (numberPAffinities > 1) {
    dplyr::bind_rows(lapply(1:(numberPAffinities - 1), function(i) {
      parameterChoices %>% dplyr::mutate(
        paffSeed = paffSeed + i
      )
    }))
  }
)

# Initial conditions:
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
# speciesAffinityDictionaryChoice = pc$saff, [11]
# speciesAffinitySeedChoice = pc$saffSeed, [12]
# patchAffinityDictionaryChoice = pc$paff, [13]
# patchAffinitySeedChoice = pc$paffSeed, [14]

parameterChoices <- parameterChoices %>% dplyr::select(
  dplyr::starts_with("pp"),
  dplyr::starts_with("dyn"),
  dplyr::starts_with("events"),
  dplyr::starts_with("initconds"),
  dplyr::starts_with("dispersal"),
  dplyr::starts_with("distance"),
  dplyr::starts_with("saff"),
  dplyr::starts_with("paff"),
  dplyr::everything()
) %>% dplyr::arrange(
  dplyr::across(dplyr::everything())
)

# Going to try to get clever to spread out the pool generation.
parameterChoices <- parameterChoices %>% dplyr::mutate(
  firsts = !duplicated(interaction(ppSeed, dynSeed))
) %>% dplyr::group_by(
  firsts
) %>% dplyr::mutate(
  priority = ifelse(
    # If you're a first, increase priority so that there are two runs
    # between all cores and each new pool creation roughly.
    firsts, (dplyr::row_number() - 1) * cores * 2, dplyr::row_number()
  )
) %>% dplyr::arrange(
  priority
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

  if (cores > 1)
    parallel::stopCluster(clust)
}
