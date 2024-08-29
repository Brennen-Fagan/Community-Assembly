directory <- '.'
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-7a-SimulationFunction.R"))

library(parallel)
library(doParallel)
library(foreach)
library(iterators)

# Setup Notes: ################################################################
# Seeds already used from 6a:
#   1:10 (poolpatch), 1:10 (dynamics), 1:13 (events)

# Requested Combinations:
# poolpatchDictionaryOrigin %>% dplyr::filter(
#   NumberEnvironments %in% c(1, 2), NSpecies == 200,
#   SpeciesAffinities == "evensplit_01",
#   PatchAffinities %in% c("rep_0.5", "rep_0", "rep_1", "rep_0.25", "rep_0.75")
# ) # IDs 26, 28, 56, 58, 86, 88, 116, 118, 146, 148

# Only 1 type of dynamics.

# Complete Extirpation with doubled events, recognising the affinity decreases
# likelihood of success.

# parameterChoices <- matrix(byrow = TRUE, c(
#   # For each Pool and Dynamics, 2 Histories.
#   # 1 Patch Systems first.
#   # poolpatch, Seed, dynamics, Seed, events, Seed, dispersal, distance
#            26,   11,        1,   11,      2,   14,        NA,        2,
#            26,   11,        1,   11,      2,   15,        NA,        2,
#            26,   12,        1,   12,      2,   16,        NA,        2,
#            26,   12,        1,   12,      2,   17,        NA,        2,
#            26,   11,        1,   11,      2,   14,        NA,        3,
#            26,   11,        1,   11,      2,   15,        NA,        3,
#            26,   12,        1,   12,      2,   16,        NA,        3,
#            26,   12,        1,   12,      2,   17,        NA,        3,
#   # ETC.
# ))
# But we now need to automate this.

# Parameters: #################################################################
parameterChoices <- dplyr::bind_rows(
  # 1 Patch Systems first.
  expand.grid(pp = c(26, 56, 86, 116, 146), dyn = 1,
              events = 2, dispersal = NA, distance = c(2, 3)),
  # 2 Patch Systems with varying dispersal.
  expand.grid(pp = c(28, 58, 88, 118, 148), dyn = 1,
              events = 2, dispersal = c(NA, 15, 10), distance = c(2, 3))
)
# Each sim type should have 2 pools with associated matrices, and each pool
# should in turn have 2 histories assigned to it to identify some stochasticity.
# So we'll duplicate the existing parameter choices and then assign seeds.
parameterChoices <- dplyr::bind_rows(
  parameterChoices, parameterChoices
) %>% dplyr::group_by(
  pp
) %>% dplyr::mutate(
  # ppSeed = c(11, 12, 11, 12, 13, 14, 13, 14, 13, 14,...)
  ppSeed = dplyr::cur_group_id()*2-1 + 9
) %>% dplyr::ungroup(
) %>% dplyr::mutate(
  ppSeed = ppSeed + duplicated(.),
  dynSeed = ppSeed
)

parameterChoices <- dplyr::bind_rows(
  parameterChoices, parameterChoices
) %>% dplyr::group_by(
  pp, ppSeed
) %>% dplyr::mutate(
  # Repeats across distances and dispersals
  eventsSeed = dplyr::cur_group_id()*2-1 + 12
) %>% dplyr::ungroup(
) %>% dplyr::mutate(
  eventsSeed = eventsSeed + duplicated(.)
) %>% dplyr::select(
  dplyr::starts_with("pp"),
  dplyr::starts_with("dyn"),
  dplyr::starts_with("events"),
  dplyr::everything()
) %>% dplyr::arrange(
  dplyr::across(dplyr::everything())
)

# Run across each row of parameterChoices: ####################################
clust <- parallel::makeCluster(4)
doParallel::registerDoParallel(clust)
toExport <- unlist(lapply(
  1:5, # Non-seed columns of parameterChoices
  function(id, pcs, dicts) {
    indices <- unique(unlist(pcs[, id]))
    dictChoices <- dicts[[id]][indices, ]
    stringCols <- unlist(lapply(dictChoices, function(x) all(is.character(x))))
    dictChoices[, stringCols] %>% unlist()
  },
  pcs = parameterChoices %>% dplyr::select(-dplyr::ends_with("Seed")),
  dicts = list(poolpatchDictionaryOrigin, dynamicsDictionaryOrigin,
               eventsDictionaryOrigin, dispersalDictionaryOrigin,
               distanceDictionaryOrigin)
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
    returnResults = FALSE,
    saveResults = TRUE,
    skipIfSaveExists = TRUE
  )
}

parallel::stopCluster(clust)
