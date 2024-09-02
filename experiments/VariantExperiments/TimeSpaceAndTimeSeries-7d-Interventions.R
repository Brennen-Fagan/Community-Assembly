library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

directory <- '.'
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Interventions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-7a-SimulationFunction.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-7c-InterventionFunction.R"))

dirTag <- "TSTS_Simulations"
dirDate <- "2024-08-29"
baseTag <-  "TSTS_Simulation" # Note the distinction

runSimulations <- FALSE
source(file.path(directory, "TimeSpaceAndTimeSeries-7b-Simulations.R"))
# For Side Effects, specifically parameterChoices.

stopifnot(exists("parameterChoices"))

parameterChoices <- dplyr::bind_rows(
  parameterChoices %>% dplyr::filter(pp %in% c(26, 28)) %>% dplyr::mutate(
    ip = dplyr::case_when(
      pp == 26 ~ 46, # (0.5) -> (0), all
      pp == 28 ~ 13, # (0.5, 0.5) -> (0.5, 0), last half
      TRUE ~ as.numeric(NA) # Oops!
    )
  ),
  parameterChoices %>% dplyr::mutate(
    ip = dplyr::case_when(
      pp == 26 ~ 49, # (0.5) -> (1), all
      pp == 28 ~ 16, # (0.5, 0.5) -> (0.5, 1), last half
      pp %in% c(56, 146) ~ 45, # (0 or 1) -> (0.5), all
      pp %in% c(58, 148) ~ 12, # (0 or 1 x 2) -> (0 or 1, 0.5) last half
      pp == 86 ~ 48, # (0.25) -> (0.75), all
      pp == 88 ~ 15, # (0.25, 0.25) -> (0.25, 0.75), last half
      pp == 116 ~ 47, # (0.75) -> (0.25), all
      pp == 118 ~ 14, # (0.75, 0.75) -> (0.75, 0.25), last half
      TRUE ~ as.numeric(NA) # Oops!
    )
  )
) %>% dplyr::mutate(
  ipSeed = 1, # No random components
  it = 1, # Dead middle
  itSeed = 1, # No random components
  itDisp = "p", # Use previous values.
  itDist = "p" # Use previous values.
) %>% dplyr::arrange(
  dplyr::across(dplyr::everything())
)

# Run across each row of parameterChoices: ####################################
clust <- parallel::makeCluster(4)
doParallel::registerDoParallel(clust)

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
               eventsDictionaryOrigin, dispersalDictionaryOrigin,
               distanceDictionaryOrigin, interventionPatchDictionaryOrigin,
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
    #   poolpatchDictionaryChoice = pc[1],
    #   poolpatchSeedChoice = pc[2],
    #   dynamicsDictionaryChoice = pc[3],
    #   dynamicsSeedChoice = pc[4],
    #   eventsDictionaryChoice = pc[5],
    #   eventsSeedChoice = pc[6],
    #   dispersalDictionaryChoice = pc[7],
    #   distanceDictionaryChoice = pc[8],
    #   Date = fileDates
    # ),
    ID = file.path(
      paste0(
        dirTag, "_", pc[[1]], "-", pc[[3]],
        "_", pc[[2]], "-", pc[[4]], "_", dirDate
      ),
      paste0(
        baseTag,
        "_", pc[[1]], "-", pc[[3]], "-", pc[[5]], "-", pc[[7]], "-", pc[[8]],
        "_", pc[[2]], "-", pc[[4]], "-", pc[[6]], ".RData"
      )
    ),
    interventionPatchDictionaryChoice = pc[[9]],
    interventionPatchSeedChoice = pc[[10]],
    interventionTimeDictionaryChoice = pc[[11]],
    interventionTimeSeedChoice = pc[[12]],
    interventionDispersalDictionaryChoice = pc[[13]],
    interventionDistanceDictionaryChoice = pc[[14]],
    returnResults = FALSE,
    saveResults = TRUE,
    skipIfSaveExists = TRUE
  )
}

parallel::stopCluster(clust)
