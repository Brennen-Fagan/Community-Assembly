# Biomass script for trying to determine sensible values for JP's logistic.

datfolders <- dir(path = "Deprecated", 
                  pattern = "TSTS_Simulations_.+2024-11-30$", 
                  full.names = TRUE) # Regex

# Problems with X11
options(bitmapType = "cairo")

# Libraries: ##################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
source("TimeSpaceAndTimeSeries-9-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')

# Not a finished function!
interventionNamingScheme <- function(aff, ppa, ipt) {
  aDO <- affinityDictionaryOrigin[aff, ]
  ppDO <- poolpatchDictionaryOrigin[ppa, ]
  if (explicit <- grepl(pattern = "rep", aDO$PatchAffinities)) {
    initState <-
      paste0("(", paste( # NOT PRETTY FOR 10, MAY WANT TO JUST REPORT FUNC CALL
        vals <- retrieveFunction(aDO$PatchAffinities)(ppDO$NumberEnvironments),
        collapse = ", "), ")")
  } else {
    initState <-
      paste0(aDO$PatchAffinities, "(", ppDO$NumberEnvironments, ")")
  }
  if(is.na(ipt)) {return(initState)}
  ipDO <- interventionPatchDictionaryOrigin[ipt, ]
  if (ppDO$NumberEnvironments == 1) {
    finState <- paste0(
      "(", retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments), ")"
    )
  } else if (is.na(ipDO$InterventionLocation) ||
             !explicit ||
             !grepl(pattern = "rep", ipDO$PatchAffinities)) {
    finState <- # InterventionPercentage is a bit of a misnomer!
      paste0(ipDO$InterventionPercentage * 100, "%", ipDO$PatchAffinities)
  } else if (ipDO$InterventionLocation == 0) {# Left
    valsnew <- retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments)
    finState <- paste0("(", paste(
      valsnew[1:(ppDO$NumberEnvironments*ipDO$InterventionPercentage)],
      vals[
        (ppDO$NumberEnvironments*ipDO$InterventionPercentage + 1):
          ppDO$NumberEnvironments],
      collapse = ", ", sep = ", "), ")")
  } else if (ipDO$InterventionLocation == 1) {# Right
    valsnew <- retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments)
    finState <- paste0("(", paste(
      vals[1:(ppDO$NumberEnvironments*(1 - ipDO$InterventionPercentage))],
      valsnew[
        (ppDO$NumberEnvironments*(1 - ipDO$InterventionPercentage) + 1):
          ppDO$NumberEnvironments],
      collapse = ", ", sep = ", "), ")")
  }
  return(paste0(initState, "->", finState))
}

biomasses <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Diversity"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    
    div <- get(names)
    biomass <- div$Presence %>% dplyr::group_by(
      Time, Environment
    ) %>% dplyr::summarise(
      biomass = sum(Abundance * Size),
      biomassBasal = sum(Abundance * Size * (Type == "Basal")),
      richness = dplyr::n(),
      richnessBasal = sum(Type == "Basal"),
      .groups = "drop"
    )
    
    return(list(biomass = biomass, Ellipsis = div$Ellipsis, 
                "Dir" = dirname(x), "File" = basename(x)))
  })

biomassesFlattened <- do.call(rbind, lapply(biomasses, function(d) {
  if ("FullID" %in% names(d$Ellipsis)) {
    id <- strsplit(
      strsplit(d$Ellipsis$FullID, "_", fixed = TRUE)[[1]], # Split seeds off.
      "-", fixed = TRUE)
  } else if ("ID" %in% names(d$Ellipsis)) {
    id <- strsplit(
      strsplit(d$Ellipsis$ID, "_", fixed = TRUE)[[1]], # Split seeds off.
      "-", fixed = TRUE)
  } else {
    id <- strsplit(
      strsplit(
        strsplit(d$File, ".", fixed = TRUE)[[1]][1], # Remove .RData.
        "_", fixed = TRUE)[[1]][3:4], # Remove TSTS_Type and split seeds off.
      "-", fixed = TRUE # Separate out the id values.
    )
  }
  
  if (length(id) < 3) {
    # I.e., no intervention.
    id[[3]] <- rep(NA, 4)
    id[[4]] <- rep(NA, 2)
  }
  
  d$biomass %>% dplyr::mutate(
    PoolPatch = id[[1]][1],
    PoolPatchSeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    InitialConditions = id[[1]][4],
    InitialConditionsSeed = id[[2]][4],
    Dispersal = id[[1]][5],
    NicheDistance = id[[1]][6],
    Affinity = id[[1]][7],
    AffinitySeed = id[[2]][5],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
  )
}))

diversitiesInterventionStrings <- biomassesFlattened %>% dplyr::select(
  Affinity, PoolPatch, InterventionPatchType
) %>% dplyr::distinct(
) %>% dplyr::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    Affinity, PoolPatch, InterventionPatchType
  ))
)

biomassesFlattened <- biomassesFlattened %>% dplyr::left_join(
  diversitiesInterventionStrings,
  by = c("Affinity", "PoolPatch", "InterventionPatchType")
)

biomassesFlattened <- biomassesFlattened %>% dplyr::mutate(
  SpeciesAffinity =
    affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
)

ggplot2::ggplot(
  biomassesFlattened, 
  ggplot2::aes(
    x = Time, y = biomass, 
    group = interaction(
      PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, 
      Events, EventsSeed, InitialConditions, InitialConditionsSeed, 
      Dispersal, NicheDistance, Affinity, AffinitySeed, 
      InterventionPatchType, InterventionPatchSeed, InterventionTimeType, 
      InterventionDispersal, InterventionTimeType, InterventionDispersal, 
      InterventionNicheDistance, Intervention, SpeciesAffinity),
    color = NicheDistance
  )
) + ggplot2::geom_line(
) + ggplot2::scale_y_log10(
) + ggplot2::facet_grid(
  SpeciesAffinity ~ Intervention
) + ggplot2::coord_cartesian(
  ylim = c(1e2, 1e4)
)

ggplot2::ggplot(
  biomassesFlattened, 
  ggplot2::aes(
    x = Time, y = biomassBasal, 
    group = interaction(
      PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, 
      Events, EventsSeed, InitialConditions, InitialConditionsSeed, 
      Dispersal, NicheDistance, Affinity, AffinitySeed, 
      InterventionPatchType, InterventionPatchSeed, InterventionTimeType, 
      InterventionDispersal, InterventionTimeType, InterventionDispersal, 
      InterventionNicheDistance, Intervention, SpeciesAffinity),
    color = NicheDistance
  )
) + ggplot2::geom_line(
) + ggplot2::scale_y_log10(
) + ggplot2::facet_grid(
  SpeciesAffinity ~ Intervention
) + ggplot2::coord_cartesian(
  ylim = c(3e1, 3e3)
)

# list(Basal = 1000) seems the most reasonable, since Basal:Consumer ratio
# and biomass don't have an easy relationship. (See 0's in 1 versus in 0.)
