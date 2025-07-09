load("diversitiesFlattened9a9_subset4Richness.RData")

# Problems with X11
options(bitmapType = "cairo")

# Grey interval that we compute over, usually after intervention (~50%)
# If second number is less than 1, we lose persistent species.
end <- c(0.602, 0.9045) # Aiming for 20000 - 30000. These go ~0.0003% away.

# Libraries: ##################################################################
library(ggplot2)
library(ggpubr)
library(tidytable)
library(tidygraph)
library(ggraph)

source("TimeSpaceAndTimeSeries-9-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')
source("CalculateTrophicStructure.R") # Calculator creator.
source("toCheddar.R") # Updated function.

# Resources: ##################################################################
interventionMatrix <- matrix(
  c("(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
    "(0.25)->(0)", "(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
    "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
    "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)", "(0.75)->(1)",
    "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)"),
  byrow = TRUE, nrow = 5)

# Functions: ##################################################################
### Manipulation: #############################################################
changeAffinityLevels <- function(df) {
  df %>% tidytable::mutate(
    SpeciesAffinity = tidytable::case_when(
      SpeciesAffinity == "rep_0" ~ "100% 0",
      SpeciesAffinity == "evensplit_01" ~ "50% 0, 50% 1",
      SpeciesAffinity == "runif" ~ "Uniform(0, 1)",
      TRUE ~ SpeciesAffinity
    ),
    SpeciesAffinity = factor(SpeciesAffinity, levels = c(
      "100% 0", "50% 0, 50% 1", "Uniform(0, 1)"
    ), ordered = TRUE)
  )
}

changeInterventionLevels <- function(df) {
  df %>% tidytable::mutate(
    Intervention = factor(
      Intervention,
      levels = t(interventionMatrix)[1:prod(dim(interventionMatrix))], 
      ordered = TRUE
    ), 
    InterventionInitial = tidytable::case_when(
      Intervention %in% interventionMatrix[1, ] ~ diag(interventionMatrix)[1],
      Intervention %in% interventionMatrix[2, ] ~ diag(interventionMatrix)[2],
      Intervention %in% interventionMatrix[3, ] ~ diag(interventionMatrix)[3],
      Intervention %in% interventionMatrix[4, ] ~ diag(interventionMatrix)[4],
      Intervention %in% interventionMatrix[5, ] ~ diag(interventionMatrix)[5],
      TRUE ~ NA_character_
    ),
    InterventionInitial = factor(
      InterventionInitial,
      levels = c(
        diag(interventionMatrix)
      ), ordered = TRUE
    ),
    InterventionFinal = tidytable::case_when(
      Intervention %in% interventionMatrix[, 1] ~ diag(interventionMatrix)[1],
      Intervention %in% interventionMatrix[, 2] ~ diag(interventionMatrix)[2],
      Intervention %in% interventionMatrix[, 3] ~ diag(interventionMatrix)[3],
      Intervention %in% interventionMatrix[, 4] ~ diag(interventionMatrix)[4],
      Intervention %in% interventionMatrix[, 5] ~ diag(interventionMatrix)[5],
      TRUE ~ NA_character_
    ),
    InterventionFinal = factor(
      InterventionFinal,
      levels = c(
        diag(interventionMatrix)
      ), ordered = TRUE
    )
  )
}

# colors: #####################################################################
#                                (0),   (0.5),   (1)
colorPalette <- c(#              Cyan, Magenta, Yellow, Black
  "(0)" =        DescTools::CmykToRgb(1,   0,   0,   0.25),
  "(0)->(0.5)" = DescTools::CmykToRgb(1,   0.5, 0,   0.25),
  "(0)->(1)" =   DescTools::CmykToRgb(1,   0,   0.75, 0.25),
  "(0.5)" =      DescTools::CmykToRgb(0,   1,   0,   0.25),
  "(0.5)->(0)" = DescTools::CmykToRgb(0.5, 1,   0,   0.25),
  "(0.5)->(1)" = DescTools::CmykToRgb(0,   1,   0.75, 0.25),
  "(1)" =        DescTools::CmykToRgb(0,   0,   1,   0.25),
  "(1)->(0)" =   DescTools::CmykToRgb(0.5, 0,   1,   0.25),
  "(1)->(0.5)" = DescTools::CmykToRgb(0,   0.5, 1,   0.25)
)

linetypePalette <- c(
  "100% 0" = "solid",
  "50% 0, 50% 1" = "longdash",
  "Uniform(0, 1)" = "dotdash"
)

# renames: ####################################################################
# For presentation (i.e., "Arrival" is a working term, but not 100% accurate.)
externalNames <- c(
  "Arrival"         = "Colonisation",
  "Failed Arrival"  = "Failure",
  "Present"         = "Present",
  "Dispersal"       = "Adjacent",
  "Extinct"         = "Neutral Ext.",
  "Dynamic Loss"    = "Dynamic Ext.",
  "EndOfSimulation" = "Persistent",
  "NA"              = "NA"
)

# START HERE: ##################################################################
### Strings: ###################################################################
# Enhance readability, from 9g TablePlots
diversitiesInterventionStrings <- ColExt %>% tidytable::select(
  Affinity, PoolPatch, InterventionPatchType
) %>% tidytable::distinct(
) %>% tidytable::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    Affinity, PoolPatch, InterventionPatchType
  ))
)

### End times: #################################################################
# Work out the end times so we can truncate the simulations
# so that we are making sure our comparisons are equivalent.
endTimes <- ColExt %>% tidytable::rename(
  DispersalParam = Dispersal
) %>% tidytable::filter(
  EventType == "EndOfSimulation"
) %>% tidytable::select(
  Times, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, Affinity, AffinitySeed
) %>% tidytable::distinct(
) %>% tidytable::group_by(
  # One of these had an early stop. We "fix" it by going to its descendants.
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, Affinity, AffinitySeed
) %>% tidytable::summarise(
  Times = max(Times),
  .groups = "drop"
) %>% tidytable::mutate( # In the plots:
  Start = end[1] * Times, # Neglect anything with an out time before this.
  Stop = end[2] * Times # Neglect anything with an in time after this.
)

# Persistences: ################################################################
Pers <- ColExt %>% tidytable::rename(
  DispersalParam = Dispersal
) %>% tidytable::filter(
  EventType != "Present", Success # False Arrivals might mess this up.
) %>% tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, Affinity,
  AffinityBins, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed,
  Events, EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance
) %>% tidytable::mutate(
  InNumber = ifelse(EventType == "Arrival" | EventType == "Dispersal", 1, 0),
  InNumber = cumsum(InNumber)
) %>% tidytable::ungroup(
) %>% tidytable::pivot_wider(
  values_from = "Times",
  names_from = EventType,
  id_cols = c(
    # Species/Environment/Simulation Identifiers
    Species, Environment, SpeciesType, Size, ReproductionRate, Speed, Affinity,
    AffinityBins, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed,
    Events, EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
    NicheDistance, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
    InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
    InterventionNicheDistance,
    # And the true ID variable!
    InNumber
  ),
  values_fill = NA
) %>% tidytable::mutate(
  In = ifelse(is.na(Dispersal), Arrival, Dispersal),
  Out = ifelse(is.na(Extinct),
               ifelse(is.na(`Dynamic Loss`),
                      EndOfSimulation,
                      `Dynamic Loss`),
               Extinct),
  InType = externalNames[ifelse(is.na(Dispersal), "Arrival", "Dispersal")],
  OutType = externalNames[ifelse(is.na(Extinct),
                                 ifelse(is.na(`Dynamic Loss`),
                                        "EndOfSimulation",
                                        "Dynamic Loss"),
                                 "Extinct")],
  Persistence = Out - In,
  # Enhance Readability:
  SpeciesAffinity =
    affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
) %>% changeAffinityLevels(
) %>% tidytable::select(
  -Dispersal, -Arrival, -Extinct, -`Dynamic Loss`, -EndOfSimulation
) %>% tidytable::left_join(
  diversitiesInterventionStrings,
  by = c("Affinity", "PoolPatch", "InterventionPatchType"),
  multiple = "all"
) %>% changeInterventionLevels(
) %>% tidytable::left_join(
  endTimes
)

# Diversities: ################################################################
diversitiesRichness <- diversitiesRichness %>% changeAffinityLevels(
) %>% changeInterventionLevels(
)

### Example Networks: #########################################################
#### Load:
targetSeed <- 343
targetDir <- dir(pattern = paste0(
  "TSTS_Simulations_142486-4929_", targetSeed, "-", targetSeed, "_2025-01-21"
),
full.names = T)
stopifnot(length(targetDir) == 1)

targetFiles <- dir(targetDir, pattern = "(Sim|Int)",
                   full.names = T)
# Restrict to source simulations with affinities 0, 0.5, or 1
targetFiles <- grep(x = targetFiles,
                    pattern = "142486-4929-28-1-NA-5-(1|6|7|15|20|21|29|34|35)_",
                    value = TRUE)
targetFilesSim <- grep(x = targetFiles,
                       pattern = "_Simulation_",
                       fixed = TRUE,
                       value = TRUE)
targetFilesInt <- grep(x = targetFiles,
                       pattern = "Int",
                       fixed = TRUE,
                       value = TRUE)
# Restrict interventions to the same
targetFilesInt <- grep(x = targetFilesInt,
                       pattern = "_11(1|3|5)-",
                       fixed = FALSE,
                       value = TRUE)

targetFiles <- c(targetFilesSim, targetFilesInt)
targetDivs <- gsub(x = targetFiles,
                   pattern = "_(Simulation|Intervention)_",
                   replacement = "_Diversity_")
targetPool <- dir(targetDir, pattern = "Pool",
                  full.names = T)

targetEnvs <- lapply(targetFiles, new.env)
targetEnvs <- lapply(seq_along(targetEnvs), function(i, e, s, d) {
  load(d[[i]], envir = e[[i]])
  load(s[[i]], envir = e[[i]])
  e[[i]]$Diversity <- flattenDiversity(e[[i]]$Diversity) %>% dplyr::left_join(
    diversitiesInterventionStrings,
    by = c("Affinity", "PoolPatch", "InterventionPatchType"),
    multiple = "all"
  ) %>% dplyr::mutate(
    PoolPatchSeed = targetSeed,
    SpeciesAffinity =
      affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
  ) %>% changeAffinityLevels()
  e[[i]]
},
e = targetEnvs, d = targetDivs, s = targetFiles)
targetEnvsPool <- new.env()
load(targetPool, envir = targetEnvsPool)

#### Handle
targetTimes <- 30000

targetEnvs <- lapply(targetEnvs, function(env) {
  intervention <- #T/F
    !("EffectiveReproductionRate" %in% names(env$result$Ellipsis$Affinity))
  
  env$calculator <- with(
    env$result,
    CalculateTrophicStructure(
      Pool = targetEnvsPool$Pool,
      ReproductionRate =
        if (!intervention) {
          Ellipsis$Affinity$EffectiveReproductionRate
        } else {
          Ellipsis$Affinity$EffectiveReproductionRateIntervention
        },
      NumEnvironments = NumEnvironments,
      InteractionMatrices = targetEnvsPool$InteractionMatrices,
      EliminationThreshold = Parameters$EliminationThreshold,
      LinkThreshold = 0.01
    )
  )
  
  env$trophics <- with(
    env$result,
    lapply(targetTimes, function(t) {
      if (intervention)
        if(t < Ellipsis$Affinity$TimeIntervention / ReactionTime)
          return(NULL)
      timerow <- which.max(Abundance[, 1] / ReactionTime > t) # first match
      retval <- env$calculator(Abundance[timerow, 1], Abundance[timerow, -1])
      retval$Time <- Abundance[timerow, 1] / ReactionTime
      return(retval)
    })
  )
  
  return(env)
})

#### Graph
targetEnvs <- lapply(targetEnvs, function(env) {
  env$graphs <- lapply(env$trophics, function(TandEV) {
    time <- TandEV$Time
    EVs <- TandEV$EdgeVertexLists[[1]] # Probably accidentally wrapped 1 extra.
    list(graphs = lapply(EVs, function(patch) { # Keep one layer for the multi-patch case.
      tidygraph::tbl_graph(
        nodes = patch$Vertices, edges = patch$Edges
      )
    }), Time = time)
  })
  return(env)
})

targetEnvs <- lapply(targetEnvs, function(env) {
  env$layout <- ggraph::create_layout(
    tidygraph::to_undirected(
      env$graphs[[length(env$graphs)]]$graphs[[1]]
    ) %>% tidygraph::convert(tidygraph::to_simple),
    "backbone"
  )
  return(env)
})

targetEnvs <- lapply(targetEnvs, function(env) {
  l <- env$layout
  l_indices <- as.numeric(gsub("s", "", l$node))
  affs <- env$result$Ellipsis$Affinity$SpeciesAffinities
  if (length(unique(affs)) < 4) {
    # l$x <- affs[l_indices] + l$x/length(unique(affs)) # retain some structure
    affs <- factor(affs, ordered = TRUE, levels = sort(unique(affs)))
    shift <- seq_along(levels(affs)) - 1 # input aff returns number of set
    l$x <- l$x - min(l$x) # shift so left side starts at 0.
    l$x <- l$x / max(l$x) # scale so that it is over a unit interval.
    l$x <- l$x + shift[affs[l_indices]] # add unit interval for each aff.
    l$x <- l$x / max(l$x) # scale one more time so over unit interval again.
  } else {
    l$x <- affs[l_indices] # should be enough variation to enable visualisation
  }
  l$y <- log10(targetEnvsPool$Pool$Size[l_indices])
  env$layout <- l
  return(env)
})

#### Plot
targetEnvs <- lapply(targetEnvs, function(env) {
  timelist <- env$graphs
  env$singletonGraphs <- lapply(timelist, function(patchlist) {
    lapply(patchlist$graphs, plotGraph, mainLayout = env$layout)
  })
  return(env)
})

targetEnvsIndex <- do.call(
  rbind,
  lapply(
    targetEnvs,
    function(env)
      env$Diversity[
        1, c("PoolPatchSeed", "SpeciesAffinity",
             "NicheDistance", "Intervention")]
  )
)
targetEnvsIndex <- cbind(ID = 1:nrow(targetEnvsIndex), targetEnvsIndex)

# Main Plots: #################################################################