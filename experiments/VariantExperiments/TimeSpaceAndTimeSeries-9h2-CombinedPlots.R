load("diversitiesFlattened9a9_subset4Richness.RData")
load("ColExt9a9_flat.RData")

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
library(ggforce)

source("TimeSpaceAndTimeSeries-9-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')
source("flattenDiversity.R") # Flatten the data for the networks.
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

### colors: ###################################################################
# Algorithmic: first index = 100%, second index = 50%
#              (0) => 1,0,0, (0.5) => 0,1,0, (1) => 0,0,1
#              (0.25) => 0.5,0.5,0, (0.75) => 0,0.5,0.5
colorPaletteAlg <- function(intervention) {
  intervention <- as.numeric(strsplit(
    gsub(intervention, pattern = "[(|)]", replacement = ""), 
    split = "->")[[1]])
  x <- intervention[1]
  y <- if(length(intervention) == 2) {
    intervention[2] 
  } else {
    intervention[1]
  }
  DescTools::CmykToRgb(
    min(max(0, (0.5-x)/0.5) + 0.5*max(0, (0.5-y)/0.5), 1), 
    min(max(0, (0.5 - abs(x - 0.5))/0.5) 
        + 0.5*max(0, (0.5 - abs(y - 0.5))/0.5), 1), 
    min(max(0, (x-0.5)/0.5)+ 0.5*max(0, (y-0.5)/0.5), 1),
    0.25
  )
}

colorPalette <- sapply(interventionMatrix, colorPaletteAlg)

linetypePalette <- c(
  "100% 0" = "solid",
  "50% 0, 50% 1" = "longdash",
  "Uniform(0, 1)" = "dotdash"
)

### renames: ##################################################################
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

# Functions: ##################################################################
### Manipulation: #############################################################
changeAffinityLevels <- function(df) {
  df |> tidytable::mutate(
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
  df |> tidytable::mutate(
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

unifyAffinityBins <- function(., n = 5) {
  tidytable::separate_wider_delim(
    .,
    col = "AffinityBins", names = c("Left", "Right"), delim = ","
  ) |> tidytable::mutate(
    Left = 
      round(as.numeric(gsub(pattern = "^[(]", replacement = "", x = Left))*n)/n,
    Right = 
      round(as.numeric(gsub(pattern = "\\]$", replacement = "", x = Right))*n)/n,
    AffinityBins = ifelse(
      is.na(Right), as.character(Left),
      paste0("(", Left, ", ", Right, "]")
    )
  )
}

### Plotting: #################################################################
plotMeanAndInner <- function(
  data, CIs = c(0.5, 0.95),
  facets = as.formula(Intervention ~ SpeciesAffinity)
) {
  # Correct for a problem with how to handle NAs by converting to strings
  data$Subset <- ifelse(is.na(data$Subset), "NA", data$Subset)
  
  # Create base with particular attention to grouping structure.
  baseplot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Time, y = Value,
      group = interaction(
        Subset,
        Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
      ),
      color = Intervention, fill = Intervention, linetype = SpeciesAffinity
    )
  )
  
  # Plot each CI overlaid. Overlaying => the innermost have the darkest alpha.
  for (CI in CIs) {
    baseplot <- baseplot + ggplot2::geom_ribbon(
      data = data |> tidytable::mutate(
        Time = round(Time, digits = -2)
      ) |> tidytable::group_by(
        Time, Subset,
        Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
      ) |> tidytable::summarise(
        top = quantile(Value, probs = CI+(1-CI)/2, na.rm = TRUE),
        bot = quantile(Value, probs = (1-CI)-(1-CI)/2, na.rm = TRUE)
      ), mapping = ggplot2::aes(
        x = Time, ymin = bot, ymax = top,
        group = interaction(
          Subset,
          Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
        ),
        fill = Intervention,
        color = Intervention
      ), inherit.aes = FALSE,
      alpha = 0.25, linewidth = 0.25 #, linetype = "dotted"
    )
  }
  
  # Add an average line and handle the meta-details.
  baseplot <- baseplot + ggplot2::geom_line(
    data = data |> tidytable::mutate(
      Time = round(Time, digits = -2)
    ) |> tidytable::group_by(
      Time, Subset,
      Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
    ) |> tidytable::summarise(
      Value = mean(Value, na.rm = TRUE)
    )
  ) + ggplot2::facet_grid(
    facets
  ) + ggplot2::scale_color_manual(
    values = colorPalette, aesthetics = c("color", "fill"),
    name = "Local Land-use"
  ) + ggplot2::scale_linetype_manual(
    name = "Species Preferences",
    values = linetypePalette
  ) + ggplot2::theme_minimal(
  ) + ggplot2::guides(
    color = if (length(CIs)>0) {"none"} else {"legend"},
    fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
  )
  
  return(baseplot)
}

# Note: unlike the other plots, this one is called inside the environments.
#       The other plots are called inside the createOverviewFigure function.
plotGraph <- function(graph, mainLayout, legends = FALSE) {
  obj <- ggraph::ggraph(
    graph = graph %N>% mutate(
      nodesize = (log10(N)+5)/10+1
    ),
    layout = graph %N>% data.frame(
    ) |> select(node) |> left_join(
      mainLayout |> select(x, y, node)
    )
  ) + ggraph::geom_edge_diagonal(
    mapping = aes(
      color = Type,
      linetype = Type,
      alpha = log10(effectNormalised),
      end_cap = circle(node2.nodesize+2, 'pt')
    ),
    arrow = arrow(length = unit(2, 'mm')), 
    show.legend = legends
  ) + ggraph::geom_node_point(
    mapping = aes(
      color = Type,
      size = nodesize
    ), 
    show.legend = legends
    # ) + ggplot2::geom_hline(
    #   yintercept = -1, linetype = "dashed", color = "black"
  ) + theme_minimal(
  ) + ylab(
    "Log10(Size)"
  ) + xlab(
    "Land-use Preference"
  ) + scale_color_manual(
    values = c("limegreen", "goldenrod2")
  ) + ggraph::scale_edge_color_manual(
    values = rev(c("limegreen", "goldenrod2"))
  ) + scale_size(
    range = c(0.5, 4)
    # limits = c(10^-5, 10^5)#, trans = "log10"
  ) + coord_cartesian(
    x = c(0, 1), y = c(-2, 0.25),
    clip = "off"
  )
  
  # if (!legends) {
  #   obj <- obj + ggplot2::theme(
  #     legend.position = "none"
  #   )
  # }
  
  return(obj)
}

# START HERE: ##################################################################
### Strings: ###################################################################
# Enhance readability, from 9g TablePlots
diversitiesInterventionStrings <- ColExt |> tidytable::select(
  Affinity, PoolPatch, InterventionPatchType
) |> tidytable::distinct(
) |> tidytable::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    Affinity, PoolPatch, InterventionPatchType
  ))
)

### End times: #################################################################
# Work out the end times so we can truncate the simulations
# so that we are making sure our comparisons are equivalent.
endTimes <- ColExt |> tidytable::rename(
  DispersalParam = Dispersal
) |> tidytable::filter(
  EventType == "EndOfSimulation"
) |> tidytable::select(
  Times, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, Affinity, AffinitySeed
) |> tidytable::distinct(
) |> tidytable::group_by(
  # One of these had an early stop. We "fix" it by going to its descendants.
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, Affinity, AffinitySeed
) |> tidytable::summarise(
  Times = max(Times),
  .groups = "drop"
) |> tidytable::mutate( # In the plots:
  Start = end[1] * Times, # Neglect anything with an out time before this.
  Stop = end[2] * Times # Neglect anything with an in time after this.
)

# Persistences: ################################################################
Pers <- ColExt |> tidytable::rename(
  DispersalParam = Dispersal
) |> tidytable::filter(
  EventType != "Present", Success # False Arrivals might mess this up.
) |> tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, Affinity,
  AffinityBins, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed,
  Events, EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance
) |> tidytable::mutate(
  InNumber = ifelse(EventType == "Arrival" | EventType == "Dispersal", 1, 0),
  InNumber = cumsum(InNumber)
) |> tidytable::ungroup(
) |> tidytable::pivot_wider(
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
) |> tidytable::mutate(
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
) |> changeAffinityLevels(
) |> tidytable::select(
  -Dispersal, -Arrival, -Extinct, -`Dynamic Loss`, -EndOfSimulation
) |> tidytable::left_join(
  diversitiesInterventionStrings,
  by = c("Affinity", "PoolPatch", "InterventionPatchType"),
  multiple = "all"
) |> changeInterventionLevels(
) |> tidytable::left_join(
  endTimes
) |> unifyAffinityBins()

# Diversities: ################################################################
diversitiesRichness <- diversitiesRichness |> changeAffinityLevels(
) |> changeInterventionLevels(
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
  e[[i]]$Diversity <- flattenDiversity(e[[i]]$Diversity) |> dplyr::left_join(
    diversitiesInterventionStrings,
    by = c("Affinity", "PoolPatch", "InterventionPatchType"),
    multiple = "all"
  ) |> dplyr::mutate(
    PoolPatchSeed = targetSeed,
    SpeciesAffinity =
      affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
  ) |> changeAffinityLevels()
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
    ) |> tidygraph::convert(tidygraph::to_simple),
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
### Plot 2:####################################################################
# a=>b&c
newplot2_dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot2_indices <- targetEnvsIndex |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)

newplot2_dataC <- Pers |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
) |> tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) |> tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) |> tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, AffinityBins,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, DispersalParam, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) |> tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
)

##### a: ######################################################################
newplot2_a <- plotMeanAndInner(
  newplot2_dataA |> tidytable::filter(
    Intervention %in% c("(0)", "(0.5)", "(1)")
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = newplot2_dataA |> tidytable::filter(
    PoolPatchSeed == targetSeed,
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    abs(Time - targetTimes) == min(abs(Time - targetTimes))
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 3),
  fill = ggplot2::guide_legend(ncol = 3)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), ylim = c(0, 42), expand = FALSE
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      newplot2_indices$ID[1]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 7, ymax = 17
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      newplot2_indices$ID[2]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 18, ymax = 28
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      newplot2_indices$ID[3]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 29, ymax = 39
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(newplot2_dataA$Start),
  xmax = max(newplot2_dataA$Stop),
  ymin = 0, ymax = max(newplot2_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.3, 0.09),
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
newplot2_b <- ggplot2::ggplot(
  newplot2_dataA |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = 0, ymax = max(newplot2_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9)
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::geom_jitter(
  alpha = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Local Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Local Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::annotate(
  "text", x = c(1.5, 4.5), y = 5, label = c("High\nMatch", "Low\nMatch")
)

##### c: ######################################################################
newplot2_c <- ggplot2::ggplot(
  newplot2_dataC,
  ggplot2::aes(
    y = Persistence,
    x = Intervention,
    color = Intervention,
    group = interaction(Intervention, SpeciesType),
    fill = SpeciesType
  )
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = -Inf, ymax = Inf,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9), show.legend = FALSE, 
  scale = "count", draw_quantiles = 0.5
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("limegreen", "goldenrod2")
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  tag = "c)", x = "Habitat"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
) + ggplot2::facet_grid(
  SpeciesType ~ .
) + ggplot2::scale_x_discrete(
  breaks = c("(0)", "(0.5)", "(1)")
) + ggplot2::labs(
  x = "Local Land-use"
)

newplot2 <- ggpubr::ggarrange(
  plotlist = list(
    newplot2_a,
    newplot2_b,
    newplot2_c
  ), nrow = 1, widths = c(0.5, 0.27, 0.23)
)

ggplot2::ggsave(plot = newplot2, filename = "Figure2_Prototype6.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

### Plot 3:####################################################################
newplot3_dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesAffinity != "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot3_dataB <- Pers |> tidytable::filter(
  SpeciesAffinity != "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
) |> tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) |> tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) |> tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, AffinityBins,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, DispersalParam, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) |> tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
)

##### a: ######################################################################
newplot3_a <- ggplot2::ggplot(
  newplot3_dataA |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) |> tidytable::summarise(
    Value = mean(Value)
  ),
  ggplot2::aes(
    x = Intervention,
    y = Value,
    color = Intervention
  )
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = 0, ymax = max(newplot2_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9), scale = "count"
) + ggplot2::geom_boxplot(
  notch = TRUE, outlier.size = 0,
  position = ggplot2::position_dodge(0.9),
  width = 0.28
) + ggplot2::geom_jitter(
  alpha = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Local Land-use"
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Local Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesAffinity ~ .
)

##### b: ######################################################################
newplot3_b <- ggplot2::ggplot(
  newplot3_dataB,
  ggplot2::aes(
    y = Persistence,
    x = Intervention,
    color = Intervention,
    group = interaction(Intervention, SpeciesType, AffinityBins),
    fill = SpeciesType
  )
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 6,
  ymin = -Inf, ymax = Inf,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(0.9), show.legend = FALSE, 
  scale = "count", draw_quantiles = 0.5
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("limegreen", "goldenrod2")
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  tag = "b)", x = "Habitat"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::facet_grid(
  SpeciesType ~ SpeciesAffinity
) + ggplot2::labs(
  x = "Local Land-use\n(Binned and Ordered by Preference)"
)

newplot3 <- ggpubr::ggarrange(
  plotlist = list(
    newplot3_a,
    newplot3_b
  ), nrow = 2, widths = c(0.5, 0.5)
)

ggplot2::ggsave(plot = newplot3, filename = "Figure3_Prototype4.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

### Plot 4:####################################################################
# Need to contrast with 2a (Richness). Long and short time scales.
# So starting from the top again, we want to construct the image to lead 
# naturally to the comparison of the two different time spans for the same
# statistic, which probably looks more like differences of various sorts.
# Leaning on some of what we had originally set off to do, we might be able
# to characterise it as 
#     time before-after short term, 
#     counterfactual before-after short term,
#     time before-after long term (slope 0),
#     counterfactual before-after long term.
# Regardless, we probably need to re-run things so we have consistent 
# comparisons that we are making, specifically for the time before-after short.
# So:
#  a => b + c, with d "contained within" a
#  d => e + f
# long time a: regular time; b: temporal comparison; c: counterfactual compare
# short time d: rescaled time; e: temporal compare; f: counterfactual compare
#
# Then it might be a good idea to summarise the b vs c and e vs f in the text
# via the differences.

newplot4_dataA <- diversitiesRichness |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot4_dataC <- newplot4_dataA |> tidytable::group_by(
  PoolPatch:AffinitySeed
) |> tidytable::mutate(
  Time = Time - Time[which.min(Time[!is.na(InterventionPatchSeed)])]
) |> tidytable::...(
  # Time 0 since intervention is not "observed" in the original runs, but
  # we have a workaround, since it should have been observed identically
  # across the different intervention runs!
)


##### a: ######################################################################
newplot4_a <- plotMeanAndInner(
  newplot4_dataA |> tidytable::filter(
    Intervention %in% c("(0)", "(0)->(0.5)", "(0.5)", "(0)->(1)", "(1)")
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 3),
  fill = ggplot2::guide_legend(ncol = 3)
) + ggplot2::coord_cartesian(
  xlim = c(0, 31000), ylim = c(0, 42), expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 16300,
  xmax = 16400,
  ymin = 0, ymax = 31,
  fill = "grey",
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.3, 0.09),
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
  # ) + ggforce::facet_zoom(
  #   xlim = c(16000, 17000),
  #   shrink = FALSE
) + ggplot2::geom_vline(
  xintercept = c(16300, 30000), color = "black", linetype = "dashed"
)


##### b: ######################################################################
# Careful with 16300; obtained by checking the table of times to see what was a
# safe largest value before intervention.
newplot4_b <- rbind(
  # Temporal Substitution
  newplot4_dataA |> tidytable::group_by(
    SpeciesAffinity, Intervention, PoolPatchSeed
  ) |> tidytable::filter(
    Time == max(16300, min(Time)) | Time == 30000
  ) |> tidytable::group_by(
    SpeciesAffinity, PoolPatchSeed, 
    Intervention, InterventionInitial, InterventionFinal
  ) |> tidytable::arrange(
    Time
  ) |> tidytable::summarise(
    Value = Value[2] - Value[1],
    Method = "Temporal"
  ),
  # Counterfactual comparison TODO, ONLY 100% VALID AFTER RE-RUNS WITH FIXED 
  #                                 POOL PREFERENCE ASSIGNMENTS.
  newplot4_dataA |> tidytable::group_by(
    SpeciesAffinity, Intervention, PoolPatchSeed
  ) |> tidytable::filter(
    Time == 30000
  ) |> tidytable::group_by(
    SpeciesAffinity, PoolPatchSeed, 
    # InterventionFinal # What if you had always been in your final state?
    InterventionInitial # What if you had stayed in your initial state?
  ) |> tidytable::mutate(
    Value = Value - Value[InterventionInitial == InterventionFinal],
    Method = "Counterfactual"
  ) |> tidytable::select(
    SpeciesAffinity, PoolPatchSeed,
    Intervention, InterventionInitial, InterventionFinal,
    Value, Method
  )
) |> tidytable::filter(
  InterventionInitial == "(0)"
) |> ggplot2::ggplot(
  ggplot2::aes(y = Value, 
               # x = InterventionInitial, 
               # group = interaction(Method, InterventionInitial),
               # color = InterventionInitial
               x = InterventionFinal,
               group = interaction(Method, InterventionFinal),
               color = InterventionFinal
  )
) + ggplot2::geom_violin(
  scale = "count", show.legend = FALSE
  # ) + ggplot2::geom_boxplot(
) + ggplot2::facet_grid(
  InterventionInitial ~ .,#Method,
  scales = "free_x"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  tag = "b)", 
  x = "Final Local Land-Use",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank(),
  panel.border = ggplot2::element_rect(
    linetype = "dashed", color = "black", fill = NA
  )
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_y_continuous(
  sec.axis = ggplot2::sec_axis(
    trans = function(x) x, name = "Initial Local Land-Use", 
    breaks = NULL, labels = NA
  )
)

##### c: ######################################################################

newplot4_c <- ggplot2::ggplot(
  # Mutate the data to start at the time of intervention, which is unique
  # for each PoolPatchSeed (since it is defined from the Events which are
  # paired precisely with each PoolPatchSeed).
  # The trick is that the intervention time isn't in the non-intervention case.
  newplot4_dataC |> tidytable::filter(
    Intervention %in% c(
      "(0)", 
      "(0)->(0.5)", 
      "(0)->(1)"
    )
  ) |> tidytable::group_by(
    PoolPatch:AffinitySeed
  ) |> tidytable::mutate(
    Time = Time - Time[which.min(Time[!is.na(InterventionPatchSeed)])]
  ), 
  ggplot2::aes(x = Time, y = Value, color = Intervention, 
               group = interaction(PoolPatchSeed, Intervention))
) + ggplot2::geom_line(
  show.legend = FALSE
) + ggplot2::labs(
  y = "Richness",
  x = "Time since Intervention"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 3),
  fill = ggplot2::guide_legend(ncol = 3)
) + ggplot2::coord_cartesian(
  xlim = c(-3, 23), ylim = c(0, 31), expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 0,
  xmax = 20,
  ymin = 0, ymax = 31,
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "c)"
) + ggplot2::theme(
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::geom_vline(
  xintercept = c(0, 20),
  linetype = "dotted"
)

##### d: ######################################################################

newplot4_d <- rbind(
  # Temporal Substitution
  newplot4_dataC |> tidytable::group_by(
    SpeciesAffinity, Intervention, PoolPatchSeed
  ) |> tidytable::filter(
    Time >= 0, Time < 20 # Prevent 3 from the same run.
  ) |> tidytable::group_by(
    SpeciesAffinity, PoolPatchSeed,
    Intervention, InterventionInitial, InterventionFinal
  ) |> tidytable::arrange(
    Time
  ) |> tidytable::summarise(
    Value = Value - Value[1],
    Method = "Temporal"
  ),
  # Counterfactual comparison TODO, ONLY 100% VALID AFTER RE-RUNS WITH FIXED
  #                                 POOL PREFERENCE ASSIGNMENTS.
  newplot4_dataC |> tidytable::group_by(
    SpeciesAffinity, Intervention, PoolPatchSeed
  ) |> tidytable::filter(
    Time >= 0, Time < 20 # Prevent 3 from the same run.
  ) |> tidytable::group_by(
    SpeciesAffinity, PoolPatchSeed, # Time, if times were matched up perfectly.
    # InterventionFinal # What if you had always been in your final state?
    InterventionInitial # What if you had stayed in your initial state?
  ) |> tidytable::mutate(
    Value = Value - Value[InterventionInitial == InterventionFinal],
    Method = "Counterfactual"
  ) |> tidytable::select(
    SpeciesAffinity, PoolPatchSeed,
    Intervention, InterventionInitial, InterventionFinal,
    Value, Method
  )
) |> tidytable::filter(
  InterventionInitial == "(0)"
) |> ggplot2::ggplot(
  ggplot2::aes(y = Value,
               # x = InterventionInitial,
               # group = interaction(Method, InterventionInitial),
               # color = InterventionInitial
               x = InterventionFinal,
               group = interaction(Method, InterventionFinal),
               color = InterventionFinal
  )
) + ggplot2::geom_violin(
  scale = "count", show.legend = FALSE
  # ) + ggplot2::geom_boxplot(
) + ggplot2::facet_grid(
  InterventionInitial ~ .,#Method,
  scales = "free_x"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  tag = "b)",
  x = "Final Local Land-Use",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank(),
  panel.border = ggplot2::element_rect(
    linetype = "dashed", color = "black", fill = NA
  )
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_y_continuous(
  sec.axis = ggplot2::sec_axis(
    trans = function(x) x, name = "Initial Local Land-Use",
    breaks = NULL, labels = NA
  )
)


### Plot 5: ###################################################################
# Pseudo-relaxation time of the system from the intervention to its new final
# state, characterised as the difference between counterfactual always in final
# state and the intervention to the final state.
# Because of the rescalings, if we index it by the original time of intervention
# and then compare, we should be able to observe roughly the same relaxation
# time (so long as we don't perform the second rescaling to have interventions
# fixed at 0.5?).