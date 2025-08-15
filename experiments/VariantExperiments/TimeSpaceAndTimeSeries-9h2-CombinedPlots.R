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
source("generateNetworks.R")

# Resources: ##################################################################
interventionMatrix <- matrix(c(
  "(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
  "(0.25)->(0)", "(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
  "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
  "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)", "(0.75)->(1)",
  "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)"
  ),
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
    name = "Habitat's Land-use"
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
      #color = node1.Type, # but then exploit+ between consumers is orange.
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
    values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
  ) + ggraph::scale_edge_color_manual(
    values =
      c("Exploit+" = "darkgreen", "Exploit-" = "burlywood4")
    # c("Basal" = "darkgreen", "Consumer" = "burlywood4")

  ) + scale_size(
    range = c(0.5, 4)
    # limits = c(10^-5, 10^5)#, trans = "log10"
  ) + coord_cartesian(
    x = c(0, 1), y = c(-2, 0.5),
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

# Most efficient to call the function this creates only once!
# This way we can perform (closer to) the bare minimum number of loads.
# This function takes in a specification (essentially, the information from
# the diversities objects, but not the measurements in those objects).
# The specification can be of multiple rows.
# newplot2_a_seed <- "358"; newplot2_a_time <- 25000
newplot2_a_seed <- "374"; newplot2_a_time <- 25000
specification <- diversitiesRichness |> tidytable::select(c(
  # Which network:
  "Time", "Environment1",
  # Which File (Base):
  "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
  "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
  "Dispersal", "NicheDistance", "Affinity", "AffinitySeed",
  # Which File (Intervention):
  "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
  "InterventionTimeSeed", "InterventionDispersal", "InterventionNicheDistance",
  # Ease of Use
  "SpeciesAffinity", "Intervention"
)) |> tidytable::filter(
  # 2A
  (
    SpeciesAffinity == "100% 0" &
      NicheDistance == "5" &
      Intervention %in% c("(0)", "(0.5)", "(1)") &
      PoolPatchSeed %in% as.character(newplot2_a_seed) &
      Time == newplot2_a_time
    # ) | (

  )
) |> tidytable::distinct(
)

exampleNetworks <- generateNetworks(specification)

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

newplot2_dataAS <- diversitiesRichness |> tidytable::filter(
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  !is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot2_indices <- exampleNetworks$Index |> tidytable::filter(
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
  rbind(
    newplot2_dataA |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    newplot2_dataA |> tidytable::filter(
      PoolPatchSeed == newplot2_a_seed,
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time))
    ) |> tidytable::mutate(
      Value = -100
    )
  ), CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = newplot2_dataA |> tidytable::filter(
    PoolPatchSeed == newplot2_a_seed,
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time))
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), ylim = c(0, 42), expand = FALSE
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    exampleNetworks$Envs[[newplot2_indices$ID[1]]]$singletonGraphs[[1]] +
      ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      ) + ggplot2::coord_cartesian(
        xlim = c(NA, NA), ylim = c(-2, 0.5)
      ) # Easiest to probably just not worry about comparing between.
  ),
  xmin = 30500, xmax = 40000, ymin = 7, ymax = 17
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    exampleNetworks$Envs[[newplot2_indices$ID[2]]]$singletonGraphs[[1]] +
      ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      ) + ggplot2::coord_cartesian(
        xlim = c(NA, NA), ylim = c(-2, 0.5)
      ) # Easiest to probably just not worry about comparing between.
  ),
  xmin = 30500, xmax = 40000, ymin = 18, ymax = 28
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    exampleNetworks$Envs[[newplot2_indices$ID[3]]]$singletonGraphs[[1]] +
      ggplot2::theme_void(
      ) + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white")
      ) + ggplot2::coord_cartesian(
        xlim = c(NA, NA), ylim = c(-2, 0.5)
      ) # Easiest to probably just not worry about comparing between.
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
) + ggplot2::geom_segment(
  data = data.frame(
    x = newplot2_a_time+250,
    y = c(10, 22, 39),
    xend = 30500,
    yend = c(11, 22, 36),
    Intervention = c("(0)", "(0.5)", "(1)")
  ),
  mapping = ggplot2::aes(
    x = x, y = y, xend = xend, yend = yend, color = Intervention
  ),
  inherit.aes = FALSE,
  arrow = arrow(length = unit(0.03, "npc"))
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.5, 0.09),
  plot.tag.position = c(0.025, 0.95)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
) + ggplot2::annotate(
  "text", x = 36500, y = c(16, 38), size = 3, lineheight = 0.7,
  label = c("Fully\nAdapted", "Poorly\nAdapted")
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
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
) + ggplot2::geom_line(
  data = ~ summarise(group_by(.x, Intervention), Value = mean(Value)),
  color = "black", group = 1
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 0.95)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::annotate(
  "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted")
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
  scale = "count", draw_quantiles = 0.5, linewidth = 1.3
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("darkgreen", "burlywood4")
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  tag = "c)", x = "Habitat"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 0.95)
) + ggplot2::facet_grid(
  factor(SpeciesType, levels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::scale_x_discrete(
  breaks = c("(0)", "(0.5)", "(1)")
) + ggplot2::labs(
  x = "Habitat's Land-use"
)

newplot2 <- ggpubr::ggarrange(
  plotlist = list(
    newplot2_a,
    newplot2_b,
    newplot2_c
  ), nrow = 1, widths = c(0.5, 0.27, 0.23)
)

newplot2 <- newplot2 + ggplot2::annotate(
  "curve",
  x = 0.35, y = 0.97, xend = c(0.57, 0.87), yend = 0.97,
  curvature = -0.075,
  arrow = arrow(length = unit(0.03, "npc"))
)

ggplot2::ggsave(plot = newplot2, filename = "Figure2_Prototype6.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

### 2 Supplement: #############################################################
##### a: ######################################################################
newplot2_as <- plotMeanAndInner(
  rbind(
    newplot2_dataAS |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    newplot2_dataAS |> tidytable::filter(
      PoolPatchSeed == newplot2_a_seed,
      Intervention %in% c("(0.25)", "(0.75)"),
      abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time)),
      !is.na(Subset)
    ) |> tidytable::mutate(
      Value = -100
    )
  ), CIs = 0.75, facets = as.formula(
    factor(Subset, levels = c("Consumer_0", "Basal_0"), ordered = TRUE) ~ .
  )
) + ggplot2::geom_point(
  data = newplot2_dataAS |> tidytable::filter(
    PoolPatchSeed == newplot2_a_seed,
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time))
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 5),
  fill = ggplot2::guide_legend(ncol = 5)
) + ggplot2::coord_cartesian(
  xlim = c(0, 31000), ylim = c(0, 42), expand = FALSE
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
  legend.position = c(0.5, 0.9),
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
)

##### b: ######################################################################
newplot2_bs <- ggplot2::ggplot(
  newplot2_dataAS |> tidytable::filter(
    Time > Start, Time < Stop
  ) |> tidytable::group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity, Subset
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
  # ) + ggplot2::geom_jitter(
  #   alpha = 0.25
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
  # ) + ggplot2::annotate(
  #   "text", x = c(1.5, 4.5), y = 5, label = c("Well\nAdapted", "Poorly\nAdapted")
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"), ordered = TRUE) ~ .
)

newplot2s <- ggpubr::ggarrange(
  plotlist = list(
    newplot2_as,
    newplot2_bs
  ), nrow = 1, widths = c(0.5, 0.4)
)

ggplot2::ggsave(plot = newplot2s, filename = "Figure2s1_Prototype2.png",
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
    Time > Start, Time < Stop, SpeciesAffinity == "50% 0, 50% 1"
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
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesAffinity ~ .
)
newplot3_b <- ggplot2::ggplot(
  newplot3_dataA |> tidytable::filter(
    Time > Start, Time < Stop, SpeciesAffinity == "Uniform(0, 1)"
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
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1)
) + ggplot2::labs(
  y = "Avg. Richness",
  x = "Habitat's Land-use"
) + ggplot2::guides(
  color = "none",
  fill = "none"
) + ggplot2::coord_cartesian(
  ylim = c(0, 42), expand = FALSE
) + ggplot2::facet_wrap(
  SpeciesAffinity ~ .
)

# Iteration 10 will have the actual affinities, rather than the affinitybins
# available, which will allow us to look at P/CDFs rather than bar charts.
# It still lets us think about how to quantify the distribution of affinities.
# I'm thinking Persistence as a weight, then by species aggregation. That way
# we get something like if I pick a random simulation, a random time, and then
# a random species, the plot shows the probability we would get a certain
# land-use preference out. (Weight by abundance as well for individuals, but
# that skews even more heavily towards basal species.)

newplot3_inset1 <- ggplot2::ggplot(
  newplot3_dataB |> tidytable::filter(SpeciesAffinity == "50% 0, 50% 1"),
  ggplot2::aes(
    x = AffinityBins,
    weight = Persistence,
    fill = Intervention
  )
) + ggplot2::geom_bar(
  show.legend = FALSE
) + ggplot2::facet_grid(
  . ~ Intervention
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_void(
) + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
) + ggplot2::coord_cartesian(
  expand = FALSE
)
newplot3_inset2 <- ggplot2::ggplot(
  newplot3_dataB |> tidytable::filter(SpeciesAffinity == "Uniform(0, 1)"),
  ggplot2::aes(
    x = AffinityBins,
    weight = Persistence,
    fill = Intervention
  )
) + ggplot2::geom_bar(
  show.legend = FALSE
) + ggplot2::facet_grid(
  . ~ Intervention
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat's Land-use"
) + ggplot2::theme_void(
) + ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
) + ggplot2::coord_cartesian(
  expand = FALSE
)

newplot3 <- ggpubr::ggarrange(
  plotlist = list(
    newplot3_a + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(newplot3_inset1),
      xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
    ),
    newplot3_b + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(newplot3_inset2),
      xmin = 0.55, xmax = 5.45, ymin = 30, ymax = 40
    )
  ), nrow = 1
)

ggplot2::ggsave(plot = newplot3, filename = "Figure3_Prototype6.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

##### b: ######################################################################
newplot3_bs <- ggplot2::ggplot(
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
  values = c("darkgreen", "burlywood4")
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::facet_grid(
  factor(SpeciesType, levels = c("Consumer", "Basal"), ordered = TRUE) ~
    SpeciesAffinity
) + ggplot2::labs(
  x = "Habitat's Land-use,\nsubdivided by Species Land-use Preference"
) + ggplot2::geom_text(
  data = rbind(
    data.frame(
      x = c(1:5 - 0.22, 1:5 + 0.22),
      y = 12000,
      lab = c(rep("0", 5), rep("1", 5)),
      SpeciesAffinity = "50% 0, 50% 1"
    ),
    data.frame(
      # Approximately the "right" spacing when blown up to large scales...
      x = c(1:5 - 0.36, 1:5 - 0.18, 1:5 - 0, 1:5 + 0.18, 1:5 + 0.36),
      y = 12000,
      lab = c(rep("0.1", 5), rep("0.3", 5), rep("0.5", 5),
              rep("0.7", 5), rep("0.9", 5)),
      SpeciesAffinity = "Uniform(0, 1)"
    )
  ),
  inherit.aes = FALSE,
  mapping = ggplot2::aes(
    x = x, y = y, label = lab
  )
)

# newplot3 <- ggpubr::ggarrange(
#   plotlist = list(
#     newplot3_a,
#     newplot3_b
#   ), nrow = 2, widths = c(0.5, 0.5)
# )

ggplot2::ggsave(plot = newplot3_bs, filename = "Figure3s1_Prototype2.png",
                units = "cm", width = 6.5*6, height = 6.5*2)

##### Why the losses: #########################################################
# Turnover amongst different groups?
# diversitiesTimeBC

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
  # SpeciesAffinity == "50% 0, 50% 1",
  # SpeciesAffinity == "Uniform(0, 1)",
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(endTimes |> dplyr::select(-Times))

newplot4_dataB <- diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  InterventionInitial == "(0.5)",
  SpeciesAffinity == "100% 0",
  !is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  Weight = 1, # for loess in geom_smooth
  .groups = "drop"
) |> tidytable::right_join(
  tidytable::expand(
    diversitiesRichness,
    tidytable::nesting(
      SpeciesAffinity, Intervention, # SpeciesAffinity not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
  ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesAffinity = ifelse(is.na(SpeciesAffinity), "100% 0", SpeciesAffinity)
) |> tidytable::filter(
  InterventionInitial == "(0.5)"
)

##### a: ######################################################################
newplot4_a <- plotMeanAndInner(
  rbind(
    newplot4_dataA |> tidytable::filter(
      Intervention %in% c("(0)", "(0.5)->(0)", "(0.5)", "(0.5)->(1)", "(1)")
    ),
    # We want to appear in the legend but not on the plot!
    newplot4_dataA |> tidytable::filter(
      PoolPatchSeed == newplot2_a_seed,
      Intervention %in% c("(0.5)->(0.25)", "(0.5)->(0.75)"),
      abs(Time - newplot2_a_time) == min(abs(Time - newplot2_a_time))
    ) |> tidytable::mutate(
      Value = -100
    )
  ),
  CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 7),
  fill = ggplot2::guide_legend(ncol = 7)
) + ggplot2::coord_cartesian(
  xlim = c(0, 30500), ylim = c(0, 42), expand = FALSE
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = 16300,
  xmax = 16500,
  ymin = 0, ymax = 42,
  fill = "grey",
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  # legend.position = c(0.35, 0.09),
  plot.tag.position = c(0.025, 1)
) + ggplot2::scale_x_continuous(
  breaks = (0:3)*10000
  # ) + ggforce::facet_zoom(
  #   xlim = c(16000, 17000),
  #   shrink = FALSE
)

##### Temporal Vs. Counterfactual Statistics: #################################
temporalVCounterfactualStats <- rbind(
  # Temporal Substitution
  diversitiesRichness |> tidytable::filter(
    NicheDistance == "5",
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) |> tidytable::left_join(endTimes |> dplyr::select(-Times)) |> tidytable::group_by(
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
  diversitiesRichness |> tidytable::filter(
    NicheDistance == "5",
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) |> tidytable::left_join(endTimes |> dplyr::select(-Times)) |> tidytable::group_by(
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
  # ) |> tidytable::filter(
  #   InterventionFinal == "(0)"
) |> tidytable::group_by(
  SpeciesAffinity:InterventionFinal
) |> tidytable::mutate(
  Deviation = Value - Value[Method=="Counterfactual"]
) |> tidytable::group_by(
  SpeciesAffinity, Intervention:InterventionFinal, Method
) |> tidytable::summarise(
  Mean = mean(Value),
  StDev = sd(Value),
  Bias = mean(Deviation),
  MeanAbsDev = mean(abs(Deviation)),
  PADGT1 = sum(abs(Deviation) > 1)/tidytable::n(),
  PADGT3 = sum(abs(Deviation) > 3)/tidytable::n(),
  PADGT5 = sum(abs(Deviation) > 5)/tidytable::n()
)

temporalVCounterfactualStats |> pivot_wider(
  names_from = Method, values_from = StDev,
  id_cols = c(SpeciesAffinity,
              Intervention, InterventionInitial, InterventionFinal)
) |> filter(
  Counterfactual != 0
) |> mutate(
  out = (Temporal - Counterfactual)
) |> pull(out) |> quantile(probs = (seq(from = 0.05, by = 0.05, to = 0.95)))

temporalVCounterfactualStats |> filter(
  Method != "Counterfactual"
) |> pull(Bias) |> summary()
temporalVCounterfactualStats |> filter(
  Method != "Counterfactual"
) |> pull(Bias) |> quantile(probs = (seq(from = 0.05, by = 0.05, to = 0.95)))

##### b: ######################################################################
# Could probably add in the counterfactual as well, but might be too messy?
newplot4_b <- ggplot2::ggplot(
  newplot4_dataB,
  aes(x = Time, y = Value,
      group = interaction(SpeciesAffinity, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  factor(Subset, levels = c("Consumer_0", "Basal_0"),
         labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::theme(
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

colorPalette_0.5 <- colorPalette[
  grepl(x = names(colorPalette), pattern = "^[(]0.5[])]")
  ]
newplot4_b_smooths <- ggplot2::ggplot_build(
  newplot4_b
)$data[[2]] |> dplyr::group_by(
  group
) |> dplyr::filter(
  x == max(x)
) |> dplyr::ungroup(
) |> dplyr::mutate(
  Subset = rev(levels(factor(newplot4_dataB$Subset)))[PANEL],
  Intervention = names(colorPalette_0.5)[
    apply(outer(colour, colorPalette_0.5, `==`), 1, which)
    ],
  yshift = y + c(-1, +2.5, -2, -4, -1, -3, -2, +2.5)
)

newplot4_b <-
  # Witchcraft from stackoverflow.com/a/6675163
  # works by pre-building the plot and then extracting coordinates.
  newplot4_b + ggplot2::coord_cartesian(
    xlim = c(0, 14), clip = "off"
  ) + ggplot2::geom_segment(
    data = newplot4_b_smooths,
    mapping = ggplot2::aes(x = x+1, y = yshift, xend = x, yend = y,
                           color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + ggplot2::geom_label(
    data = newplot4_b_smooths,
    mapping = ggplot2::aes(x = x+2.5, y = yshift,
                           label = Intervention, color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  )

newplot4 <- ggpubr::ggarrange(
  plotlist = list(
    newplot4_a,
    newplot4_b
  ), nrow = 1, common.legend = TRUE #, widths = c(0.5, 0.27, 0.23)
) + ggplot2::annotate(
  "curve", x = 0.29, y = 0.1, xend = 0.5, yend = 0.15,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
)

ggplot2::ggsave(plot = newplot4, filename = "Figure4_Prototype2.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

##### table: ##################################################################
diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  InterventionChange = abs(
    as.numeric(gsub(InterventionInitial, pattern = "[(]|[)]", replacement = ""))
    - as.numeric(gsub(InterventionFinal, pattern = "[(]|[)]", replacement = ""))
  ),
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  .groups = "drop"
) |> with(table(Intervention, sign(Value), SpeciesAffinity))

##### bs1: ####################################################################
newplot4_bs <- diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  InterventionInitial %in% c("(0)", "(1)"),
  SpeciesAffinity == "100% 0",
  !is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  Weight = 1, # for loess in geom_smooth
  .groups = "drop"
) |> tidytable::right_join(
  tidytable::expand(
    diversitiesRichness,
    tidytable::nesting(
      SpeciesAffinity, Intervention, # SpeciesAffinity not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
    # )
  ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesAffinity = ifelse(is.na(SpeciesAffinity), "100% 0", SpeciesAffinity)
) |> tidytable::filter(
  InterventionInitial %in% c("(0)", "(1)")
) |> ggplot2::ggplot(
  aes(x = Time, y = Value,
      group = interaction(SpeciesAffinity, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  # show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  # tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  strip.text.x = ggplot2::element_blank()
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  InterventionInitial +
    factor(Subset, levels = c("Consumer_0", "Basal_0"),
           labels = c("Consumer", "Basal"), ordered = TRUE) ~ .
) + ggplot2::theme(
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

colorPalette_01 <- colorPalette[
  grepl(x = names(colorPalette), pattern = "^[(](0|1)[])]")
  ]
newplot4_bs_smooths <- ggplot2::ggplot_build(
  newplot4_bs
)$data[[2]] |> dplyr::group_by(
  group
) |> dplyr::filter(
  x == max(x)
) |> dplyr::ungroup(
) |> dplyr::mutate(
  Subset = rev(levels(factor(newplot4_dataB$Subset)))[
    ((as.numeric(PANEL) - 1) %% 2) + 1
    ],
  Intervention = names(colorPalette_01)[
    apply(outer(colour, colorPalette_01, `==`), 1, which)
    ],
  InterventionInitial = substr(Intervention, 1, 3),
  yshift = y + c(+2.5, 0, -5, -10,
                 -15, -10, -5, +2,
                 +2, -5, -10, -15,
                 -6, -2, 0, 0)
)

newplot4_bs <-
  # Witchcraft from stackoverflow.com/a/6675163
  # works by pre-building the plot and then extracting coordinates.
  newplot4_bs + ggplot2::coord_cartesian(
    xlim = c(0, 13), clip = "off"
  ) + ggplot2::geom_segment(
    data = newplot4_bs_smooths,
    mapping = ggplot2::aes(x = x+1, y = yshift, xend = x, yend = y,
                           color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + ggplot2::geom_label(
    data = newplot4_bs_smooths,
    mapping = ggplot2::aes(x = x+2, y = yshift,
                           label = Intervention, color = Intervention),
    show.legend = FALSE,
    inherit.aes = FALSE
  )

ggplot2::ggsave(plot = newplot4_bs, filename = "Figure4s1_Prototype1.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

##### bs2: ####################################################################
newplot4_bs2 <- diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  # InterventionInitial %in% c("(0)", "(1)"),
  SpeciesAffinity == "100% 0",
  !is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  Weight = 1, # for loess in geom_smooth
  .groups = "drop"
) |> tidytable::right_join(
  tidytable::expand(
    diversitiesRichness,
    tidytable::nesting(
      SpeciesAffinity, Intervention, # SpeciesAffinity not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
  )
  # ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesAffinity =
    ifelse(is.na(SpeciesAffinity), "100% 0", SpeciesAffinity)
) |> tidytable::filter(
  # InterventionInitial %in% c("(0)", "(1)")
  Subset %in% c("Basal_0", "Consumer_0"),
  !is.na(InterventionInitial), !is.na(InterventionFinal),
  InterventionInitial != InterventionFinal
) |> ggplot2::ggplot(
  aes(x = Time, y = Value,
      group = interaction(SpeciesAffinity, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  # show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x",
  show.legend = FALSE
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  subtitle =
    "Columns = Final Land-use, Rows = Initial Land-use and Species Type",
  # tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  InterventionInitial +
    factor(Subset, levels = c("Consumer_0", "Basal_0"),
           labels = c("C0", "B0"),
           ordered = TRUE) ~ InterventionFinal
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  # strip.text.x = ggplot2::element_blank()
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

ggplot2::ggsave(plot = newplot4_bs2, filename = "Figure4s2_Prototype1.png",
                units = "cm", width = 6.5*5, height = 6.5*4)

##### bs3: ####################################################################
newplot4_bs3 <- diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  # InterventionInitial %in% c("(0)", "(1)"),
  SpeciesAffinity == "50% 0, 50% 1",
  !is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  Weight = 1, # for loess in geom_smooth
  .groups = "drop"
) |> tidytable::right_join(
  tidytable::expand(
    diversitiesRichness,
    tidytable::nesting(
      SpeciesAffinity, Intervention, # SpeciesAffinity not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
  )
  # ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesAffinity =
    ifelse(is.na(SpeciesAffinity), "50% 0, 50% 1", SpeciesAffinity)
) |> tidytable::filter(
  # InterventionInitial %in% c("(0)", "(1)")
  Subset %in% c("Basal_0", "Consumer_0",
                "Basal_1", "Consumer_1"),
  !is.na(InterventionInitial), !is.na(InterventionFinal),
  InterventionInitial != InterventionFinal
) |> ggplot2::ggplot(
  aes(x = Time, y = Value,
      group = interaction(SpeciesAffinity, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  # show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x",
  show.legend = FALSE
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  subtitle =
    "Columns = Final Land-use, Rows = Initial Land-use and Species Type",
  # tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  InterventionInitial +
    factor(Subset, levels = c("Consumer_0", "Basal_0",
                              "Consumer_1", "Basal_1"),
           labels = c("C0", "B0", "C1", "B1"),
           ordered = TRUE) ~ InterventionFinal
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  # strip.text.x = ggplot2::element_blank()
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

ggplot2::ggsave(plot = newplot4_bs3, filename = "Figure4s3_Prototype1.png",
                units = "cm", width = 6.5*5, height = 6.5*4)

##### bs4: ####################################################################
newplot4_dataBS4 <- diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  # InterventionInitial %in% c("(0)", "(1)"),
  SpeciesAffinity == "Uniform(0, 1)",
  !is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed,
  InterventionInitial, InterventionFinal, Subset
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  Weight = 1, # for loess in geom_smooth
  AffinityBins = gsub(pattern = "(Basal|Consumer)[_]", replacement = "",
                      x = Subset, perl = TRUE),
  SpeciesGuild = gsub(pattern = "(?=_).+", replacement = "",
                      x = Subset, perl = TRUE),
  .groups = "drop"
  #####
) |> unifyAffinityBins(
) |> tidytable::mutate(
  Subset = paste0(SpeciesGuild, "_", AffinityBins)
)

newplot4_bs4 <- newplot4_dataBS4 |> tidytable::right_join(
  tidytable::expand(
    newplot4_dataBS4,
    tidytable::nesting(
      SpeciesAffinity, Intervention, # SpeciesAffinity not working???
      InterventionInitial, InterventionFinal,
      Subset
    )
  )
  # ) |> tidytable::filter(Subset %in% c("Basal_0", "Consumer_0"))
) |> tidytable::mutate(
  Time = ifelse(is.na(Time), 0, Time),
  Value = ifelse(is.na(Value), 0, Value),
  Weight = ifelse(is.na(Weight), 1e9, Weight), # Unclear has an effect.
  SpeciesAffinity =
    ifelse(is.na(SpeciesAffinity), "Uniform(0, 1)", SpeciesAffinity)
) |> tidytable::filter(
  # InterventionInitial %in% c("(0)", "(1)")
  Subset %in% c("Basal_(0, 0.2]", "Consumer_(0, 0.2]",
                "Basal_(0.2, 0.4]", "Consumer_(0.2, 0.4]",
                "Basal_(0.4, 0.6]", "Consumer_(0.4, 0.6]",
                "Basal_(0.6, 0.8]", "Consumer_(0.6, 0.8]",
                "Basal_(0.8, 1]", "Consumer_(0.8, 1]"),
  !is.na(InterventionInitial), !is.na(InterventionFinal),
  InterventionInitial != InterventionFinal,
  InterventionInitial == "(0.5)"
) |> ggplot2::ggplot(
  aes(x = Time, y = Value,
      group = interaction(SpeciesAffinity, Intervention),
      # color = InterventionInitial
      # color = InterventionFinal
      color = Intervention
  )
) + ggplot2::geom_hline(
  yintercept = 0, color = "black"
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_smooth(
  # show.legend = FALSE,
  ggplot2::aes(weight = Weight),
  method = "loess",
  formula = "y~x",
  show.legend = FALSE
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  title = "Short Time Scales",
  subtitle =
    "Columns = Final Land-use, Rows = Initial Land-use and Species Type",
  # tag = "b)",
  x = "Time since Land-use Change",
  y = "Impact - Control (Richness)"
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::facet_grid(
  InterventionInitial +
    factor(Subset#,
           # levels = c("Consumer_0", "Basal_0",
           #                    "Consumer_1", "Basal_1"),
           # labels = c("C0", "B0", "C1", "B1"),
           # ordered = TRUE
           ) ~ InterventionFinal
) + ggplot2::theme(
  plot.tag.position = c(0.01, 1),
  # strip.text.x = ggplot2::element_blank()
  plot.background = ggplot2::element_rect(linetype = "solid")
  # panel.border = ggplot2::element_rect(linetype = "solid", fill = NA)
)

ggplot2::ggsave(plot = newplot4_bs4, filename = "Figure4s4_Prototype2.png",
                units = "cm", width = 6.5*5, height = 6.5*4)

### Plot 5: ###################################################################
# Progression of network change as we undergo intervention. As a base plot
# we use the richness changes of two experiments that are trading places.
# Both of these should probably be 100% 0. Maybe (0) -> (0.5) and (0.5) -> (0).
# Because of the size of the plots, we can't actually show them along the
# richness plots directly. We could potentially put labeled points instead,
# and then time staggered facets. I think we'll need higher resolution evals
# in order to capture the level of detail we're describing in the main text.
# We also need to convert the existing code for creating the networks into
# more general code, since we'll need to make a few here as well...

# Try something like the below to identify a good option:
# diversitiesRichness |> tidytable::filter(
#   NicheDistance == "5",
#   (PoolPatchSeed %in% as.character(383)),#:386)),
#   Metric == "Alpha Hill:0",
#   SpeciesAffinity == "100% 0",
#   Intervention %in% c("(0)->(0.5)", "(0.5)->(0)", "(0)", "(0.5)"),
#   Time > 16000, Time < 18000,
#   is.na(Subset)
# ) |> ggplot(
#   aes(x = Time, y = Value, color = Intervention)
# ) + ggplot2::geom_line(
# ) + ggplot2::facet_wrap(
#   PoolPatchSeed ~ .
# ) + coord_cartesian(
#   xlim = c(16600, 17300)
# )

newplot5_a_Specification <- rbind(diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(383)),#:386)),
  Metric == "Alpha Hill:0",
  SpeciesAffinity == "100% 0",
  Intervention %in% c("(0)", "(0.5)"),
  Time %in% c(16700),
  is.na(Subset)
), diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(383)),#:386)),
  Metric == "Alpha Hill:0",
  SpeciesAffinity == "100% 0",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)"),
  Time %in% c(16720, 16800, 16900, 17100),
  is.na(Subset)
))

newplot5_a_Networks <- generateNetworks(newplot5_a_Specification)

newplot5_a <- diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(383)),#:386)),
  Metric == "Alpha Hill:0",
  SpeciesAffinity == "100% 0",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)", "(0)", "(0.5)"),
  Time > 16000, Time < 18000,
  is.na(Subset)
) |> ggplot(
  aes(x = Time, y = Value, color = Intervention)
) + ggplot2::geom_line(
  # show.legend = FALSE
) + coord_cartesian(
  xlim = c(16600, 17300),
  ylim = c(6, 28)
) + ggplot2::geom_rect(
  data = newplot5_a_Specification |> tidytable::group_by(
    Time
  ) |> tidytable::summarise(
    xmin = Time - 2, xmax = Time + 2,
    ymin = 5, ymax = 27
  ),
  mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  color = "grey", inherit.aes = FALSE, alpha = 0.2
) + ggplot2::geom_point(
  show.legend = FALSE,
  data = newplot5_a_Specification
) + ggplot2::geom_vline(
  xintercept = (
    newplot5_a_Networks$Envs[[3]]$result$Ellipsis$Affinity$TimeIntervention
    / newplot5_a_Networks$Envs[[3]]$result$ReactionTime
  ), linetype = "dashed"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = ""#"Habitat's Land-use"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 4),
  fill = ggplot2::guide_legend(ncol = 4)
) + ggplot2::labs(
  tag = "a)"
  #   y = "Richness "
) + ggplot2::theme_minimal(
) + ggplot2::theme(
  legend.position = c(0.5, 1.2),
  plot.tag.position = c(0, 1),
  axis.title.y = ggplot2::element_blank()
)

# For some reason, this is returning a list of two plots rather than a single
# plot when used with ncol or nrow.
newplot5_b <- ggarrange(plotlist = ggarrange(
  plotlist = lapply(
    seq_along(newplot5_a_Networks$Envs),
    function(i, e) {
      e <- e[[i]]
      g <- e$singletonGraphs[[1]] + ggplot2::theme_void(
      ) + ggplot2::annotate(
        "text", x = -0.05, y = 0.5, label = e$Row$Time, hjust = 0
      ) + ggplot2::theme(
        panel.border = ggplot2::element_rect(
          # Ordinal in different order than color palette, and indexing
          # defaults to as.numeric for ordinal rather than as.character.
          color = colorPalette[as.character(e$Row$Intervention)],
          fill = NA
        ),
        plot.title = ggplot2::element_text(size = 13, hjust = 1)
      ) + ggplot2::ggtitle(
        if (i %in% c(1, 2, 6, 7)) {e$Row$Intervention} else {""}
      )

      return(g)
    },
    e = newplot5_a_Networks$Envs[c(1, 3:6, 2, 7:10)]
  ),
  ncol = 5
), nrow = 2, labels = list("b)", "c)"),
font.label = list(face = "plain"),
vjust = 1.4, hjust = 0)


newplot5 <- ggpubr::ggarrange(
  plotlist = list(
    newplot5_a,
    newplot5_b
  ), nrow = 2, heights = c(0.2/0.9, 0.7/0.9)
)

# Decorate with additional indicators.
newplot5 <-
  newplot5 + ggplot2::annotate(
  "curve", x = 0.2, y = 0.865, xend = 0.1, yend = 0.73,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
) + ggplot2::annotate(
  "curve", x = 0.23, y = 0.865, xend = 0.3, yend = 0.73,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
  curvature = -0.3
) + ggplot2::annotate(
  "curve", x = 0.33, y = 0.865, xend = 0.5, yend = 0.73,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
  curvature = -0.3
) + ggplot2::annotate(
  "curve", x = 0.45, y = 0.865, xend = 0.7, yend = 0.73,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
  curvature = -0.3
) + ggplot2::annotate(
  "curve", x = 0.7, y = 0.865, xend = 0.9, yend = 0.73,
  arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
  curvature = -0.2
) + ggplot2::annotate(
  "segment", linetype = "dashed", x = 0.215, y = 0.865, xend = 0.2, yend = 0.745
) + ggplot2::annotate(
  "segment", linetype = "dashed", x = 0.2, y = 0.745, xend = 0.2, yend = 0
)

ggplot2::ggsave(plot = newplot5, filename = "Figure5_Prototype2.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

##### 5s: #####################################################################
# Chris wanted to know if we could generalise figure 5 further. Given it is a
# one-off, I'm not entirely certain about how to do this, but I can draft some
# ideas. First, the figure 4 supplements characterise how the system changes in
# response to the varying land-uses and preferences. Chris is wanting to know
# then how the network structure is changing over time.

# One possibility to play with: biofabric, sorted by descending size within
# components.
# ggraph::ggraph(
#   layout = ggraph::create_layout(
#     graph = newplot5_a_Networks$Envs[[3]]$graphs[[1]] %N>% mutate(
#       nodesize = (log10(N)+5)/10+1
#     ), layout = "fabric", sort.by = rev(Size))
# ) + ggraph::geom_node_range(
#   mapping = aes(
#     color = Type,
#     size = nodesize
#   )
# ) + ggraph::geom_edge_span(
#   mapping = aes(
#     color = Type,
#     #color = node1.Type, # but then exploit+ between consumers is orange.
#     linetype = Type,
#     alpha = log10(effectNormalised)
#   ),
#   arrow = arrow(length = unit(2, 'mm'))
# ) + theme_minimal(
# ) + scale_color_manual(
#   values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
# ) + ggraph::scale_edge_color_manual(
#   values =
#     c("Exploit+" = "darkgreen", "Exploit-" = "burlywood4")
#   # c("Basal" = "darkgreen", "Consumer" = "burlywood4")
#
# )
# Not obviously effective, and the sorting doesn't work too well.
# So it's not quite plug-n-play.

##### 5s Idea 2: ##############################################################
# Another interesting option: we plot vertical kde's for each intervention at
# "-1", 10, 100, 200, and 400 time units from the intervention across all sims
# with that land-use/intervention and species affinity setup. Then, *inside the
# kde's* we plot all of points with little-no sizes along the axis (or maybe
# along the edge of the kde?). More imporantly, we plot the *edges* with a fixed
# alpha and size so that we can see the spread of interactions across the kdes.

####### Mockup 1: #############################################################
newplot5_as_Specification <- diversitiesRichness |> tidytable::filter(
  NicheDistance == "5",
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  # SpeciesAffinity == "100% 0",
  # SpeciesAffinity == "50% 0, 50% 1",
  SpeciesAffinity == "Uniform(0, 1)",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)"),
  is.na(Subset)
) |> tidytable::group_by(
  PoolPatch:InterventionFinal
) |> tidytable::filter(
  Time %% 1 != 0 | # I.e., time recorded just before the intervention!
    dplyr::lag(Time) %% 1 != 0 | # First time after the intervention (v9 !!)
    dplyr::lag(Time, n = 11) %% 1 != 0 | # 10 * 10 + 1 time steps
    dplyr::lag(Time, n = 21) %% 1 != 0 |
    dplyr::lag(Time, n = 41) %% 1 != 0
) |> tidytable::mutate(
  Time2 = tidytable::case_when(
    Time %% 1 != 0 ~ -1,
    dplyr::lag(Time) %% 1 != 0 ~ 10,
    dplyr::lag(Time, n = 2) %% 1 != 0 ~ 100,
    dplyr::lag(Time, n = 3) %% 1 != 0 ~ 200,
    dplyr::lag(Time, n = 4) %% 1 != 0 ~ 400
  )
) |> tidytable::ungroup(
)

newplot5_as_Networks <- generateNetworks(newplot5_as_Specification)

newplot5_kdes <- lapply(
  newplot5_as_Networks$Envs, function(e) {
    e$trophics$EdgeVertexLists[[1]][[1]]$Vertices |> tidytable::select(
      node, Type, Size, N
    ) |> cbind(
      e$Row |> tidytable::select(Time, Time2, PoolPatch:InterventionFinal)
      ) |> tidytable::mutate(
        AffinityVals = e$result$Ellipsis$Affinity$SpeciesAffinities[
          as.numeric(substring(node, 2))
          ]
      )
  }
) |> tidytable::bind_rows(
# ) |> tidytable::group_by(
#   PoolPatch:InterventionFinal
# ) |> tidytable::arrange(
#   Time
# ) |> tidytable::mutate( # fix for having not done it ahead of time...
#   Time2 = tidytable::case_when(
#     Time == unique(Time)[1] ~ -1,
#     Time == unique(Time)[2] ~ 10,
#     Time == unique(Time)[3] ~ 100,
#     Time == unique(Time)[4] ~ 200,
#     Time == unique(Time)[5] ~ 400
#   )
# ) |> tidytable::ungroup(
)

# newplot5_graph <- lapply(newplot5_as_Networks$Envs, function(e)
#   e$graphs[[1]] %N>% select(
#     node, Type, Size
#   ) %N>% mutate(
#     e$Row |> select(Time, PoolPatch:InterventionFinal)
#   ) %E>% mutate(
#     e$Row |> select(Time, PoolPatch:InterventionFinal)
#   )
# ) |> tidygraph::bind_graphs(
# ) %N>% tidygraph::group_by( # Only supports manual specification
#   PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
#   EventsSeed, InitialConditions, InitialConditionsSeed, Dispersal,
#   NicheDistance, Affinity, AffinitySeed, InterventionPatchType,
#   InterventionPatchSeed, InterventionTimeType, InterventionTimeSeed,
#   InterventionDispersal, InterventionNicheDistance, Intervention,
#   SpeciesAffinity, InterventionInitial, InterventionFinal
# ) %N>% tidygraph::arrange(
#   Time
# ) %N>% tidygraph::mutate( # fix for having not done it ahead of time...
#   Time2 = dplyr::case_when(
#     Time == unique(Time)[1] ~ -1,
#     Time == unique(Time)[2] ~ 10,
#     Time == unique(Time)[3] ~ 100,
#     Time == unique(Time)[4] ~ 200,
#     Time == unique(Time)[5] ~ 400
#   )
# ) %N>% tidygraph::ungroup(
# ) %E>% tidygraph::group_by( # Only supports manual specification
#   PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
#   EventsSeed, InitialConditions, InitialConditionsSeed, Dispersal,
#   NicheDistance, Affinity, AffinitySeed, InterventionPatchType,
#   InterventionPatchSeed, InterventionTimeType, InterventionTimeSeed,
#   InterventionDispersal, InterventionNicheDistance, Intervention,
#   SpeciesAffinity, InterventionInitial, InterventionFinal
# ) %E>% tidygraph::arrange(
#   Time
# ) %E>% tidygraph::mutate( # fix for having not done it ahead of time...
#   Time2 = dplyr::case_when(
#     Time == unique(Time)[1] ~ -1,
#     Time == unique(Time)[2] ~ 10,
#     Time == unique(Time)[3] ~ 100,
#     Time == unique(Time)[4] ~ 200,
#     Time == unique(Time)[5] ~ 400
#   )
# ) %E>% tidygraph::ungroup(
# )

# ggraph::ggraph(
#   ggraph::create_layout(newplot5_graph, "manual",
#                         y = newplot5_graph %N>% tidygraph::pull(Size), x = -0.5)
ggplot2::ggplot(
) + ggplot2::geom_density(
  data = newplot5_kdes,
  mapping = ggplot2::aes(
    y = Size, fill = Type, color = Intervention
  ),
  trim = TRUE
  # ) + ggplot2::geom_rug(
  #   ggplot2::aes(color = Type)
# ) + ggraph::geom_edge_arc(
#   mapping = aes(
#     color = Type,
#     #color = node1.Type, # but then exploit+ between consumers is orange.
#     linetype = Type,
#     alpha = log10(effectNormalised),
#     end_cap = circle(2, 'pt')
#   ),
#   arrow = arrow(length = unit(2, 'mm')),
#   alpha = 0.2,
#   show.legend = FALSE
# ) + ggraph::geom_node_point(
#   mapping = aes(
#     color = Type
#   ),
#   show.legend = FALSE
#   # ) + ggplot2::geom_hline(
#   #   yintercept = -1, linetype = "dashed", color = "black"
) + ggplot2::facet_grid(
  Intervention + SpeciesAffinity ~ Time2
) + ggplot2::scale_y_log10(
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
) + ggplot2::theme_minimal(
)

ggplot2::ggsave(
  ggplot2::ggplot(
  ) + ggplot2::geom_density(
    data = newplot5_kdes,
    mapping = ggplot2::aes(
      y = Size, fill = Type, color = Intervention
    ),
    trim = TRUE
  ) + ggplot2::facet_grid(
    Intervention + SpeciesAffinity ~ Time2
  ) + ggplot2::scale_y_log10(
  ) + ggplot2::scale_color_manual(
    values = colorPalette,
    name = "Habitat Land-use"
  ) + ggplot2::scale_fill_manual(
    values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
  ) + ggplot2::theme_minimal(
  ),
  # filename = "Figure5s1_Prototype1.png", # 100% 0
  # filename = "Figure5s2_Prototype1.png", # 50% 0, 50% 1
  filename = "Figure5s3_Prototype1.png", # Uniform(0, 1)
  units = "cm", width = 6.5*3, height = 6.5*2
)

ggplot2::ggsave(
  ggplot2::ggplot(
  ) + ggplot2::geom_density_2d(
  # ) + ggplot2::geom_bin_2d(
    data = newplot5_kdes,
    mapping = ggplot2::aes(
      x = AffinityVals, y = Size,
      # fill = Type,
      color = Intervention,
      group = interaction(Type, Intervention)
    ),
    # bins = 10
    alpha = 0.4,
    contour_var = "count",
    adjust = 0.7
    # trim = TRUE
  ) + ggplot2::geom_hline(
    yintercept = 0.1, color = "red", show.legend = FALSE
  ) + ggplot2::facet_grid(
    Intervention + SpeciesAffinity ~ Time2
  ) + ggplot2::scale_y_log10(
  ) + ggplot2::scale_color_manual(
    values = colorPalette,
    name = "Habitat Land-use"
  ) + ggplot2::scale_fill_manual(
    values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
  # ) + ggplot2::scale_fill_viridis_c(
  ) + ggplot2::theme_minimal(
  ) + ggplot2::xlab(
    "Land-use Type"
  ),
  filename = "Figure5s4_Prototype1.png", # Uniform(0, 1)
  units = "cm", width = 10*3, height = 10*2
)

### Plot 6: ###################################################################
# Pseudo-relaxation time of the system from the intervention to its new final
# state, characterised as the difference between counterfactual always in final
# state and the intervention to the final state.
# Because of the rescalings, if we index it by the original time of intervention
# and then compare, we should be able to observe roughly the same relaxation
# time (so long as we don't perform the second rescaling to have interventions
# fixed at 0.5?).
