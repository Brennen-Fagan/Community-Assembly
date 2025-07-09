# Plotting with both ColExt and Divs
load("ColExt9a9_flat.RData")
load("diversitiesFlattened9a9_subset4.RData")

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
source("flattenDiversity.R")
source("CalculateTrophicStructure.R") # Calculator creator.
source("toCheddar.R") # Updated function.

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
      levels = c(
        "(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)->(0)", "(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)", "(0.75)->(1)",
        "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)"
      ), ordered = TRUE
    )
  )
}

### Plotting: #################################################################
plotMeanAndInner <- function(
  data, CIs = c(0.5, 0.95),
  facets = as.formula(Intervention ~ SpeciesAffinity)
) {
  data$Subset <- ifelse(is.na(data$Subset), "NA", data$Subset)

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
  for (CI in CIs) {
    baseplot <- baseplot + ggplot2::geom_ribbon(
      data = data %>% tidytable::mutate(
        Time = round(Time, digits = -2)
      ) %>% tidytable::group_by(
        Time, Subset,
        Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
      ) %>% tidytable::summarise(
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
  baseplot <- (baseplot
               # ) + ggplot2::geom_smooth(
               #   se = FALSE,
               #   # method = "loess", span = 0.5
               #   method = "gam", formula =
               #     y ~ s(x, bs = "tp", k = 20) + s(PoolPatchSeed, bs = "re")
  ) + ggplot2::geom_line(
    data = data %>% tidytable::mutate(
      Time = round(Time, digits = -2)
    ) %>% tidytable::group_by(
      Time, Subset,
      Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
    ) %>% tidytable::summarise(
      Value = mean(Value, na.rm = TRUE)
    )
  ) + ggplot2::facet_grid(
    facets
  ) + ggplot2::scale_color_manual(
    values = colorPalette, aesthetics = c("color", "fill"),
    name = "Habitat Land-use"
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

plotValueChart <- function(
  data, facets = as.formula(Intervention ~ SpeciesAffinity)
) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = interaction(InType,
                                 OutType,
                                 sep = "\n"),#"",
                 y = ChartValue,
                 fill = interaction(InType,
                                    OutType,
                                    sep = "/"))
  ) + ggplot2::geom_col(
    show.legend = FALSE
  ) + ggplot2::facet_grid(
    facets
  ) + ggplot2::theme_minimal(
  ) + ggplot2::labs(
    x = "Enter/Exit Type", y = "Count"#, fill = "In/Out"
  )
}

# Note: unlike the other plots, this one is called inside the environments.
#       The other plots are called inside the createOverviewFigure function.
plotGraph <- function(graph, mainLayout, legends = FALSE) {
  obj <- ggraph::ggraph(
    graph = graph %N>% mutate(
      nodesize = (log10(N)+5)/10+1
    ),
    layout = graph %N>% data.frame(
    ) %>% select(node) %>% left_join(
      mainLayout %>% select(x, y, node)
    )
  ) + ggraph::geom_edge_diagonal(
    mapping = aes(
      color = Type,
      linetype = Type,
      alpha = log10(effectNormalised),
      end_cap = circle(node2.nodesize+2, 'pt')
    ),
    # linewidth = 2,
    arrow = arrow(length = unit(2, 'mm'))#,
    # start_cap = circle(1, 'mm'),
    # end_cap = circle(5, 'mm')
  ) + ggraph::geom_node_point(
    mapping = aes(
      color = Type,
      size = nodesize
    )
  ) + ggplot2::geom_hline(
    yintercept = -1, linetype = "dashed", color = "black"
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

  if (!legends) {
    obj <- obj + ggplot2::theme(
      legend.position = "none"
    )
  }

  return(obj)
}

# Primary figure generator for Figures 2 and 3 currently.
createOverviewFigure <- function(
  .divs, # diversitiesAll
  .ps, # Pers
  .ces, # ColExt
  .ets, # endTimes
  .gs, # Singleton Graphs
  ... # commonfilters
) {
  ggpubr::ggarrange(
    # Richness
    plotMeanAndInner(
      .divs %>% tidytable::filter(
        ...,
        Metric == "Alpha Hill:0",
        is.na(Subset)
      ),
      CIs = 0.75,
      facets = as.formula(. ~ .)
    ) + ggplot2::labs(
      y = "Richness"
    ) + ggplot2::guides(
      linetype = "none",
      color = ggplot2::guide_legend(ncol = 3),
      fill = ggplot2::guide_legend(ncol = 3)
    ) + ggplot2::theme(
      # legend.position = "bottom",
      # legend.direction = "horizontal"
      legend.position = c(0.4, 0.09)
    ) + ggplot2::geom_rect(
      data = data.frame(
        1 # 1 rectangle per row, so dummy df to prevent overplotting
      ),
      xmin = min(.ets$Start),
      xmax = max(.ets$Stop),
      ymin = 0, ymax = 45,
      fill = "grey",
      alpha = 0.2,
      inherit.aes = FALSE
    ),
    # plotValueChart(
    #   rbind(
    #     .ps %>% tidytable::filter(
    #       ...,
    #       Persistence > 0,
    #       InType != externalNames["Dispersal"],
    #       In < Stop, Out > Start
    #     ) %>% tidytable::group_by(
    #       SpeciesAffinity, InType, OutType, Intervention
    #     ) %>% tidytable::summarise(
    #       ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
    #     ) %>% dplyr::mutate( # Tidytable renders as character again!
    #       Intervention =
    #         factor(Intervention,
    #                levels = rev(c("(0)", "(0)->(0.5)", "(0)->(1)",
    #                               "(0.5)",
    #                               "(1)")),
    #                ordered = TRUE)
    #     ),
    #     .ces %>% tidytable::filter(
    #       ...,
    #       !Success | EventType == "Present",
    #       Times > Start, Times < Stop
    #     ) %>% tidytable::mutate(
    #       InType = externalNames[
    #         ifelse(EventType == "Arrival", "Failed Arrival", "Present")
    #         ],
    #       OutType = externalNames["NA"]
    #     ) %>% tidytable::group_by(
    #       SpeciesAffinity, InType, OutType, Intervention
    #     ) %>% tidytable::summarise(
    #       ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
    #     ) %>% dplyr::mutate( # Tidytable renders as character again!
    #       Intervention =
    #         factor(Intervention,
    #                levels = rev(c("(0)", "(0)->(0.5)", "(0)->(1)",
    #                               "(0.5)",
    #                               "(1)")),
    #                ordered = TRUE)
    #     )
    #   ),
    #   facets = as.formula(Intervention ~ .)
    # ) + ggplot2::theme(
    #   legend.position = "bottom"
    # ) + ggplot2::guides(
    #   fill = ggplot2::guide_legend(nrow = 2)
    # ) + ggplot2::ylab(
    #   "Average Count in a Simulation"
    # ),
    ggpubr::ggpar(
      ggpubr::ggarrange(
        plotlist = lapply(seq_along(.gs), function(i) {
          g <- .gs[[i]]
          if (i != length(.gs)) {
            g <- g + ggplot2::labs(
              x = NULL
            )
          }
          if (i != round(length(.gs)/2)) {
            g <- g + ggplot2::labs(
              y = ""
            )
            # } else if (!length(.gs) %% 2) { # Even
            #   g <- g + ggplot2::theme(
            #     axis.title.y.left =
            #   )
          }
          return(g)
        }),
        ncol = 1
      ),
      legend = "none"),
    ggplot2::ggplot(
      .ps %>% tidytable::filter(
        ...,
        Persistence > 0,
        In < Stop, Out > Start
      ),
      ggplot2::aes(
        x = Size,
        # color = AffinityBins,
        y = Persistence
      )
    ) + ggplot2::geom_hex(
    ) + ggplot2::theme_minimal(
    ) + ggplot2::scale_x_log10(
    ) + ggplot2::scale_y_log10(
    ) + ggplot2::scale_fill_viridis_c(
      trans = "log10", direction = -1
    ) + ggplot2::theme(
      legend.position = "bottom"
    ) + ggplot2::facet_grid(
      factor(Intervention, levels = rev(levels(Intervention))) ~ .
    ) + ggplot2::geom_vline(
      xintercept = 0.10, color = "red"
    ),
    nrow = 1, widths = c(1, 0.6, 0.6), #common.legend = TRUE
    labels = "auto"
  )
}

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

# Diversities: #################################################################
interventionMatrix <- matrix(
  c("(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
    "(0.25)->(0)", "(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
    "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
    "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)", "(0.75)->(1)",
    "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)"),
  byrow = TRUE, nrow = 5)

diversitiesAll <- diversitiesAll %>% changeAffinityLevels(
) %>% changeInterventionLevels(
) %>% tidytable::mutate(
  #TODO consider moving this into changeInterventionLevels.
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

# Colonisation/Extirpation: ####################################################
ColExt <- ColExt %>% tidytable::rename(
  DispersalParam = Dispersal
) %>% tidytable::mutate(
  SpeciesAffinity =
    affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
) %>% changeAffinityLevels(
) %>% tidytable::left_join(
  diversitiesInterventionStrings,
  by = c("Affinity", "PoolPatch", "InterventionPatchType"),
  multiple = "all"
) %>% changeInterventionLevels(
) %>% tidytable::left_join(
  endTimes %>% tidytable::select(-Times),
  by = c(
    "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed", "Events",
    "EventsSeed", "InitialConditions", "InitialConditionsSeed",
    "DispersalParam", "NicheDistance", "Affinity", "AffinitySeed"
  )
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

# Define Species Color Scale from AffinityBins: ################################
# Can't be done before hand because the affinity bins are calculated differently
# for each of the runs.
SpeciesPaletteDF <- data.frame(
  Bins = sort(levels(Pers$AffinityBins))
) %>% tidytable::separate_wider_delim(
  Bins, delim = c(","),
  cols_remove = FALSE
) %>% tidytable::mutate(
  Bins1 = gsub(pattern = "(", replacement = "", Bins1, fixed = TRUE),
  Bins2 = gsub(pattern = "]", replacement = "", Bins2, fixed = TRUE),
  Bins2 = ifelse(is.na(Bins2), Bins1, Bins2),
  Bins1 = as.numeric(Bins1),
  Bins2 = as.numeric(Bins2),
  ColorNumber = (Bins1 + Bins2)/2,
  color =
    # Define the Color Palette to reasonable resolution
    colorRampPalette(colorPalette[c("(0)", "(0.5)", "(1)")])(length(Bins))[
      # Apply Color Palette judiciously, stackoverflow.com/a/35656396
      findInterval(ColorNumber, seq(0, 1, length.out = length(Bins)))
      ]
) %>% tidytable::arrange(
  Bins
)

SpeciesPalette <- c(SpeciesPaletteDF$color)
names(SpeciesPalette) <- SpeciesPaletteDF$Bins

# Plots: #####################################################################
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

plot2 <- createOverviewFigure(
  .divs = diversitiesAll,
  .ps = Pers,
  .ces = ColExt,
  .ets = endTimes,
  .gs = rev(lapply(
    targetEnvs[targetEnvsIndex %>% dplyr::filter(
      SpeciesAffinity == "rep_0",
      NicheDistance == "5",
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ) %>% dplyr::pull(ID)
    ],
    function(env) env$singletonGraphs[[1]][[1]]
  )),
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)

ggplot2::ggsave(
  plot2,
  filename = "Figure2_Prototype5.png",
  units = "px", height = 3200, width = 4800 # Words are too small
  # height = 1600, width = 2400 # Words over-run
)

plot3 <- createOverviewFigure(
  .divs = diversitiesAll,
  .ps = Pers,
  .ces = ColExt,
  .ets = endTimes,
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)"),
  !(PoolPatchSeed %in% c("341", "342"))
)

ggplot2::ggsave(
  plot3,
  filename = "Figure3_Prototype3.png",
  units = "px", height = 3200, width = 4800 # Words are too small
  # height = 1600, width = 2400 # Words over-run
)


# Should make proper by using a left_join with endtimes
# This is also an expensive operation to do on everything so should be
# filtered down to a small subset if possible.
diversitiesAll %>% dplyr::filter(
  Time > 20000, Time < 30000,
  Metric == "Alpha Hill:0", !is.na(Subset), NicheDistance == 5
) %>% tidyr::separate(
  Subset, into = c("Type", "Preference"), sep = "_"
) %>% group_by(
  Time, Environment1, Environment2, Metric, Type,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal
) %>% summarise(
  # Collapse over preferences
  Value = sum(ifelse(is.na(Value), 0, Value)),
  .groups = "drop"
) %>% dplyr::group_by(
  Environment1, Environment2, Metric, Type,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal
) %>% dplyr::summarise(
  q0025 = quantile(Value, 0.025, na.rm = TRUE),
  q025 = quantile(Value, 0.25, na.rm = TRUE),
  q050 = median(Value, na.rm = TRUE),
  q075 = quantile(Value, 0.75, na.rm = TRUE),
  q0975 = quantile(Value, 0.975, na.rm = TRUE)
) %>% tidyr::pivot_wider(
  names_from = Type, values_from = q0025:q0975
) -> temp
plot4 <- temp %>% ggplot2::ggplot(
  ggplot2::aes(x = q050_Basal, y = q050_Consumer, color = Intervention)
) + ggplot2::geom_point(alpha = 5/44) + ggplot2::facet_wrap(
  .~SpeciesAffinity
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::stat_ellipse() + ggplot2::theme_minimal()



# New Plot 2:##################################################################
newplot2_filtration <- function(.) {tidytable::filter(
  .,
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)}

newplot2_dataA <- diversitiesAll %>% newplot2_filtration(
) %>% tidytable::filter(
  Metric == "Alpha Hill:0",
  is.na(Subset)
) %>% tidytable::left_join(endTimes %>% dplyr::select(-Times))

newplot2_dataB <- diversitiesAll %>% newplot2_filtration(
) %>% tidytable::left_join(
  endTimes %>% dplyr::select(-Times)
) %>% tidytable::filter(
  Time > Start, Time < Stop,
  Metric == "Alpha Abundance",
  !is.na(Subset)
) %>% tidytable::separate_wider_delim(
  cols = Subset, names = c("Guild", "AffinityBin"), delim = "_"
) %>% tidytable::group_by(
  Environment1, Environment2, Metric, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal, Guild, Time
) %>% tidytable::summarise(# Across Affinity Bins
  Value = sum(Value, na.rm = TRUE), .groups = "drop_last"
) %>% tidytable::summarise(# Across Time
  Mean = mean(Value),
  q025 = quantile(Value, probs = 0.25),
  q075 = quantile(Value, probs = 0.75)
) %>% tidytable::pivot_wider(
  names_from = Guild, values_from = c(Mean, q025, q075)
)

newplot2_dataC <- Pers %>% newplot2_filtration(
) %>% tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) %>% tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) %>% tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, AffinityBins,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, DispersalParam, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) %>% tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
)

newplot2_dataD <- diversitiesAll %>% newplot2_filtration(
) %>% tidytable::left_join(
  endTimes %>% dplyr::select(-Times)
) %>% tidytable::filter(
  Time > Start, Time < Stop,
  Metric == "TimeJaccard",
  is.na(Subset)
) %>% tidytable::group_by(
  Environment1, Environment2, Metric, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal, Start, Stop
) %>% tidytable::summarise(# Across Time
  Value = mean(Value, na.rm = TRUE), .groups = "drop"
)

newplot2_a <- plotMeanAndInner(
  newplot2_dataA, CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = newplot2_dataA %>% tidytable::filter(
    PoolPatchSeed == targetSeed,
    abs(Time - targetTimes) == min(abs(Time - targetTimes))
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 3),
  fill = ggplot2::guide_legend(ncol = 3)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), expand = FALSE
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot2_filtration %>% dplyr::pull(ID))[1]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 7, ymax = 17
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot2_filtration %>% dplyr::pull(ID))[2]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 18, ymax = 28
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot2_filtration %>% dplyr::pull(ID))[3]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
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
  legend.position = c(0.25, 0.09),
  plot.tag.position = c(0.05, 1)
)

newplot2_b <- ggplot2::ggplot(
  newplot2_dataB,
  ggplot2::aes(
    x = Mean_Basal,
    y = Mean_Consumer,
    color = Intervention)
) + ggplot2::geom_point(
  show.legend = FALSE
  # ) + ggplot2::geom_errorbar(
  #   inherit.aes = FALSE,
  #   mapping = ggplot2::aes(
  #     x = Mean_Basal, ymin = q025_Consumer, ymax = q075_Consumer,
  #     color = Intervention
  #   )
  # ) + ggplot2::geom_errorbarh(
  #   inherit.aes = FALSE,
  #   mapping = ggplot2::aes(
  #     y = Mean_Consumer, xmin = q025_Basal, xmax = q075_Basal,
  #     color = Intervention
  #   )
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::coord_fixed(
) + ggplot2::labs(
  x = "Basal\nAbundance",
  y = "Consumer Abundance"
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
)

# Size as a function of Total amount of time spent in simulation
# Color is intervention type.
newplot2_c <- ggplot2::ggplot(
  newplot2_dataC,
  ggplot2::aes(
    # x = Persistence,
    y = Size,
    x = Intervention,
    # color = Intervention,
    fill = SpeciesType
  )
) + ggplot2::geom_hline(
  yintercept = 0.10, linetype = "dashed"
) + ggplot2::geom_violin(
  position = "identity",
  scale = "count", adjust = 1/10,
  show.legend = FALSE
) + ggplot2::scale_y_log10(
  # ) + ggplot2::scale_color_manual(
  #   values = colorPalette, #aesthetics = c("color", "fill"),
  #   name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("limegreen", "goldenrod2")
) + ggplot2::theme_minimal(
  # ) + ggplot2::guides(
  #   color = "none",
  #   fill = ggplot2::guide_legend(
  #     title = ggplot2::element_blank(),
  #     reverse = TRUE
  #   )
  # ) + ggplot2::theme(
  #   legend.position = c(.84, .93)
) + ggplot2::labs(
  tag = "c)", x = "Habitat"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
)

# Average Temporal Beta as a function of Intervention Type.
newplot2_d <- ggplot2::ggplot(
  newplot2_dataD,
  ggplot2::aes(
    x = Intervention,
    y = Value,
    fill = Intervention
  )
) + ggplot2::geom_boxplot(
  notch = TRUE,
  # Double-checking, but no outliers it seems
  outlier.colour = "black", outlier.size = 200,
  show.legend = FALSE
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  tag = "d)", x = "Habitat", y = "Time-Jaccard"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
)

newplot2 <- ggpubr::ggarrange(
  plotlist = list(
    newplot2_a,
    newplot2_b,
    ggpubr::ggarrange(newplot2_c, newplot2_d, ncol = 1)
  ), nrow = 1, widths = c(0.55, 0.25, 0.2)
)
ggplot2::ggsave(plot = newplot2, filename = "aplot2.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

# New Plot 3:##################################################################
newplot3_filtration <- function(.) {tidytable::filter(
  .,
  SpeciesAffinity == "50% 0, 50% 1",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)}

newplot3_dataA <- diversitiesAll %>% newplot3_filtration(
) %>% tidytable::filter(
  Metric == "Alpha Hill:0",
  is.na(Subset)
) %>% tidytable::left_join(endTimes %>% dplyr::select(-Times))

newplot3_dataB <- diversitiesAll %>% newplot3_filtration(
) %>% tidytable::left_join(
  endTimes %>% dplyr::select(-Times)
) %>% tidytable::filter(
  Time > Start, Time < Stop,
  Metric == "Alpha Biomass",
  !is.na(Subset)
) %>% tidytable::separate_wider_delim(
  cols = Subset, names = c("Guild", "AffinityBin"), delim = "_"
) %>% tidytable::group_by(
  Environment1, Environment2, Metric, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal, Guild, AffinityBin
) %>% tidytable::summarise(# Across Time
  Mean = mean(Value)
) %>% tidytable::pivot_wider(
  names_from = "Guild", values_from = "Mean"
)

newplot3_dataC <- Pers %>% newplot3_filtration(
) %>% tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) %>% tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) %>% tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, AffinityBins,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, DispersalParam, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) %>% tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
) %>% tidytable::group_by(
  Environment, SpeciesType, AffinityBins,
  PoolPatch, Interactions, Events,
  InitialConditions, DispersalParam, NicheDistance,
  Affinity, InterventionPatchType,
  InterventionTimeType, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) %>% tidytable::summarise( # Sum over Appearances.
  Persistence = mean(Persistence),
  Size = 10^(mean(log10(Size))),
  .groups = "drop"
)

newplot3_a <- plotMeanAndInner(
  newplot3_dataA, CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = newplot3_dataA %>% tidytable::filter(
    PoolPatchSeed == targetSeed,
    abs(Time - targetTimes) == min(abs(Time - targetTimes))
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 3),
  fill = ggplot2::guide_legend(ncol = 3)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), expand = FALSE
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot3_filtration %>% dplyr::pull(ID))[1]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 21, ymax = 30
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot3_filtration %>% dplyr::pull(ID))[2]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 11, ymax = 20
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot3_filtration %>% dplyr::pull(ID))[3]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 1, ymax = 10
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(newplot3_dataA$Start),
  xmax = max(newplot3_dataA$Stop),
  ymin = 0, ymax = max(newplot3_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.25, 0.09),
  plot.tag.position = c(0.02, 1)
)

# Trying to support a breakdown of richness and abundance to show what is
# happening to, and because of, the species that have the wrong land-use
# preference for the given habitat.

newplot3_b <- ggplot2::ggplot(
  newplot3_dataB ,
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = Intervention,
    shape = AffinityBin
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::coord_fixed(
) + ggplot2::labs(
  x = "Basal Abundance",
  y = "Consumer Abundance"
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
)

newplot3_c <- ggplot2::ggplot(
  newplot3_dataC,
  ggplot2::aes(
    x = Persistence,
    y = Size,
    shape = AffinityBins,
    color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_hline(
  yintercept = 10^-1,
  linetype = "dashed"
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "c)"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
)

newplot3 <- ggpubr::ggarrange(
  plotlist = list(
    newplot3_a,
    ggpubr::ggarrange(newplot3_b, newplot3_c, ncol = 1)
  ), nrow = 1, widths = c(0.55, 0.45)
)
ggplot2::ggsave(plot = newplot3, filename = "aplot3.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

# New Plot 4: #################################################################
newplot4_filtration <- function(.) {tidytable::filter(
  .,
  SpeciesAffinity == "Uniform(0, 1)",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)}

newplot4_dataA <- diversitiesAll %>% newplot4_filtration(
) %>% tidytable::filter(
  Metric == "Alpha Hill:0",
  is.na(Subset)
) %>% tidytable::left_join(endTimes %>% dplyr::select(-Times))

newplot4_dataB <- diversitiesAll %>% newplot4_filtration(
) %>% tidytable::left_join(
  endTimes %>% dplyr::select(-Times)
) %>% tidytable::filter(
  Time > Start, Time < Stop,
  Metric == "Alpha Biomass",
  !is.na(Subset)
) %>% tidytable::separate_wider_delim(
  cols = Subset, names = c("Guild", "AffinityBin"), delim = "_"
) %>% tidytable::group_by(
  Environment1, Environment2, Metric, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal, Guild, AffinityBin
) %>% tidytable::summarise(# Across Time
  Mean = mean(Value)
) %>% tidytable::pivot_wider(
  names_from = "Guild", values_from = "Mean"
) %>% tidytable::separate(
  col = "AffinityBin", into = c("Left", "Right"), sep = ","
) %>% tidytable::mutate(
  Left = round(as.numeric(gsub(pattern = "^.", replacement = "", x = Left))*5)/5,
  Right = round(as.numeric(gsub(pattern = ".$", replacement = "", x = Right))*5)/5,
  AffinityBin = paste0("(", Left, ", ", Right, "]")
)

newplot4_dataC <- Pers %>% newplot4_filtration(
) %>% tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) %>% tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) %>% tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, AffinityBins,
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, DispersalParam, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) %>% tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
) %>% tidytable::group_by(
  Environment, SpeciesType, AffinityBins,
  PoolPatch, Interactions, Events,
  InitialConditions, DispersalParam, NicheDistance,
  Affinity, InterventionPatchType,
  InterventionTimeType, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
) %>% tidytable::summarise( # Sum over Appearances.
  Persistence = mean(Persistence),
  Size = 10^(mean(log10(Size))),
  .groups = "drop"
) %>% tidytable::separate(
  col = "AffinityBins", into = c("Left", "Right"), sep = ","
) %>% tidytable::mutate(
  Left = round(as.numeric(gsub(pattern = "^.", replacement = "", x = Left))*5)/5,
  Right = round(as.numeric(gsub(pattern = ".$", replacement = "", x = Right))*5)/5,
  AffinityBins = paste0("(", Left, ", ", Right, "]")
)

newplot4_a <- plotMeanAndInner(
  newplot4_dataA, CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = newplot4_dataA %>% tidytable::filter(
    PoolPatchSeed == targetSeed,
    abs(Time - targetTimes) == min(abs(Time - targetTimes)),
    Intervention == "(0)"
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 3),
  fill = ggplot2::guide_legend(ncol = 3)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), expand = FALSE
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot4_filtration %>% dplyr::pull(ID))[3]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 6, ymax = 13
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(newplot4_dataA$Start),
  xmax = max(newplot4_dataA$Stop),
  ymin = 0, ymax = max(newplot4_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.25, 0.09),
  plot.tag.position = c(0.02, 1)
)

# Trying to support a breakdown of richness and abundance to show what is
# happening to, and because of, the species that have the wrong land-use
# preference for the given habitat.

newplot4_b <- ggplot2::ggplot(
  newplot4_dataB ,
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = Intervention,
    shape = AffinityBin
  )
) + ggplot2::geom_point(
  show.legend = TRUE
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::coord_fixed(
) + ggplot2::labs(
  x = "Basal Abundance",
  y = "Consumer Abundance"
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme(
  plot.tag.position = c(0.02, 1)
) + ggplot2::guides(
  color = "none", fill = "none"
) + ggplot2::labs(
  shape = "Land-use\nPreference"
)

newplot4_c <- ggplot2::ggplot(
  newplot4_dataC,
  ggplot2::aes(
    x = Persistence,
    y = Size,
    shape = AffinityBins,
    color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_hline(
  yintercept = 10^-1,
  linetype = "dashed"
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "c)"
) + ggplot2::theme(
  plot.tag.position = c(0.02, 1)
)

newplot4 <- ggpubr::ggarrange(
  plotlist = list(
    newplot4_a,
    ggpubr::ggarrange(newplot4_b, newplot4_c, ncol = 1)
  ), nrow = 1, widths = c(0.55, 0.45)
)
ggplot2::ggsave(plot = newplot4, filename = "aplot4.png",
                units = "cm", width = 6.5*3, height = 6.5*2)

# New Plot 5: #################################################################
newplot5_filtration <- function(.) {tidytable::filter(
  .,
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)}

newplot5_dataA <- diversitiesAll %>% newplot5_filtration(
) %>% tidytable::filter(
  Metric == "Alpha Hill:0",
  is.na(Subset)
) %>% tidytable::left_join(endTimes %>% dplyr::select(-Times))

newplot5_dataB <- diversitiesAll %>% newplot5_filtration(
) %>% tidytable::filter(
  Metric == "TimeJaccard",
  is.na(Subset)
) %>% tidytable::left_join(endTimes %>% dplyr::select(-Times))

newplot5_dataC <- tidytable::bind_rows(
  Pers %>% newplot2_filtration(),
  Pers %>% newplot5_filtration() %>% tidytable::filter(Intervention != "(0)")
) %>% tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) %>% tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
)

newplot5_a <- plotMeanAndInner(
  newplot5_dataA, CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_point(
  data = newplot5_dataA %>% tidytable::filter(
    PoolPatchSeed == targetSeed
  ) %>% tidytable::group_by(
    Intervention
  ) %>% tidytable::mutate(
    TimeDist = abs(Time + 1e-6 - targetTimes), # Preference a side.
    IsMin = TimeDist == min(TimeDist)
  ) %>% tidytable::filter(
    IsMin
  )
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::guides(
  linetype = "none",
  color = ggplot2::guide_legend(ncol = 3),
  fill = ggplot2::guide_legend(ncol = 3)
) + ggplot2::coord_cartesian(
  xlim = c(0, 40000), expand = FALSE
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot5_filtration %>% dplyr::pull(ID))[1]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 7, ymax = 17
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot5_filtration %>% dplyr::pull(ID))[2]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 18, ymax = 28
) + ggplot2::annotation_custom(
  ggplot2::ggplotGrob(
    targetEnvs[[
      (targetEnvsIndex %>% newplot5_filtration %>% dplyr::pull(ID))[3]
      ]]$singletonGraphs[[1]][[1]] + ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  ),
  xmin = 30500, xmax = 40000, ymin = 29, ymax = 39
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(newplot5_dataA$Start),
  xmax = max(newplot5_dataA$Stop),
  ymin = 0, ymax = max(newplot5_dataA$Value),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  tag = "a)"
) + ggplot2::theme(
  legend.position = c(0.25, 0.09),
  plot.tag.position = c(0.05, 1)
)

newplot5_ba <- plotMeanAndInner(
  newplot5_dataB, CIs = 0.75, facets = as.formula(. ~ .)
) + ggplot2::geom_rect(
  data = data.frame(
    1 # 1 rectangle per row, so dummy df to prevent overplotting
  ),
  xmin = min(newplot5_dataB$Start),
  xmax = max(newplot5_dataB$Stop),
  ymin = 0, ymax = max(newplot5_dataB$Value, na.rm = TRUE),
  fill = "grey",
  alpha = 0.2,
  inherit.aes = FALSE
) + ggplot2::labs(
  y = "Time-Jaccard",
  tag = "b)"
) + ggplot2::theme(
  legend.position = "none",
  plot.tag.position = c(0.05, 1)
) + ggplot2::coord_cartesian(
  xlim = c(0, 30000), ylim = c(0, 0.5), expand = FALSE
)
newplot5_bb <- ggplot2::ggplot(
  newplot5_dataB %>% tidytable::filter(
    Time > Start, Time < Stop
  ),
  ggplot2::aes(
    y = Value, color = Intervention
  )
) + ggplot2::geom_density(
  adjust = 1/2, show.legend = FALSE
) + ggplot2::coord_cartesian(
  ylim = c(0, 0.5), expand = FALSE
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::scale_x_continuous(
  name = "", labels = function(x) rep("", length(x))
) + ggplot2::theme(
  axis.title.y = ggplot2::element_blank(),
  axis.text.y = ggplot2::element_blank(),
  axis.ticks.y = ggplot2::element_blank()#,
  # axis.text.x = ggplot2::element_blank()
)
newplot5_b <- ggpubr::ggarrange(
  newplot5_ba, newplot5_bb,
  nrow = 1, widths = c(0.9, 0.1), align = "h"
)

ggplot2::ggplot(
  newplot5_dataC,
  ggplot2::aes(
    x = Intervention,
    fill = interaction(InType, OutType)
  )) + ggplot2::geom_bar(
  )

# Statistics: #################################################################
# SI FIGURE?
# Pers %>% newplot2_filtration(
# ) %>% tidytable::filter(
#   In < Stop, Out > Start # Not things outside of [Start, Stop]
# ) %>% tidytable::mutate(
#   # Shorten intervals for equivalent comparisons.
#   InType = ifelse(In < Start, "Persistent", InType),
#   OutType = ifelse(Out > Stop, "Persistent", OutType),
#   In = ifelse(In < Start, Start, In),
#   Out = ifelse(Out > Stop, Stop, Out),
#   Persistence = Out - In
# ) %>% ggplot2::ggplot(
#   ggplot2::aes(x = PoolPatchSeed, fill = OutType)
# ) + ggplot2::geom_bar(
# ) + ggplot2::facet_grid(
#   . ~ Intervention
# )
# Change in outtypes for Figure 2.
Pers %>% newplot2_filtration(
) %>% tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) %>% tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) %>% with(table(OutType, Intervention))

# Change in amount of abundance between basals.
# Note we are pairing time-averages of simulations, then dividing, then
# averaging, but this ignores the internal (averaged-out) time dynamics, which
# judging by the inner 50% of values over time, appears to be quite variable.
newplot2_dataB %>% tidytable::select(
  -Affinity, -AffinitySeed, -InterventionInitial, -InterventionFinal
) %>% tidytable::ungroup()  %>% tidytable::pivot_wider(
  values_from = c(q025_Basal, Mean_Basal, q075_Basal,
                  q025_Consumer, Mean_Consumer, q075_Consumer),
  names_from = Intervention
) %>% tidytable::distinct(
) %>% tidytable::mutate(
  RatioB_1_05 = `Mean_Basal_(1)` / `Mean_Basal_(0.5)`,
  RatioB_0_1 = `Mean_Basal_(0)` / `Mean_Basal_(1)`,
  RatioC_05_1 = `Mean_Consumer_(0.5)` / `Mean_Consumer_(1)`,
  RatioC_0_1 = `Mean_Consumer_(0)` / `Mean_Consumer_(1)`
) %>% tidytable::summarise(
  MinB_1_05 = min(RatioB_1_05),
  RatioB_1_05 = mean(RatioB_1_05),
  MaxB_1_05 = max(RatioB_1_05),
  MinB_0_1 = min(RatioB_0_1),
  RatioB_0_1 = mean(RatioB_0_1),
  MaxB_0_1 = max(RatioB_0_1),
  MinC_05_1 = min(RatioC_05_1),
  RatioC_05_1 = mean(RatioC_05_1),
  MaxC_05_1 = max(RatioC_05_1),
  MinC_0_1 = min(RatioC_0_1),
  RatioC_0_1 = mean(RatioC_0_1),
  MaxC_0_1 = max(RatioC_0_1)
)

# Show that the Fig3,4 persistences are fairly transient in comparison.
ggplot2::ggplot(
  Pers %>% newplot2_filtration(
  ) %>% tidytable::filter(
    In < Stop, Out > Start # Not things outside of [Start, Stop]
  ) %>% tidytable::mutate(
    # Shorten intervals for equivalent comparisons.
    InType = ifelse(In < Start, "Persistent", InType),
    OutType = ifelse(Out > Stop, "Persistent", OutType),
    In = ifelse(In < Start, Start, In),
    Out = ifelse(Out > Stop, Stop, Out),
    Persistence = Out - In
  ) %>% tidytable::group_by(
    Species, Environment, SpeciesType, Size, ReproductionRate, Speed, AffinityBins,
    PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
    InitialConditions, InitialConditionsSeed, DispersalParam, NicheDistance,
    Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
    InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
    InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
  ) %>% tidytable::summarise( # Sum over Appearances.
    Persistence = sum(Persistence),
    .groups = "drop"
  ) %>% tidytable::group_by(
    Environment, SpeciesType, AffinityBins,
    PoolPatch, Interactions, Events,
    InitialConditions, DispersalParam, NicheDistance,
    Affinity, InterventionPatchType,
    InterventionTimeType, InterventionDispersal,
    InterventionNicheDistance, Intervention, SpeciesAffinity, Start, Stop
  ) %>% tidytable::summarise( # Sum over Appearances.
    Persistence = mean(Persistence),
    Size = 10^(mean(log10(Size))),
    .groups = "drop"
  ) -> temp,
  ggplot2::aes(
    x = Persistence,
    y = Size,
    shape = AffinityBins,
    color = Intervention
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::geom_hline(
  yintercept = 10^-1,
  linetype = "dashed"
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
)

temp %>% group_by(
  SpeciesType, Intervention, AffinityBins
) %>% summarise(
  min(Persistence), max(Persistence)
)

newplot3_dataC %>% group_by(
  SpeciesType, Intervention, AffinityBins
) %>% summarise(
  min(Persistence), max(Persistence)
)

# Decide how to statistically test the data.
# diversitiesAll %>% tidytable::filter(
#   NicheDistance == "5",
#   Intervention %in% c("(0)", "(0.5)", "(1)"),
#   (PoolPatchSeed %in% as.character(343:386)),
#   Metric == "Alpha Hill:0"
# ) %>% ari

diff_1000_5050 <- diversitiesAll %>% tidytable::filter(
  NicheDistance == "5",
  Intervention %in% c("(0)"),
  SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1"),
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) %>% tidytable::group_by(Time, PoolPatchSeed) %>% tidytable::mutate(
  Value = ifelse(SpeciesAffinity == "100% 0", # Reference
                 -Value, Value)
) %>% tidytable::summarise(
  Value = sum(Value)
)
ggplot2::ggplot(
  diff_1000_5050,
  ggplot2::aes(
    x = Time, y = Value#, group = PoolPatchSeed
  )
) + ggplot2::geom_hex()
nlme::gls(
  Value ~ 1,
  data = diff_1000_5050,
  correlation = nlme::corAR1(0,form = ~ Time|PoolPatchSeed)
) %>% summary

diff_nopref_unif <- tidytable::left_join(
  diversitiesAll %>% tidytable::filter(
    NicheDistance == "5",
    Intervention %in% c("(0)", "(0.5)", "(1)"),
    SpeciesAffinity %in% c("Uniform(0, 1)"),
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) %>% tidytable::rename(
    Value1 = Value,
    Affinity1 = SpeciesAffinity,
    Intervention1 = Intervention
  ),
  diversitiesAll %>% tidytable::filter(
    NicheDistance == "5",
    Intervention %in% c("(0.5)"),
    SpeciesAffinity %in% c("100% 0"),
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) %>% tidytable::rename(
    Value2 = Value,
    Affinity2 = SpeciesAffinity
  ) %>% tidytable::select(
    -tidytable::starts_with("Intervention"),
    -Affinity, -AffinitySeed
  )
) %>% tidytable::mutate(
  Value = Value2 - Value1,
  Group = paste(Intervention1, PoolPatchSeed)
)
diff_nopref_unif %>% dplyr::group_by(Intervention1) %>% dplyr::group_map(
  function(x, y) nlme::gls(
    Value ~ 1,
    data = x,
    correlation = nlme::corAR1(0, form = ~ Time|Group)
  ) %>% summary()
)

diff_unif_5v0 <- tidytable::left_join(
  diversitiesAll %>% tidytable::filter(
    NicheDistance == "5",
    Intervention %in% c("(0)", "(1)"),
    SpeciesAffinity %in% c("Uniform(0, 1)"),
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) %>% tidytable::rename(
    Value1 = Value,
    Affinity1 = SpeciesAffinity,
    Intervention1 = Intervention
  ),
  diversitiesAll %>% tidytable::filter(
    NicheDistance == "5",
    Intervention %in% c("(0.5)"),
    SpeciesAffinity %in% c("Uniform(0, 1)"),
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) %>% tidytable::rename(
    Value2 = Value,
    Affinity2 = SpeciesAffinity
  ) %>% tidytable::select(
    -tidytable::starts_with("Intervention"),
    -Affinity, -AffinitySeed
  )
) %>% tidytable::mutate(
  Value = Value2 - Value1,
  Group = paste(Intervention1, PoolPatchSeed)
)
diff_unif_5v0 %>% dplyr::group_by(Intervention1) %>% dplyr::group_map(
  function(x, y)
    list(y, nlme::gls(
      Value ~ 1,
      data = x,
      correlation = nlme::corAR1(0, form = ~ Time|Group)
    ) %>% summary())
)

diff_unif_to0 <- tidytable::left_join(
  diversitiesAll %>% tidytable::filter(
    NicheDistance == "5",
    Intervention %in% c("(0.5)->(0)", "(1)->(0)"),
    SpeciesAffinity %in% c("Uniform(0, 1)"),
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) %>% tidytable::rename(
    Value1 = Value,
    Affinity1 = SpeciesAffinity,
    Intervention1 = Intervention
  ) %>% tidytable::mutate(
    Time1 = round(Time, -1)
  ),
  diversitiesAll %>% tidytable::filter(
    NicheDistance == "5",
    Intervention %in% c("(0)"),
    SpeciesAffinity %in% c("Uniform(0, 1)"),
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset)
  ) %>% tidytable::mutate(
    Time = round(Time, -1)
  ) %>% tidytable::rename(
    Time1 = Time, # So that they match up when there are multiple times in left.
    Value2 = Value,
    Affinity2 = SpeciesAffinity
  ) %>% tidytable::select(
    -tidytable::starts_with("Intervention"),
    -Affinity, -AffinitySeed
  )
) %>% tidytable::mutate(
  Value =  Value1 - Value2,
  Group = paste(Intervention1, PoolPatchSeed)
) %>% tidytable::filter(
  Time >= 20000, Time <= 30000
)
diff_unif_to0  %>% dplyr::group_by(Intervention1) %>% dplyr::group_map(
  function(x, y)
    list(y, nlme::gls(
      Value ~ 1,
      data = x,
      correlation = nlme::corAR1(0, form = ~ Time|Group)
    ) %>% summary())
)

# MIGHT REQUIRE VIKING: ########################################################
# Can we look for clusters? We remove transient species since they don't have
# time to change the ecosystem generally.
# %>% tidytable::pull(Persistence) %>% sort() %>% plot()
temp <- Pers %>% tidytable::filter(
  SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)",
                      "(0)->(0.5)", "(0)->(1)",
                      "(0.5)->(0)", "(0.5)->(1)",
                      "(1)->(0)", "(1)->(0.5)"),
  (PoolPatchSeed %in% as.character(343:386)),
  In < Stop, Out > Start, Persistence >= 100
) %>% tidytable::separate(
  col = "AffinityBins", into = c("Left", "Right"), sep = ","
) %>% tidytable::mutate(
  Left = round(as.numeric(gsub(pattern = "^.", replacement = "", x = Left))*5)/5,
  Right = round(as.numeric(gsub(pattern = ".$", replacement = "", x = Right))*5)/5,
  AffinityBins = paste0("(", Left, ", ", Right, "]")
)
# Then we use something like the ECDF over sizes to characterise the ecosystems.
# (Size because this characterises the food chain. Abundance is also possibly
#  reasonable because it relates to size, but its also affected by preference.)
# The trick is to do this through time and Pers doesn't have direct access.
temp2 <- temp %>% dplyr::group_by(
  Environment,
  dplyr::across(PoolPatchSeed:InterventionNicheDistance),
  SpeciesAffinity, Intervention
) %>% dplyr::group_map(
  .f = function(values, key) {
    times <- range(values$Start, values$Stop)
    # print(times)
    times <- seq(times[1], times[2], by = 1000) # Reduce by for accuracy.
    ecdfs <- vector("list", length = length(times))
    for (i in seq_along(times)) {
      time <- times[i]
      vals <- values %>% tidytable::filter(
        In < time, time < Out
      ) %>% tidytable::pull(
        Size
      )
      if (length(vals) > 0) {
        # Why consumer ecdf(0) if there isn't anything?
        # log(0) = -Inf, so it makes some sense on the natural scale,
        # and it also means that we penalise more heavily the difference
        # between no consumers and large consumers than no and small.
        ecdfs[[i]] <- list(
          Basals = if (sum(vals < 10^-1) > 0) {
            ecdf(vals[vals < 10^-1])
          } else {
            function(x) {rep(NA, length(x))}
          },
          Consumers = if (sum(vals >= 10^-1) > 0) {
            ecdf(vals[vals >= 10^-1])
          } else {
            ecdf(0)
          }
        )
      } else {
        ecdfs[[i]] <- list(
          Basals = function(x) {rep(NA, length(x))},
          Consumers = ecdf(0)
        )
      }
    }
    return(list(
      Key = key, Times = times, ECDFs = ecdfs
    ))
  }
)

# Spot Check
with(list(pick = 10),
     expand.grid(
       Index = 1:length(temp2[[pick]]$Times),
       Sizes = 10^(seq(-2.1, 0.6, by = 0.1))
     ) %>% mutate(
       Times = temp2[[1]]$Times[Index],
       ECDF = unlist(mapply(
         function(i, s) (
           temp2[[pick]]$ECDFs[[i]]$Basal(s) +
             temp2[[pick]]$ECDFs[[i]]$Consumers(s)
           ),
         i = Index, s = Sizes
       ))
     ) %>% ggplot2::ggplot(
       ggplot2::aes(
         x = Times,
         y = Sizes,
         color = ECDF,
         group = Times)
     ) + ggplot2::geom_line(
     ) + ggplot2::scale_color_steps2(
       midpoint = 1
     ) + ggplot2::scale_y_log10(
     ) + ggplot2::coord_cartesian(
       ylim = 10^c(-2, 0.5),
       xlim = c(20000, 30000)
     )
)


# We need to assign indices to all of the ECDFS, then take their distances
# and then attempt to cluster them to see if we get any signal.
# We'll start with ks-metric.
index_ecdf <- 0
index_coords <- data.frame()
for (index_group in seq_along(temp2)) {
  temp2[[index_group]]$Index <-
    1:length(temp2[[index_group]]$Times) + index_ecdf
  index_ecdf <- temp2[[index_group]]$Index[length(temp2[[index_group]]$Index)]
  index_coords <- rbind(
    index_coords,
    data.frame(
      Index = temp2[[index_group]]$Index,
      First = index_group,
      Second = 1:length(temp2[[index_group]]$Times)
    ))
}

# Can't directly call ks.test without recomputing things a lot.
matrix_ecdf_dist <- outer(
  X = 1:index_ecdf, Y = 1:index_ecdf,
  FUN = Vectorize(function(i, j) {
    if (i <= j) {return(0)} # Save work, and diag should be 0.
    it <- temp2[[index_coords[i,]$First]]$ECDFs[[index_coords[i,]$Second]]
    jt <- temp2[[index_coords[j,]$First]]$ECDFs[[index_coords[j,]$Second]]

    # Logs help respect the different scales that basals and consumers are on.
    # Note: this might be sufficient with emdw, but splitting Bs and Cs means
    # we have to pay equal attention to both sides.
    emdB <- if (length(environment(it$Basal)$x) == 0 ||
                length(environment(jt$Basal)$x) == 0) {
      NaN
    } else {
      emdist::emdw(
        A = log10(environment(it$Basal)$x),
        B = log10(environment(jt$Basal)$x),
        wA = diff(c(0, it$Basal(environment(it$Basal)$x))),
        wB = diff(c(0, jt$Basal(environment(jt$Basal)$x)))
      )
    }
    emdC <-  if (length(environment(it$Consumer)$x) == 0 ||
                 length(environment(jt$Consumer)$x) == 0) {
      NaN
    } else {
      emdist::emdw(
        A = log10(environment(it$Consumer)$x),
        B = log10(environment(jt$Consumer)$x),
        wA = diff(c(0, it$Consumer(environment(it$Consumer)$x))),
        wB = diff(c(0, jt$Consumer(environment(jt$Consumer)$x)))
      )
    }
    # KS - Maximal CDF Difference
    # testvals <- unique(
    #   environment(it)$x,
    #   environment(jt)$x
    # )
    # return(max(abs(it(testvals) - jt(testvals))))
    # Earth Mover - Rearrange PDF
    return(emdB + emdC)
  }))
matrix_ecdf_dist <- matrix_ecdf_dist + t(matrix_ecdf_dist)
# Not sure about best storage method.
# matrix_ecdf_dist <- Matrix::sparseMatrix(
#   dims = c(index_ecdf, index_ecdf),
#   symmetric = TRUE
#   )

color_codes <- do.call(rbind, lapply(
  seq_along(temp2),
  function(i) with(temp2[[i]]$Key,
                   data.frame(
                     First = i,
                     Key = paste(Intervention, SpeciesAffinity, sep = ", ")))
  )) %>% dplyr::mutate(
    KeyIndex = as.numeric(factor(Key)),
    Color = colorspace::qualitative_hcl(length(unique(Key)))[KeyIndex]
  )
index_coords <- index_coords %>% dplyr::left_join(color_codes)

matrix_ecdf_dist_nona <- matrix_ecdf_dist
matrix_ecdf_dist_nona[is.na(matrix_ecdf_dist_nona)] <- Inf
  # 1e3 * max(matrix_ecdf_dist, na.rm = T)
matrix_ecdf_dist_omitna <- matrix_ecdf_dist[
  rowSums(is.na(matrix_ecdf_dist)) != ncol(matrix_ecdf_dist) - 1,
  colSums(is.na(matrix_ecdf_dist)) != nrow(matrix_ecdf_dist) - 1
  ] # stackoverflow.com/a/6437778

dianaout <- cluster::diana(matrix_ecdf_dist_omitna, diss = TRUE, keep.diss = TRUE)
agnesout <- cluster::agnes(matrix_ecdf_dist_omitna, diss = TRUE, keep.diss = TRUE)
lapply(1:20, function(i) cluster::pam(matrix_ecdf_dist_nona, diss = TRUE, k = i)$silinfo$avg.width) %>% do.call(what = rbind) %>% plot
lapply(2:20, function(i) data.frame(index = i, value = mean((cluster::silhouette(cutree(dianaout, k = i), dianaout$diss))[, 3]))) %>% do.call(what = rbind) %>% points(col = "blue")
lapply(2:20, function(i) data.frame(index = i, value = mean((cluster::silhouette(cutree(agnesout, k = i), agnesout$diss))[, 3]))) %>% do.call(what = rbind) %>% points(col = "red")

plot(cluster::silhouette(cutree(dianaout, k = 4), dianaout$diss))
plot(cluster::silhouette(cutree(agnesout, k = 4), agnesout$diss))

ggplot2::ggplot(
  ggplot2::aes(y = Group, x = Key, color = Key),
  data = cbind(
    index_coords[-which(rowSums(is.na(matrix_ecdf_dist)) == ncol(matrix_ecdf_dist) - 1),],
    Group = factor(cutree(agnesout, k = 7))
  )
) + ggplot2::geom_point(position = "jitter")
# Suggests that I probably need to be swapping my statistics.
# I expected something with all basals to be quite distinct from a mixed
# basal-consumer system. Meanwhile the Uniform(0, 1) seem to enter multiple
# states from this depiction (which might be true if there is a build-up or
# some form of cycle) but it seems odd it would have an something similar to an
# all basal state!
# The problem is how to capture this difference.
# The obvious answer is to explicitly consider the two separately
# (which necessitates thinking about how to combine the two together again).
# Maybe functional diversity literature has an answer (for when a function is
# entirely missing)?
# If we're comparing ecdfs, setting all of the mass of the missing bit at -Inf
# or 0 makes sense; this makes a system with one small consumer more like the
# no consumer system than one with a variety of consumers or large consumers.
# These aren't realisations of a joint distribution. They are possibly
# realisation from a Dirichlet Process, for which the parameters change between
# simulations (and the question is how many DPs would we need, and do we need
# more than one for a single set of parameters).

# Probably should switch to silhouette widths, and compare with cluster::diana.

# plot(cluster::pam(matrix_ecdf_dist_nona, diss = TRUE, k = 7, keep.diss = TRUE, pamonce = 6), which.plots = 1,
#      col.p = index_coords$Color)
# # KS - Probably not good enough to detect differences.
# # Clustering algorithm or the distance metric?
# # Hard to say, we can see the 3 main clusters emerging as we might expect
# # (red is 0.5, green is 0, blue is 1), 0.5 and 0 are closer to each other than
# # to 1. But we're not seeing separation within these main clusters when we
# # suspect there might be.

with(
  list(dat = diversitiesAll %>% tidytable::filter(
    NicheDistance == "5",
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset),
    Time != round(Time) | lag(Time) != round(lag(Time))
  )),
  dat %>% group_by(
    PoolPatchSeed, Intervention, SpeciesAffinity
  ) %>% summarise(
    Value = diff(Value),
    Time = diff(Time)
  ) %>% ggplot2::ggplot(
    ggplot2::aes(
      x = Time,
      y = Value,
      color = Intervention,
      shape = SpeciesAffinity,
      linetype = SpeciesAffinity
    )
  ) + ggplot2::geom_point(
  ) + ggplot2::geom_smooth(
  ) + ggplot2::facet_wrap(
    . ~ SpeciesAffinity
  )
)

with(
  diversitiesRichness %>% tidytable::filter(
    NicheDistance == "5",
    (PoolPatchSeed %in% as.character(343:386)),
    Metric == "Alpha Hill:0",
    is.na(Subset),
    Time != round(Time) | lag(Time) != round(lag(Time))
  ) %>% changeAffinityLevels() %>% changeInterventionLevels(
  ) %>% tidytable::mutate(
    #TODO consider moving this into changeInterventionLevels.
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
    ),
    InterventionDist = abs(
      as.numeric(gsub(InterventionFinal, pattern = "[(|)]", replacement = "")) -
        as.numeric(gsub(InterventionInitial, pattern = "[(|)]", replacement = ""))
      )
  ) %>% group_by(
      PoolPatchSeed, Intervention, InterventionDist, SpeciesAffinity
    ) %>% summarise(
      Value = round(diff(Value)),
      Time = diff(Time)
    ),
    table(InterventionDist, SpeciesAffinity, sign(Value))
  )
