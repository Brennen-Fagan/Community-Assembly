# Plotting with both ColExt and Divs
load("ColExt9a9_flat.RData")
load("diversitiesFlattened9a9_subset3.RData")

# Problems with X11
options(bitmapType = "cairo")

# Grey interval that we compute over, usually after intervention (~50%)
# If second number is less than 1, we lose persistent species.
end <- c(0.6, 0.904) # Aiming for 20000 - 30000

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
  "Uniform(0, 1)" = "dotted"
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
        "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
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
    color = "none",
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
  ) + lims(
    x = c(0, 1), y = c(-2, 0.25)
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
diversitiesAll <- diversitiesAll %>% changeAffinityLevels(
) %>% changeInterventionLevels(
) %>% tidytable::mutate(
  #TODO consider moving this into changeInterventionLevels.
  InterventionInitial = tidytable::case_when(
    Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)") ~ "(0)",
    Intervention %in% c("(0.5)->(0)", "(0.5)", "(0.5)->(1)") ~ "(0.5)",
    Intervention %in% c("(1)->(0)", "(1)->(0.5)", "(1)") ~ "(1)",
    TRUE ~ NA_character_
  ),
  InterventionInitial = factor(
    InterventionInitial,
    levels = c(
      "(0)", "(0.5)",  "(1)"
    ), ordered = TRUE
  ),
  InterventionFinal = tidytable::case_when(
    Intervention %in% c("(0)", "(0.5)->(0)", "(1)->(0)") ~ "(0)",
    Intervention %in% c("(0)->(0.5)", "(0.5)",  "(1)->(0.5)") ~ "(0.5)",
    Intervention %in% c("(0)->(1)", "(0.5)->(1)", "(1)") ~ "(1)",
    TRUE ~ NA_character_
  ),
  InterventionFinal = factor(
    InterventionFinal,
    levels = c(
      "(0)", "(0.5)",  "(1)"
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



# Panels:
newplot2_filtration <- function(.) {tidytable::filter(
  .,
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  (PoolPatchSeed %in% as.character(343:386))
)}

newplot2_dataA <- diversitiesAll %>% newplot2_filtration() %>% tidytable::filter(
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
) + ggplot2::theme(
  legend.position = c(0.2, 0.09)
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
        axis.ticks.x = ggplot2::element_blank()
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
        axis.ticks.x = ggplot2::element_blank()
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
        axis.ticks.x = ggplot2::element_blank()
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
)

newplot2_b <- ggplot2::ggplot(
  newplot2_dataB,
  ggplot2::aes(
    x = Mean_Basal,
    y = Mean_Consumer,
    color = Intervention)
  ) + ggplot2::geom_point(
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
    x = "Avg. Total Basal Abundance",
    y = "Avg. Total Consumer Abundance"
  ) + ggplot2::theme_minimal(
  ) + ggplot2::scale_color_manual(
    values = colorPalette, aesthetics = c("color", "fill"),
    name = "Habitat Land-use"
  )

# Size as a function of Total amount of time spent in simulation
# Color is intervention type.

# Average Temporal Beta as a function of Intervention Type.
