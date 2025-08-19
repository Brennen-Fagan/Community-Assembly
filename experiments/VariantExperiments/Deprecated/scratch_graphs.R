# Trying to build my panel 2s of the networks, but such that if I need to,
# I only have to change the folder that I am pointing at.
# source("flattenDiversity.R")

targetDir <- dir(pattern = "TSTS_Simulations_142486-4929_343-343_2025-01-21",
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
    SpeciesAffinity =
      affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
  )
  e[[i]]
},
e = targetEnvs, d = targetDivs, s = targetFiles)
targetEnvsPool <- new.env()
load(targetPool, envir = targetEnvsPool)

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

plotGraph <- function(graph, mainLayout) {
  ggraph::ggraph(
    graph = graph,
    layout = graph %N>% data.frame(
    ) %>% select(node) %>% left_join(
      mainLayout %>% select(x, y, node)
    )
  ) + ggraph::geom_node_point(
    mapping = aes(
      color = Type,
      size = log10(N)
    )
  ) + ggraph::geom_edge_diagonal(
    mapping = aes(
      color = Type,
      linetype = Type,
      alpha = log10(effectNormalised)
    ),
    # linewidth = 2,
    arrow = arrow(length = unit(2, 'mm')),
    start_cap = circle(5, 'mm'),
    end_cap = circle(5, 'mm')
  ) + ggplot2::geom_hline(
    yintercept = -1, linetype = "dashed", color = "black"
  ) + theme_minimal(
  ) + ylab(
    "Log10(Size)"
  ) + xlab(
    "Land-use Preference"
  ) + scale_color_manual(
    values = c("limegreen", "goldenrod2")
  ) + lims(
    x = c(0, 1), y = c(-2, 0.25)
  )
}

targetEnvs <- lapply(targetEnvs, function(env) {
  timelist <- env$graphs
  env$singletonGraphs <- lapply(timelist, function(patchlist) {
    lapply(patchlist$graphs, plotGraph, mainLayout = env$layout)
  })
  return(env)
})
ggpubr::ggarrange(plotlist = unlist(unlist(singletonGraphs, recursive = FALSE), recursive = FALSE), common.legend = TRUE)
# To plot, I'm thinking Richness over Time, two lines, then three vertical lines
# indicating where we will be showing the different system configuration. As a
# companion, the two sets of systems as food webs at each point in time, with
# a red line to denote the difference in basal/consumer, size as the y-axis,
# and an x-axis that is a jittered affinity. Color the nodes by their actual
# affinity. Select the three time points so right after stabilisation,
# a very short time after intervention/land-use-change, and far afterwards.
# I should also add a marker of the land-use on the bottom of the graph plots.
# We'll go with the richness over time as short and wide on the bottom.

# diversitiesExample <- tidytable::bind_rows(
#   lapply(loadenvs[1:4], function(env) env$Diversity)
# ) %>% changeAffinityLevels(
# ) %>% changeInterventionLevels(
# )
#
# # Conversion appears to be straightforward.
#

#
# # Not as diverse as I might like obviously, so maybe there are better choices.
# singletonRichness <- ggplot2::ggplot(
#   diversitiesExample %>% tidytable::filter(
#     Metric == "Alpha Hill:0", is.na(Subset)
#   ), ggplot2::aes(
#     x = Time, y = Value, color = Intervention,
#     group = interaction(Intervention, SpeciesAffinity)
#   )
# ) + ggplot2::geom_line(
# ) + ggplot2::geom_vline(
#   xintercept = targetTimes
# ) + ggplot2::theme_minimal(
# )
#

#


# Would need to add in blanks for the missing (because nonintervention) plots.
#
