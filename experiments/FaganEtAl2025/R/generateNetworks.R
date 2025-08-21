generateNetworks <- function(
  Specification,
  dIS = diversitiesInterventionStrings,
  aDO = affinityDictionaryOrigin
) {
  # Goal: Pass data.frame, return list of networks where a network corresponds
  #       to a row of the data.frame.
  # NOTE: This will need to be updated for Version 10. Currently Version 9.
  # Tested singleton and multiple versions; appears to be working fine.
  columns <- c(
    # Which network:
    "Time", "Environment1",
    # Which File (Base):
    "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
    "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
    "Dispersal", "NicheDistance", "Affinity", "AffinitySeed",
    # Which File (Intervention):
    "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
    "InterventionTimeSeed", "InterventionDispersal", "InterventionNicheDistance"
  )
  stopifnot(
    columns %in% names(Specification),
    length(unique( # Only 1 Unique Length amongst the columns
      unlist(lapply(columns, function(clm) length(Specification[[clm]])))
    )) == 1
  )

  ## Load Targets: #############################################################
  targetDirs <- sapply(
    with(Specification, paste0(
      "TSTS_Simulations_",
      PoolPatch, "-", Interactions, "_",
      PoolPatchSeed, "-", InteractionsSeed, "_"
    )),
    dir,
    path = ".",
    full.names = TRUE
  )

  targetFiles <- dir(targetDirs, pattern = "(Sim|Int)", full.names = T)

  ### Retrieval Base: ##########################################################
  targetFiles <- lapply(
    with(Specification, paste0(
      PoolPatch, "-", Interactions, "-", Events, "-",
      InitialConditions, "-", Dispersal, "-", NicheDistance, "-", Affinity, "_",
      PoolPatchSeed, "-", InteractionsSeed, "-", EventsSeed, "-",
      InitialConditionsSeed, "-", AffinitySeed# , "_"
    )),
    function(pat, ...) unique(grep(pat, ...)),
    x = targetFiles,
    value = TRUE
  )

  targetFiles <- lapply(
    seq_along(targetFiles), function(i, f, s) {
      grep(
        if (!is.na(s[i, ]$InterventionPatchType)) {
          with(s[i, ], paste0(
            InterventionPatchType, "-", InterventionTimeType, "-",
            InterventionDispersal, "-", InterventionNicheDistance, "_",
            InterventionPatchSeed, "-", InterventionTimeSeed
          ))
        } else {
          "_Simulation_"
        },
        x = f[[i]],
        fixed = TRUE,
        value = TRUE
      )
    },
    f = targetFiles, s = Specification
  )

  stopifnot(unlist(lapply(targetFiles, length)) == 1)
  targetFiles <- unlist(targetFiles)

  # Families of targets:
  # targetFiles <- c(targetFilesSim, targetFilesInt)
  # Need to fold together back in order of Specification's rows.
  targetDivs <- gsub(x = targetFiles,
                     pattern = "_(Simulation|Intervention)_",
                     replacement = "_Diversity_")
  targetPools <- sapply(dirname(targetFiles),
                        dir,
                        pattern = "Pool",
                        full.names = T)

  ### Perform Load: ############################################################
  # Load Diversities and prepare associations.
  targetDivsUN <- unique(targetDivs)
  targetDivsU <- lapply(targetDivsUN, function(d) {
    div <- load(d)
    div <- get(div)
    flattenDiversity(div) |> dplyr::left_join(
      dIS,
      by = c("Affinity", "PoolPatch", "InterventionPatchType"),
      multiple = "all"
    ) |> dplyr::mutate(# First Number in the Second Group (split = _)
      PoolPatchSeed =
        gsub(div$Ellipsis$ID,
             pattern = ".*?((?<=_)[0-9]+).*", # Capture a group, and everything.
             replacement = "\\1", # Replace it with just the captured group.
             perl = TRUE),
      SpeciesAffinity = aDO$SpeciesAffinities[as.numeric(Affinity)]
    ) |> changeAffinityLevels()
  })
  names(targetDivsU) <- targetDivsUN

  # Repeat for Pools
  targetPoolsUN <- unique(targetPools)
  targetPoolsU <- lapply(targetPoolsUN, function(d) {
    pool <- new.env()
    load(d, envir = pool)
    return(pool)
  })
  names(targetPoolsU) <- targetPoolsUN
  targetPoolGroups <- apply(
    outer(targetPoolsUN, targetPools, FUN = function(x, y) x==y),
    2, which
  )

  # Implement Associations
  # Can make more memory efficient by reducing at this step just to the needed
  # temporal slices.
  # Can also reduce loads if this becomes a bottleneck (why?) by treating the
  # load(s[[i]]) as with the diversities/pools above.
  targetEnvs <- replicate(length(targetFiles), new.env())
  targetEnvs <- lapply(
    seq_along(targetEnvs), function(i, e, s, d, p, r) {
      # load(d[[i]], envir = e[[i]])
      load(s[[i]], envir = e[[i]]) # loads results.
      e[[i]]$Diversity <- targetDivsU[[d[i]]] # string/name indexing.
      e[[i]]$Pool <- targetPoolsU[[p[i]]]
      e[[i]]$Row <- r[i, ]
      return(e[[i]])
    },
    e = targetEnvs,
    s = targetFiles,
    d = targetDivs,
    p = targetPools,
    r = Specification
  )

  # Reconstruct environment
  targetEnvs <- lapply(targetEnvs, function(env) {
    intervention <- #T/F
      !("EffectiveReproductionRate" %in% names(env$result$Ellipsis$Affinity))

    env$calculator <- with(
      env$result,
      CalculateTrophicStructure(
        Pool = env$Pool$Pool,
        ReproductionRate =
          if (!intervention) {
            Ellipsis$Affinity$EffectiveReproductionRate
          } else {
            Ellipsis$Affinity$EffectiveReproductionRateIntervention
          },
        NumEnvironments = NumEnvironments,
        InteractionMatrices = env$Pool$InteractionMatrices,
        EliminationThreshold = Parameters$EliminationThreshold,
        LinkThreshold = 0.01
      )
    )

    env$trophics <- with(
      env$result,
      if (intervention &&
          env$Row$Time < Abundance[1, 1] / ReactionTime) {
        NULL
      } else {
        timerow <-
          which.max(Abundance[, 1] / ReactionTime > env$Row$Time) # first match
        retval <- env$calculator(Abundance[timerow, 1], Abundance[timerow, -1])
        retval$Time <- Abundance[timerow, 1] / ReactionTime
        retval
      }
    )

    return(env)
  })

  ## CREATE NETWORKS: #########################################################
  # Should be easy from here out?
  # May need to handle rearrangement at the end to make sure promise above kept.

  targetEnvs <- lapply(targetEnvs, function(env) {
    env$graphs <- lapply(env$trophics$EdgeVertexLists[[1]], function(patch) {
      tidygraph::tbl_graph(
        nodes = patch$Vertices, edges = patch$Edges
      )
    })
    return(env)
  })

  # Use the target pool groups to create layouts.
  targetEnvsLayouts <- lapply(
    unique(targetPoolGroups),
    function(i) {
      indices <- which(i == targetPoolGroups)
      # Within patches.
      env_undirecteds <- lapply(targetEnvs[indices], function(env) {
        env_undirected <- lapply(env$graphs, tidygraph::to_undirected)
        env_undirected <- lapply(env$graphs, tidygraph::select, "node")
        if (length(env_undirected) > 1) {
          for (j in 2:length(env_undirected)) {
            env_undirected[[1]] <- tidygraph::graph_join(
              env_undirected[[1]], env_undirected[[j]], by = "node"
            )
          }
        }
        return(env_undirected)
      })
      # Repeat across targets.
      env_undirected <- env_undirecteds[[1]][[1]]
      if (length(env_undirecteds) > 1) {
        for (j in 2:length(env_undirecteds)) {
          env_undirected <- tidygraph::graph_join(
            env_undirected, env_undirecteds[[j]][[1]], by = "node"
          )
        }
      }

      if ((env_undirected %>% mutate(
        Components = tidygraph::graph_component_count(),
        Nodes = tidygraph::graph_order(),
        NotSingletons = Components != Nodes
      ) %>% pull(NotSingletons))[1]
      ) {
        # need to undirect again (because recording can be either direction?)
        env_undirected <- tidygraph::to_undirected(env_undirected)
        env_undirected <- env_undirected |> tidygraph::convert(tidygraph::to_simple)
        # layout <- ggraph::create_layout(env_undirected, "backbone")
        layout <- ggraph::create_layout(
          env_undirected,
          "igraph", algorithm = "bipartite",
          types = targetPoolsU[[i]]$Pool$Type[
            env_undirected %N>% tidygraph::pull(
              node
            ) |> substring(2) |> as.numeric()
            ] == "Basal")
      } else {
        layout <- ggraph::create_layout(env_undirected, "auto")
      }

      return(layout)
    }
  )

  targetEnvs <- lapply(seq_along(targetEnvs), function(i, envs, grps, layouts) {
    envs[[i]]$layout <- layouts[[grps[i]]]
    return(envs[[i]])
  },
  envs = targetEnvs, grps = targetPoolGroups, layouts = targetEnvsLayouts
  )

  # Adjust for meaningful x and y instead if possible.
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
    l$y <- log10(env$Pool$Pool$Size[l_indices])
    env$layout <- l
    return(env)
  })

  #### Plot
  targetEnvs <- lapply(targetEnvs, function(env) {
    timelist <- env$graphs
    env$singletonGraphs <- lapply(env$graphs, function(patchgraph) {
      plotGraph(patchgraph, mainLayout = env$layout)
    })
    return(env)
  })

  targetEnvsIndex <- do.call(
    rbind,
    lapply(
      targetEnvs,
      function(env)
        cbind(
          Time = env$Row$Time,
          env$Diversity[
            1, c("PoolPatchSeed", "SpeciesAffinity",
                 "NicheDistance", "Intervention")]
        )
    )
  )
  targetEnvsIndex <- cbind(ID = 1:nrow(targetEnvsIndex), targetEnvsIndex)

  return(list(
    Index = targetEnvsIndex,
    Envs = targetEnvs
  ))
}
