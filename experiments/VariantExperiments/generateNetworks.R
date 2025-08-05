generateNetworks <- function(
  Specification,
  dIS = diversitiesInterventionStrings,
  aDO = affinityDictionaryOrigin
) {
  # Goal: Pass data.frame, return list of networks where a network corresponds
  #       to a row of the data.frame.
  # NOTE: This will need to be updated for Version 10. Currently Version 9.
  # Singleton version appears to be working fine. Need to deal with > 1.
  # TODO: TEST.
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
  targetDirs <- with(Specification, dir(
    pattern = paste0(
      "TSTS_Simulations_",
      PoolPatch, "-", Interactions, "_",
      PoolPatchSeed, "-", InteractionsSeed, "_"
    ),
    full.names = TRUE
  ))

  targetFiles <- dir(targetDirs, pattern = "(Sim|Int)", full.names = T)

  ### Retrieval Base: ##########################################################
  targetFiles <- with(Specification, grep(
    x = targetFiles,
    pattern = paste0(
      PoolPatch, "-", Interactions, "-", Events, "-",
      InitialConditions, "-", Dispersal, "-", NicheDistance, "-", Affinity, "_",
      PoolPatchSeed, "-", InteractionsSeed, "-", EventsSeed, "-",
      InitialConditionsSeed, "-", AffinitySeed, "_"
    ),
    value = TRUE
  ))
  targetFilesSim <- grep(
    x = targetFiles,
    pattern = "_Simulation_",
    fixed = TRUE,
    value = TRUE
  ) # Done.

  ### Retrieval Intervention: ##################################################
  targetFilesInt <- grep(
    x = targetFiles,
    pattern = "Int",
    fixed = TRUE,
    value = TRUE
  ) # Check the right intervention asked for.

  targetFilesInt <- grep(
    x = targetFilesInt,
    pattern = with(Specification, paste0(
      InterventionPatchType, "-", InterventionTimeType, "-",
      InterventionDispersal, "-", InterventionNicheDistance, "_",
      InterventionPatchSeed, "-", InterventionTimeSeed
    )),
    fixed = TRUE,
    value = TRUE
  )

  # Families of targets:
  targetFiles <- c(targetFilesSim, targetFilesInt)
  targetDivs <- gsub(x = targetFiles,
                     pattern = "_(Simulation|Intervention)_",
                     replacement = "_Diversity_")
  targetPools <- dir(targetDirs, pattern = "Pool",
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
      PoolPatchSeed = grep(div$Ellipsis$ID,
                           pattern = "(?<=_)[0-9]+", value = TRUE, perl = TRUE),
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

  # Implement Associations
  # Can make more memory efficient by reducing at this step just to the needed
  # temporal slices.
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


  targetEnvs <- lapply(targetEnvs, function(env) {
    # print(env$graphs[[1]]$graphs[[1]])
    env_undirected <- lapply(env$graphs, tidygraph::to_undirected)
    if (length(env_undirected) > 1) {
      for (i in 2:length(env_undirected)) {
        env_undirected[[1]] <- tidygraph::graph_join(
          env_undirected[[1]], env_undirected[[i]]
        )
      }
    }
    if ((env_undirected[[1]] %>% mutate(
      Components = tidygraph::graph_component_count(),
      Nodes = tidygraph::graph_order(),
      NotSingletons = Components != Nodes
    ) %>% pull(NotSingletons))[1]
    ) {
      env_undirected[[1]] <- env_undirected[[1]] |> tidygraph::convert(tidygraph::to_simple)
      env$layout <- ggraph::create_layout(env_undirected[[1]], "backbone")
    } else {
      env$layout <- ggraph::create_layout(env_undirected[[1]], "auto")
    }
    return(env)
  })

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
        env$Diversity[
          1, c("PoolPatchSeed", "SpeciesAffinity",
               "NicheDistance", "Intervention")]
    )
  )
  targetEnvsIndex <- cbind(ID = 1:nrow(targetEnvsIndex), targetEnvsIndex)

  return(list(
    Index = targetEnvsIndex,
    Envs = targetEnvs
  ))
}
