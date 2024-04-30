# I've been asked to create visualisations of food webs for each environment
# through time. This is non-trivial.
# At present, my idea is to compute the "maximal" food web from the data first.
# This way, each species has a well identified position.
# One awkward thing to resolve is if it is easier to have the species on some
# form of grid or if it is better to have lines connecting them between patches.
# Regardless, we can then remove (or maybe dim?) species that are not present in
# a patch.

# In an ideal world, I think it would like an arrangement of foodwebs with a
# timeline thatcan be "played" like a movie forwards and paused. When paused,
# one could click on a node to see its sources of gains and losses.
# For a line or ring of ten patches:
#  1 - 2 - 3 - 4 - 5
# 10 - 9 - 8 - 7 - 6
# to keep patches next to each other.
# For 2 patches, obviously 1 - 2.

# I'm thinking dim on the grid is best bet.
# #############################################################################

# Problems with X11
options(bitmapType = "cairo")

systype <- match.arg(
  "Simulation",
  c("Simulation", "Intervention")
)

# Intervention specific version
# LOAD POOLPATCHDYNAMICS
# LOAD INTERVENTION

set <- "18-1"
tag <- paste0(set, "-2-25-3_6-6-6",
              if (systype == "Intervention") "_40-1-15-3_1-1")
dir <- "TSTS_Simulations_18-1_6-6_2024-03-08"
load(file.path(dir, paste0("TSTS_",systype, "_", tag, ".RData")))
load(file.path(dir, paste0("TSTS_PoolPatchDynamics_", set, ".RData")))
load(file.path(dir, paste0("TSTS_Diversity_", tag, ".RData")))
load(file.path(dir, paste0("TSTS_DivAbund_", tag, ".RData")))

library(dplyr)

threshold <- 0.00 # in [0, 1]

if (systype == "Intervention") {
  timesInUse <-
    result$Abundance[, 1] >= result$Ellipsis$Affinity$TimeIntervention

  times <- result$Abundance[timesInUse, 1]
} else {
  timesInUse <- rep(TRUE, nrow(result$Abundance))
  times <- result$Abundance[, 1]
}
jaccard <- lapply(Diversity$Diversities$beta, as.data.frame)
jaccard <- lapply(jaccard, function(x) {x$Jaccard <- as.numeric(x$Jaccard); x})
jaccard <- jaccard %>% dplyr::bind_rows() %>% dplyr::filter(
  Env1 == 1#, Time > result$Ellipsis$Affinity$TimeIntervention
)
braycurtis <- DivAbund$DivAbund$beta
braycurtis$BrayCurtis <- as.numeric(braycurtis$BrayCurtis)
braycurtis <- braycurtis %>% dplyr::filter(
  Env1 == 1#, Time > result$Ellipsis$Affinity$TimeIntervention
)
if (result$Ellipsis$Timescale == "Simulation") {
  times <- times / result$ReactionTime
  #jaccard$Time <- jaccard$Time / result$ReactionTime
  #braycurtis$Time <- braycurtis$Time / result$ReactionTime
}


# Nodes: ######################################################################
# Esp. Node Sizes: (Essentially the *red*uced *Com*munity from earlier.)
environs <- lapply(
  1:result$NumEnvironments,
  function(i, ns) {
    ### Abundance and Interaction Matrices: ###################################
    # All candidates for an environment.
    indices <- 1 + 1:ns + (i - 1) * ns

    # Reduce to the appropriate submatrix.
    temp <- result$Abundance[timesInUse, indices]

    # Keep only columns that have species present at some point in the sim.
    keep <- apply(temp, 2,
                  function(x) any(x > result$Parameters$EliminationThreshold))
    temp <- temp[, keep]

    # Fix any "nanofoxes".
    temp <- ifelse(temp < result$Parameters$EliminationThreshold,
                   0, temp)

    # Preserve naming.
    colnames(temp) <- which(keep)

    # Apply to the interaction matrix
    interactions <- InteractionMatrices$Mats[[i]][keep, keep]
    colnames(interactions) <- which(keep)
    rownames(interactions) <- which(keep)

    ### Recreate Dispersal Matrix: ############################################
    if (systype == "Intervention") {
      dispersalDictionaryChoice <-
        as.numeric((
          (strsplit(result$Ellipsis$ID, "_")[[1]][3]) %>% strsplit("-")
          )[[1]][3])
    } else if (systype == "Simulation") {
      dispersalDictionaryChoice <-
        as.numeric((
          (strsplit(result$Ellipsis$ID, "_")[[1]][1]) %>% strsplit("-")
        )[[1]][4])
    }

    dispersalDictionary <- rbind(
      data.frame(Resistance = Inf, Configuration = "None"),
      expand.grid(
        Resistance = 10^c(0:9),
        Configuration = c("Ring", "Line", "Complete")
      ))[ifelse(is.na(dispersalDictionaryChoice),
                1, dispersalDictionaryChoice + 2), ]

    DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
      EnvironmentDistances = with(c(
        dispersalDictionary,
        Environments = result$NumEnvironments
      ), {
        if (Configuration == "None") {
          DistanceMatrix <- Matrix::sparseMatrix(
            i = Environments, j = Environments, x = 0)
        }
        if (Configuration == "Ring" || Configuration == "Line")
          DistanceMatrix <- Matrix::bandSparse(
            Environments, k = c(-1, 1),
            diagonals = list(rep(Resistance, Environments - 1),
                             rep(Resistance, Environments - 1))
          )
        if (Configuration == "Ring") {
          DistanceMatrix[Environments, 1] <- Resistance
          DistanceMatrix[1, Environments] <- Resistance
        }
        if (Configuration == "Complete") {
          DistanceMatrix <- matrix(Resistance,
                                   nrow = Environments,
                                   ncol = Environments)
          diag(DistanceMatrix) <- 0
        }
        DistanceMatrix
      }
      ),
      SpeciesSpeeds = Pool$Speed
    )

    In <- DispersalMatrix
    In[In < 0] <- 0
    Out <- DispersalMatrix
    Out[Out > 0] <- 0

    ### Return: ###############################################################
    return(list(
      Abundance = temp,
      Size = Pool$Size[keep],
      Type = Pool$Type[keep],
      Affinity = Pool$Affinity[keep],
      Matrix = interactions,
      DispersalGain = In[indices[keep] - 1, ], # -1 since time not in Dispersal
      DispersalLoss = Out[indices[keep] - 1, ],
      Intrinsic =
        if (systype == "Intervention") {
        result$Ellipsis$Affinity$EffectiveReproductionRateIntervention[
          indices[keep] - 1 # since time not in intrinsic rep rate.
          ]
        } else if (systype == "Simulation") {
          result$Ellipsis$Affinity$EffectiveReproductionRate[indices[keep] - 1]
        }
    ))
  }, ns = (ncol(result$Abundance) - 1) / result$NumEnvironments)

### Adjust for species present in other patches: ##############################

full <- sort(unique(as.numeric(
  unlist(lapply(environs, function(x) {rownames(x$Matrix)}))
)))

environs <- lapply(environs, function(x) {
  Missing <- full[!full %in% as.numeric(rownames(x$Matrix))]
  x$Missing <- Missing
  return(x)
})

# Each Row is then a set of nodes at a specific time.
# We'll generate each plot from each row, although we may want to consider the
# timing of events, but checking plot(temp[, 1]) suggests roughly linearity.
# plot(temp[, 1])
# lines(seq(from = temp[1, 1],
#           to = temp[nrow(temp), 1],
#           length.out = nrow(temp)), col = "red")


# Plotting: ###################################################################
# Game plan roughly from:
# https://archive.schochastics.net/post/ggraph-tricks-for-common-problems/
# We'll be creating a layout (of layouts) and create a graph of plots.
# Edges are interactions with other species.
# We'll also want edges from/to nowhere for intrinsic values.
# After all of this, we'll then want to make a video from it.

# We'll want to duplicate the layout and then offset them.

library(tidygraph)
library(ggraph)
library(ggpubr)
library(animation)

### Helper Functions: #########################################################
# These need to be called multiple times do to have (nontrivial) layouts.
createBaseGraph <- function(timestep, environs) {
  lapply(
    environs, function(e) {
      g <- tidygraph::as_tbl_graph(
        # Careful here: obvious way is to go from %*% to * to avoid addition
        #               in the matrix-vector product.
        #               R does its multiplications column wise though!
        #               Hence transpose t().
        #               Not doubled because tidygraph uses the opposite
        #               from-to convention from me.
        #               colSums(...) returns the correct values.
        #               (= (InteractionMatrices$Mats[[1]] %*%
        #                  (result$Abundance[timestepResult, -1][1:200]))[
        #                   as.numeric(colnames(environs[[1]]$Abundance))])
        t(e$Matrix) * e$Abundance[timestep, ]
      ) %>% tidygraph::mutate(
        Present = e$Abundance[timestep, ] > 0,
        Abundance = e$Abundance[timestep, ],
        Size = e$Size,
        Type = e$Type,
        Affinity = e$Affinity
      )

      if (sum(e$Abundance[timestep, ] > 0) > 1) {
        g <- g %>% tidygraph::activate(edges) %>% tidygraph::mutate(
          Type = ifelse(weight > 0, "Consumption", ifelse(
            to == from, "Intraspecific", "Predation"
          ))
        ) %>% tidygraph::activate(nodes)
      }

      return(g)
    }
  )
}

addIntrinsicToGraph <- function(graph, timestep, environs) {
  #KEY:
  #   ON if ever PRESENT
  #   PRESENT if currently have non-zero abundance
  #   SPECIES if you represent a real species.
  lapply(seq_along(environs), function(i, g, e) {
    # Add missing nodes
    temp <- tidygraph::bind_graphs(
      g[[i]] %>% dplyr::mutate(On = TRUE, Species = TRUE),
      tidygraph::tbl_graph(
        data.frame(name = as.character(e[[i]]$Missing))
      ) %>% dplyr::mutate(
        On = FALSE, Species = TRUE, Present = FALSE, Abundance = 0,
        Size = Pool$Size[e[[i]]$Missing],
        Affinity = Pool$Affinity[e[[i]]$Missing]
      )
    ) %>% tidygraph::arrange(name)

    # Add fake nodes
    temp <- tidygraph::bind_graphs(
      temp,
      tidygraph::tbl_graph(
        data.frame(name = paste0(".", temp %>% dplyr::pull(name)))
      ) %>% dplyr::mutate(
        On = FALSE, Species = FALSE, Present = FALSE, Abundance = 0,
        Size = Pool$Size[temp %>% dplyr::pull(name) %>% as.numeric],
        Affinity = NA
      )
    )

    # Add intrinsic rates as edges to/from the fake nodes.
    temp <- temp %>% tidygraph::bind_edges(
      data.frame(
        to = ifelse(e[[i]]$Intrinsic > 0,
                    rownames(e[[i]]$Matrix),
                    paste0(".", rownames(e[[i]]$Matrix))),
        from = ifelse(e[[i]]$Intrinsic > 0,
                      paste0(".", rownames(e[[i]]$Matrix)),
                      rownames(e[[i]]$Matrix)),
        weight = e[[i]]$Intrinsic,
        Type = ifelse(e[[i]]$Intrinsic > 0, "Growth", "Decay")
      ), node_key = "name"
    )

    # Add dispersal rates as edges to/from the fake nodes.
    # Note units need to be per capita, so we calculate as:
    #  (Dispersal %*% Abundance) / (Local Abundance)
    temp <- temp %>% tidygraph::bind_edges(
      data.frame(
        to = rownames(e[[i]]$Matrix),
        from = paste0(".", rownames(e[[i]]$Matrix)),
        weight = (e[[i]]$DispersalGain %*%
                    result$Abundance[timesInUse, ][timestep, -1])[, 1] /
          e[[i]]$Abundance[timestep, ],
        Type = "Dispersal"
      ) %>% dplyr::filter(!is.nan(weight),
                          weight != 0,
                          !is.infinite(weight)), node_key = "name"
    )

    temp <- temp %>% tidygraph::bind_edges(
      data.frame(
        to = paste0(".", rownames(e[[i]]$Matrix)),
        from = rownames(e[[i]]$Matrix),
        weight = (e[[i]]$DispersalLoss %*%
                    result$Abundance[timesInUse, ][timestep, -1])[, 1] /
          e[[i]]$Abundance[timestep, ],
        Type = "Dispersal"
      ) %>% dplyr::filter(!is.nan(weight),
                          weight != 0,
                          !is.infinite(weight)), node_key = "name"
    )

    return(temp)
  }, g = graph, e = environs)
}

thresholdGraphEdges <- function(graph, threshold) {
  lapply(graph, function(g) {
    g %>% tidygraph::activate(edges) %>% tidygraph::group_by(
      to
    ) %>% tidygraph::mutate(
      PercentContribution = abs(weight) / sum(abs(weight)),
      DestinationPresent = .N()$Present[to] | !.N()$Species[to], # Allow Decays
      SourcePresent = .N()$Present[from] | !.N()$Species[from], # Allow Growths
      Linetype = ifelse(weight > 0, "Positive", "Negative")
    ) %>% tidygraph::filter(
      PercentContribution > threshold,
      DestinationPresent
    ) %>% tidygraph::ungroup()
  })
}

computeGraphLayout <- function(graph) {
  grafNonEmpty <- graph[[which.max(
    unlist(lapply(graph, function(g) {g %N>% filter(
      Species
    ) %E>% tidygraph::filter(
      Type != "Predation",
      Type != "Intraspecific",
      Type != "Dispersal"
    ) %E>% pull(to) %>% length}))
  )]]
  lay <- ggraph::create_layout(
    tidygraph::to_undirected(
      grafNonEmpty %N>% filter(
        Species
      ) %E>% tidygraph::filter(
        Type != "Predation",
        Type != "Intraspecific",
        Type != "Dispersal"
      )
    ) ,
    "backbone"
    # "stress", y = log10(Size)
  )
  lay$y <- log10(lay$Size)
  return(lay)

  # # Create a complete graph between all real consumers,
  # # Joined to all consumers and all basals,
  # # Joined to all real nodes with their fake equivalents.
  # nodesN <- grep(graf[[1]] %N>% dplyr::pull(name), pattern = ".",
  #                fixed = TRUE, value = TRUE, invert = TRUE)
  # nodeTypes <- Pool$Type[as.numeric(nodesN)]
  # nodesTable <- table(nodeTypes)
  #
  # nodesConsumers <- create_complete(
  #   nodesTable[2]
  # ) %N>% mutate(
  #   name = nodesN[nodeTypes == "Consumer"]
  # )
  #
  # nodesBipartite <- create_bipartite(
  #   nodesTable[1], nodesTable[2]
  # ) %N>% mutate(
  #   Type = ifelse(!type, "Basal", "Consumer")
  # ) %>% group_by(
  #   Type
  # ) %>% mutate(
  #   name = nodesN[nodeTypes == Type],
  #   Size = Pool$Size[as.numeric(name)]
  # ) %>% select(
  #   -type
  # )
  #
  # # nodesPairs <- lapply(1:length(nodesN), function(i) {
  # #   create_complete(2) %N>% mutate(
  # #     name = c(nodesN[i], paste0(".", nodesN[i])),
  # #     Size = Pool$Size[as.numeric(nodesN[i])]
  # #   )
  # # }) %>% tidygraph::bind_graphs()
  #
  # nodesAll <- tidygraph::graph_join(
  #   nodesConsumers, nodesBipartite, by = "name"
  # # ) %>% tidygraph::graph_join(
  # #   nodesPairs, by = c("name", "Size")
  # ) %>% tidygraph::to_undirected()
  #
  # lay <- ggraph::create_layout(
  #   nodesAll,
  #   # "backbone"
  #   # "stress", y = log10(Size),
  #   # "auto", #y = log10(Size)
  #   # "unrooted"
  #   # "igraph", algorithm = "graphopt"
  #   "linear"
  # )
  # lay$y <- log10(lay$Size)
  # ggraph(
  #   nodesAll, layout = data.frame(lay[, 1:2])
  # ) + ggraph::geom_edge_diagonal(
  #   alpha = 0.01
  # ) + ggraph::geom_node_label(
  #   mapping = aes(color = Type, label = name)
  # )
  # # Might need to do a linear layout of consumers,
  # # tack on a stress layout with the basals,
  # # then manually add all of the fake nodes to the side.
}

anyAbundanceSoFar <- FALSE
largestTotalEffect <- 0

# Create a common layout.
timestepLayout <- which.max(rowSums(result$Abundance[timesInUse, -1] > 0))
graf4Layout <- createBaseGraph(timestepLayout, environs)
graf4Layout <- addIntrinsicToGraph(graf4Layout, timestepLayout, environs)
graf4Layout <- thresholdGraphEdges(graf4Layout, threshold)
lay <- computeGraphLayout(graf4Layout)

#
#
# # GIF
# animation::saveGIF(
#   movie.name = "test.gif",
# Video
animation::saveVideo(
  # video.name = "test.mp4",
  video.name = paste0("test_", tag, ".mp4"),
  ani.height = 1024 * 2, ani.width = 1280 * 2, interval = 0.1,
  expr = {
    for (timestep in seq(from = 1, to = nrow(environs[[1]]$Abundance),
                         # length.out = nrow(result$Events)/1000
                         length.out = nrow(result$Events)/1
    )) {
      # ani.options(ani.height = 1280 * 2, ani.width = 1024 * 2, interval = 1)

      if (!anyAbundanceSoFar) {
        if (!any(result$Abundance[timesInUse,][timestep, -1] > 0)) {
          next
        } else {
          anyAbundanceSoFar <- TRUE
        }
      }

      timestep <- round(timestep)

      graf <- createBaseGraph(timestep, environs)
      graf <- addIntrinsicToGraph(graf, timestep, environs)
      graf <- thresholdGraphEdges(graf, threshold)


      edgecolors <- c(
        # to colorbrewer2.org!
        # Extrinsics overlaid, so need to be different colors
        # Extrinsics can be dark, intrinsics can be light.
        "Consumption" = "#2c7bb6", # Positive, Extrinsic
        "Predation" = "#d7191c", # Negative, Extrinsic
        "Intraspecific" = "#ffffbf",# Negative Extrinsic Competition
        "Growth" = "#abd9e9", # Positive, Intrinsic
        "Decay" = "#fdae61", # Negative, Intrinsic
        "Dispersal" = "#000000"
      )

      graf <- lapply(graf, function(g) {
        df <- g %>% activate(
          edges
        ) %>% igraph::as_data_frame(
        ) %>% dplyr::mutate(
          Node = ifelse(
            # Not Decay or Negative Dispersal, use To, else From.
            !Type == "Decay" & !(Type == "Dispersal" & Linetype == "Negative"),
            to, from)
        ) %>% dplyr::filter(
          Node %in% (g %>% activate(nodes) %>% filter(Present) %>% pull(name))
        ) %>% dplyr::group_by(
          Node#,Linetype, Type
        ) %>% dplyr::summarise(
          TotalEffect = sum(weight),
          .groups = "drop"
        )

        fills <- df %>% dplyr::group_by(Node) %>% dplyr::group_modify(
          .f = function(.x, .y) {
            # .x$Type[which.max(abs(.x$TotalEffect))]
            data.frame(
              fill = .x$TotalEffect
            )
          }
        )

        g %N>% left_join(fills, by = c("name" = "Node"))

      })

      candidateLargestTotalEffect <-
        max(abs(unlist(lapply(graf, function(g) g %N>% pull(fill)))))
      if (largestTotalEffect < candidateLargestTotalEffect) {
        largestTotalEffect <- candidateLargestTotalEffect
        print(largestTotalEffect)
      }


      plots <- lapply(seq_along(graf), function(i) {
        g <- graf[[i]] %N>% filter(Species)

        if (i != 1) {
          corWith1 <- (full_join(
            graf[[1]] %N>% filter(Abundance > 0) %>% as.data.frame(
            ) %>% select(name, Abundance),
            graf[[i]] %N>% filter(Abundance > 0) %>% as.data.frame(
            ) %>% select(name, Abundance),
            by = "name", copy = TRUE) %>% select(
              -name
            ) %>% mutate(
              across(starts_with("Abundance"), .fns = ~ifelse(is.na(.x), 0, .x))
            ) %>% cor())[1, 2]# %>% signif(digits = 2)


          jacWith1 <- jaccard %>% dplyr::filter(
            Env2 == i
          )
          jacWith1 <- jacWith1[
            which.min(abs(jacWith1$Time - times[timestep])),
            ]$Jaccard #%>% signif(digits = 4)


          bcWith1 <- braycurtis %>% dplyr::filter(
            Env2 == i
          )
          bcWith1 <- bcWith1[
            which.min(abs(bcWith1$Time - times[timestep])),
            ]$BrayCurtis #%>% signif(digits = 4)
        }

        lastSuccessfulNeutralExtirpation <- result$Events
        if (result$Ellipsis$Timescale == "Simulation") {
          lastSuccessfulNeutralExtirpation <-
            lastSuccessfulNeutralExtirpation %>% dplyr::mutate(
              Times = Times / result$ReactionTime
            )
        }
        lastSuccessfulNeutralExtirpation <-
          lastSuccessfulNeutralExtirpation  %>% dplyr::filter(
            Success, Times <= times[timestep], Type == "Extinct"
          )
        lastSuccessfulNeutralExtirpation <-
          lastSuccessfulNeutralExtirpation[
            nrow(lastSuccessfulNeutralExtirpation),
            ]


        ggraph::ggraph(
          g, layout = data.frame(lay[, 1:2])
        ) + ggraph::geom_node_point(
          mapping = aes(
            color = Type, alpha = interaction(Present, On),
            stroke = log1p(Abundance), fill = fill, shape = factor(Affinity)
          ),
          #shape = 22,
          size = 12
        ) + ggraph::geom_edge_diagonal(
          mapping = aes(
            color = Type, filter = SourcePresent,
            linetype = Linetype, alpha = PercentContribution
          ),
          linewidth = 2,
          arrow = arrow(length = unit(2, 'mm')),
          start_cap = circle(5, 'mm'),
          end_cap = circle(5, 'mm')
        ) + scale_shape_manual(
          values = c(22, 23), name = "b:"
        ) + scale_alpha_manual(
          values = c("FALSE.FALSE" = 0, "FALSE.TRUE" = 0.25, "TRUE.TRUE" = 1),
          guide = "none"
        ) + scale_color_manual(
          name = "Node",
          values = c("Basal" = "darkgreen", "Consumer" = "yellow3")
        ) + scale_fill_gradient2(
          low = "#FF00FF", mid = "#FFFFFF", high = "#00FFFF",
          limits = c(-6, 6), breaks = c(-4,0,4),
          transform = scales::transform_pseudo_log(sigma = 1/10),
          name = "Per. Cap. Rate"
        ) + scale_edge_linetype_manual(
          name = "",
          values = c("Positive" = "dashed", "Negative" = "longdash")
        ) + scale_edge_color_manual(
          name = "Edge",
          values = edgecolors, drop = FALSE
        ) + scale_edge_alpha(
          name = "Contribution %",
          transform = "sqrt", breaks = c(1e-5, 1e-3, 1e-1), limits = c(1e-7, 1)
        ) + theme_minimal(
        ) + ylab(
          "Log10(Size)"
        ) + theme(
          text = ggplot2::element_text(size = 40),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        ) + ggtitle(
          paste("P:", i,
                "b:", if(i %in% result$Ellipsis$Affinity$PatchInterventions) {
                  result$Ellipsis$Affinity$PatchAffinitiesIntervention[i,]
                } else {
                  result$Ellipsis$Affinity$PatchAffinitiesOld[i,]
                },
                if(i == 1) {paste(
                  "ts:", timestep, "T:",
                  round(times[timestep], digits = 1),
                  "Last N. Ext.:",
                  round(lastSuccessfulNeutralExtirpation$Times,
                        digits = 1))},
                if(i != 1) {paste(
                  c("R:", "Jac:", "BC:"),
                  unlist(lapply(c(corWith1, jacWith1, bcWith1), function(x) {
                    if (x == 1) {
                      " 1         "
                    } else if (x == 0) {
                      " 0         "
                    } else if (x > 0) {
                      paste0(" ", formatC(x,preserve.width = "common",
                                          digits = 4, format = "f"))
                    } else {
                      formatC(x,preserve.width = "common",
                              digits = 4, format = "f")
                    }
                  })),
                  collapse = ", "
                )}
          )
        )
      })

      # # https://archive.schochastics.net/post/ggraph-tricks-for-common-problems/
      # pointsAsBarPlots <- lapply(graf, function(g) {
      #   df <- g %>% activate(
      #     edges
      #   ) %>% igraph::as_data_frame(
      #   ) %>% dplyr::mutate(
      #     Node = ifelse(
      #       # Not Decay or Negative Dispersal, use To, else From.
      #       !Type == "Decay" & !(Type == "Dispersal" & Linetype == "Negative"),
      #       to, from)
      #   ) %>% dplyr::filter(
      #     Node %in% (g %>% activate(nodes) %>% filter(Present) %>% pull(name))
      #   ) %>% dplyr::group_by(
      #     Node, Type, Linetype
      #   ) %>% dplyr::summarise(
      #     TotalEffect = sum(weight),
      #     .groups = "drop"
      #   )
      #
      #   limits <- range(df$TotalEffect)
      #
      #   df %>% dplyr::group_by(Node) %>% dplyr::group_map(
      #     .f = function(.x, .y) {
      #       grob <- ggplot2::ggplotGrob(
      #         ggplot2::ggplot(
      #           .x
      #         ) + ggplot2::geom_col(
      #           ggplot2::aes(
      #             x = interaction(Type, Linetype),
      #             y = TotalEffect,
      #             fill = Type
      #           )
      #         ) + ggplot2::scale_fill_manual(
      #           values = edgecolors
      #         ) + ggplot2::coord_cartesian(
      #           ylim = limits,
      #           expand = FALSE
      #         ) + ggplot2::theme(
      #           legend.position = "none",
      #           panel.background =
      #             ggplot2::element_rect(fill = "white", color = NA),
      #           line = ggplot2::element_blank(),
      #           text = ggplot2::element_blank()
      #         )
      #       )
      #       coordinates <- grob$layout[grob$layout$name == "panel", ]
      #       list(
      #         node = .y$Node,
      #         grobs = grob[coordinates$t:coordinates$b,
      #                      coordinates$l:coordinates$r]
      #       )
      #     }
      #   )
      # })
      #
      # annotationspace <- 0.14
      #
      # barplotsAsAnnotations <- lapply(pointsAsBarPlots, function(g) {
      #   lapply(g, function(bp) {
      #     node <- bp$node
      #     grob <- bp$grobs
      #     xy <- lay %>% dplyr::filter(name == node) %>% dplyr::select(x, y)
      #     ggplot2::annotation_custom(
      #       grob,
      #       xmin = xy$x - annotationspace * 2,
      #       xmax = xy$x + annotationspace * 2,
      #       ymin = xy$y - annotationspace / 10,
      #       ymax = xy$y + annotationspace / 10
      #     )
      #   })
      # })
      #
      # combinedplots <- lapply(seq_along(graf), function(i) {
      #   Reduce("+", barplotsAsAnnotations[[i]], plots[[i]])
      # })

      tempplot <- ggpubr::ggarrange(
        plotlist = plots, #combinedplots,
        ncol = 2, nrow = 1,
        common.legend = TRUE,
        legend = "bottom"
      )
      plot(tempplot)
    }
  }
)
