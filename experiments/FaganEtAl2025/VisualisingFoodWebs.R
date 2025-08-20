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

library(dplyr)

# Problems with X11
options(bitmapType = "cairo")

systype <- match.arg(
  "Intervention",
  c("Simulation", "Intervention")
)

directory <- '.' # Assumed to be VariantExperiments
source(file.path(directory, "TimeSpaceAndTimeSeries-9-Dictionaries.R"))

# E.g. TSTS_Diversity_142486-4929-28-1-NA-3-1_341-341-384-387-391
# or TSTS_Diversity_142486-4929-28-1-NA-3-1_341-341-384-387-391_111-1-p-p_1-1
# TSTS_Diversity_142486-4929-28-1-NA-7-35_341-341-384-387-410
set <- "142486-4929"; initSeeds <- "341-341"; date <- "2024-11-30"
tag <- paste0(set, "-28-1-NA-7-7_", initSeeds, "-384-387-394",
              if (systype == "Intervention") "_115-1-p-p_1-1") # usually 111:115
dir <- paste0("TSTS_Simulations_", set, "_", initSeeds, "_", date)
load(file.path(dir, paste0("TSTS_", systype, "_", tag, ".RData")))
load(file.path(dir, paste0("TSTS_PoolPatchDynamics_", set, "_", initSeeds, ".RData")))
load(file.path(dir, paste0("TSTS_Diversity_", tag, ".RData")))
# load(file.path(dir, paste0("TSTS_DivAbund_", tag, ".RData"))) # Deprecated.



threshold <- 1e-7 # in [0, 1]

if (systype == "Intervention") {
  timesInUse <-
    result$Abundance[, 1] >= result$Ellipsis$Affinity$TimeIntervention

  times <- result$Abundance[timesInUse, 1]
} else {
  timesInUse <- rep(TRUE, nrow(result$Abundance))
  times <- result$Abundance[, 1]
}

# Currently working in 1 patch environments, so I can't test how I might want to
# adapt these.
# jaccard <- lapply(Diversity$Diversities$beta, as.data.frame)
# jaccard <- lapply(jaccard, function(x) {x$Jaccard <- as.numeric(x$Jaccard); x})
# jaccard <- jaccard %>% dplyr::bind_rows() %>% dplyr::filter(
#   Env1 == 1#, Time > result$Ellipsis$Affinity$TimeIntervention
# )
# braycurtis <- DivAbund$DivAbund$beta
# braycurtis$BrayCurtis <- as.numeric(braycurtis$BrayCurtis)
# braycurtis <- braycurtis %>% dplyr::filter(
#   Env1 == 1#, Time > result$Ellipsis$Affinity$TimeIntervention
# )

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
      dispersalDictionaryChoice <-(
        (strsplit(result$Ellipsis$ID, "_")[[1]][3]) %>% strsplit("-")
      )[[1]][3]

      if (dispersalDictionaryChoice == "p") {
        # Invoke Previous.
        dispersalDictionaryChoice <-
          as.numeric((
            (strsplit(result$Ellipsis$ID, "_")[[1]][1]) %>% strsplit("-")
          )[[1]][5])

      } else if (dispersalDictionaryChoice == "NA") {
        # True NAs dealt with
        dispersalDictionaryChoice <- NA

      } else {
        # Convert to Numeric and make sure conversion worked.
        dispersalDictionaryChoice <- tryCatch({
          as.numeric(dispersalDictionaryChoice)
        }, error = function(e) {return(e)})
        stopifnot(!is.na(dispersalDictionaryChoice) || # stop if (false) NA
                    is.numeric(dispersalDictionaryChoice)) # stop if not numeric
      }
    } else if (systype == "Simulation") {
      dispersalDictionaryChoice <-
        as.numeric((
          (strsplit(result$Ellipsis$ID, "_")[[1]][1]) %>% strsplit("-")
        )[[1]][5])
    } else {
      stop("systype not detected/understood when retrieving dispersal matrix.")
    }

    dispersalDictionary <-
      dispersalDictionaryOrigin[ifelse(is.na(dispersalDictionaryChoice),
                                       1, dispersalDictionaryChoice + 2), ]

    if (result$NumEnvironments > 1) {
      DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
        EnvironmentDistances = convertDispersalDictToDistMatrix(
          dispersalDictionary,
          nEnv = result$NumEnvironments
        ),
        SpeciesSpeeds = Pool$Speed
      )
    } else {
      DispersalMatrix <- Matrix::sparseMatrix(
        i = {}, j = {}, # From documentation
        dims = c(nrow(Pool), nrow(Pool))
      )
    }

    In <- DispersalMatrix
    In[In < 0] <- 0
    Out <- DispersalMatrix
    Out[Out > 0] <- 0

    ### Return: ###############################################################
    return(list(
      Abundance = temp,
      Size = Pool$Size[keep],
      Type = Pool$Type[keep],
      Affinity = result$Ellipsis$Affinity$SpeciesAffinities[keep],
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
# timing of events, but checking time plot suggests roughly linearity.
# plot(result$Abundance[, 1])
# lines(seq(from = result$Abundance[1, 1],
#           to = result$Abundance[nrow(result$Abundance), 1],
#           length.out = nrow(result$Abundance)), col = "red")


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
# These need to be called multiple times so to have (nontrivial) layouts.
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

      if (sum(e$Abundance[timestep, ] > 0) >= 1) { # If any abundance...
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

patchBackgrounds <- c("#e66101", "#5e3c99")
plotEventRug <- function(gplot) {
  yrange <- ggplot2::layer_scales(gplot)$y$range$range
  xrange <- ggplot2::layer_scales(gplot)$x$range$range
  gplot + ggplot2::geom_rect(
    data = result$Events %>% dplyr::filter(Success) %>% dplyr::mutate(
      Environment = as.factor(Environment), Times = Times / result$ReactionTime,
      NextEventTime = lead(Times)
    ),
    mapping = ggplot2::aes(
      xmin = Times, xmax = NextEventTime, fill = Environment
    ),
    ymin = -Inf, ymax = Inf,
    alpha = 0.15, inherit.aes = FALSE
  ) + ggplot2::geom_rug(
    data = result$Events %>% dplyr::filter(
      Success, Type == "Extinct"
    ) %>% dplyr::mutate(
      Environment = as.factor(Environment), Times = Times / result$ReactionTime
    ),
    mapping = ggplot2::aes(
      x = Times, color = Environment
    ), inherit.aes = FALSE,
    sides = "b"
  ) + ggplot2::geom_rug(
    data = result$Events %>% dplyr::filter(
      Success, Type == "Arrival"
    ) %>% dplyr::mutate(
      Environment = as.factor(Environment), Times = Times / result$ReactionTime
    ),
    mapping = ggplot2::aes(
      x = Times, color = Environment
    ), inherit.aes = FALSE,
    sides = "t"
  ) + ggplot2::theme_bw() + ggplot2::annotate(
    "text",
    x = xrange[1] + diff(xrange)/20,
    y = yrange + c(-0.02, 0.02) * diff(yrange),
    label = c("Ext.", "Imm.")
  ) + ggplot2::scale_color_manual(
    values = patchBackgrounds, aesthetics = c("color", "fill")
  )
}

anyAbundanceSoFar <- FALSE
largestTotalEffect <- 0

# Create a common layout.
timestepLayout <- which.max(rowSums(result$Abundance[timesInUse, -1] > 0))
graf4Layout <- createBaseGraph(timestepLayout, environs)
graf4Layout <- addIntrinsicToGraph(graf4Layout, timestepLayout, environs)
graf4Layout <- thresholdGraphEdges(graf4Layout, threshold)
lay <- computeGraphLayout(graf4Layout)

# This might cause problems with the diet???
# temp <- apply(result$Abundance[, -1] > result$Parameters$EliminationThreshold,
#               2, which.max)
# temp <- temp[names(temp) %in% lay$name]
# temp <- temp[order(names(temp))]
# lay <- lay %>% dplyr::arrange(name) %>% dplyr::mutate(x = temp)
#
#
# # GIF
# animation::saveGIF(
#   movie.name = "test.gif",
# Video
animation::saveVideo(
  # video.name = "test.mp4",
  video.name = paste0("test2_", tag, ".mp4"),
  ani.height = 1024 * 2, ani.width = 1280 * 2, interval = 0.1,
  expr = {
    for (timestep in seq(from = 1, to = nrow(environs[[1]]$Abundance),
                         # length.out = nrow(result$Events)/1000
                         length.out = nrow(result$Events)/1
    )) {
      # ani.options(ani.height = 1280 * 2, ani.width = 1024 * 2, interval = 1)

      timestep <- round(timestep)

      if (!anyAbundanceSoFar) {
        if (!any(result$Abundance[timesInUse,][timestep, -1] > 0)) {
          next
        } else {
          anyAbundanceSoFar <- TRUE
        }
      }

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
        max(abs(unlist(lapply(graf, function(g) g %N>% pull(fill)))),
            na.rm = TRUE)
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


          # jacWith1 <- jaccard %>% dplyr::filter(
          #   Env2 == i
          # )
          # jacWith1 <- jacWith1[
          #   which.min(abs(jacWith1$Time - times[timestep])),
          #   ]$Jaccard #%>% signif(digits = 4)
          #
          #
          # bcWith1 <- braycurtis %>% dplyr::filter(
          #   Env2 == i
          # )
          # bcWith1 <- bcWith1[
          #   which.min(abs(bcWith1$Time - times[timestep])),
          #   ]$BrayCurtis #%>% signif(digits = 4)
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
        ) + ggplot2::geom_rect(
          data = data.frame(1),
          # type = "rect", # older version or ggplot, so geom_rect.
          xmin = min(lay$x),
          xmax = max(lay$x),
          ymin = min(lay$y),
          ymax = max(lay$y),
          alpha = 0.15,
          fill = if(i == 1) {
            patchBackgrounds[1]
          } else if (i == length(graf)) {
            patchBackgrounds[2]
          }
        ) + ggraph::geom_node_point(
          mapping = aes(
            #color = Type,
            alpha = interaction(Present, On),
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
        ) + ggplot2::geom_hline(
          yintercept = -1, linetype = "dashed", color = "black"
        ) + scale_shape_manual(
          values = c(22, 23), name = "b:"
        ) + scale_alpha_manual(
          values = c("FALSE.FALSE" = 0, "FALSE.TRUE" = 0.25, "TRUE.TRUE" = 1),
          guide = "none"
          # ) + scale_color_manual(
          #   name = "Node",
          #   values = c("Basal" = "darkgreen", "Consumer" = "yellow3")
        ) + scale_fill_gradient2(
          low = "#FF00FF", mid = "#FFFFFF", high = "#00FFFF",
          limits = c(-6, 6), breaks = c(-4,0,4),
          transform = scales::transform_pseudo_log(sigma = 1/10),
          name = "Per. Cap. Rate"
        ) + scale_edge_linetype_manual(
          name = "",
          values = c("Positive" = "dashed", "Negative" = "longdash"),
          guide = "none" # Legend entry not particularly clear.
        ) + scale_edge_color_manual(
          name = "Edge",
          values = edgecolors, drop = FALSE
        ) + scale_edge_alpha(
          name = "Contribution %",
          transform = "sqrt", breaks = c(1e-5, 1e-3, 1e-1),
          limits = c(threshold, 1)
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
                  c("R:"),#, "Jac:", "BC:"),
                  unlist(lapply(c(corWith1),#, jacWith1, bcWith1),
                                function(x) {
                                  if (is.na(x)) {
                                    " NA        "
                                  } else if (x == 1) {
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

      # https://archive.schochastics.net/post/ggraph-tricks-for-common-problems/
      pointsAsBarPlots <- lapply(graf, function(g) {
        # GAINS VERSION
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
        # DIET VERSION:
        df <- g %E>% igraph::as_data_frame(
        ) %>% dplyr::filter(
          Type == "Consumption"
        ) %>% dplyr::mutate(
          from = as.numeric(from)
        ) %>% dplyr::left_join(
          Pool %>% dplyr::select(
            ID, Type
          ) %>% dplyr::rename(
            Prey = Type
          ),
          by = c("from" = "ID")
        ) %>% dplyr::group_by(
          to, Prey
        ) %>% dplyr::summarise(
          Diet = sum(PercentContribution),
          .groups = "drop"
        ) %>% dplyr::rename(
          Node = to
        )

        df %>% dplyr::group_by(Node) %>% dplyr::group_map(
          .f = function(.x, .y) {
            grob <- ggplot2::ggplotGrob(
              ggplot2::ggplot(
                .x
              ) + ggplot2::geom_col(
                ggplot2::aes(x = 1, fill = Prey, y = Diet)
              ) + ggplot2::scale_fill_manual(
                values = c("Basal" = "darkgreen", "Consumer" = "goldenrod")
              ) + ggplot2::coord_cartesian(
                expand = FALSE
              ) + ggplot2::theme(
                legend.position = "none",
                panel.background =
                  ggplot2::element_rect(fill = "white", color = NA),
                line = ggplot2::element_blank(),
                text = ggplot2::element_blank()
              )
            )
            coordinates <- grob$layout[grob$layout$name == "panel", ]
            list(
              node = .y$Node,
              grobs = grob[coordinates$t:coordinates$b,
                           coordinates$l:coordinates$r]
            )
          }
        )
      }
      )

      annotationspace <- 0.14

      barplotsAsAnnotations <- lapply(pointsAsBarPlots, function(g) {
        lapply(g, function(bp) {
          node <- bp$node
          grob <- bp$grobs
          xy <- lay %>% dplyr::filter(name == node) %>% dplyr::select(x, y)
          # print(xy)
          ggplot2::annotation_custom(
            grob,
            xmin = xy$x - annotationspace * 2,
            xmax = xy$x + annotationspace * 2,
            ymin = xy$y - annotationspace / 10,
            ymax = xy$y + annotationspace / 10
          )
        })
      })

      combinedplots <- lapply(seq_along(graf), function(i) {
        Reduce("+", barplotsAsAnnotations[[i]], plots[[i]])
      })

      dashboard <- lapply(
        list(
          # Alpha Richness
          ggplot2::ggplot(
            Diversity$Diversity %>% dplyr::filter(
              Metric == "Alpha Hill:0",
              is.na(Subset)
            ) %>% mutate(
              Environment1 = as.factor(Environment1)
            ),
            ggplot2::aes(x = Time, y = Value,
                         color = Environment1, linetype = Environment1)
          ) + ggplot2::geom_line(
          ) + ggplot2::ylab(
            "Richness"
          ),
          ggplot2::ggplot(
            Diversity$Diversity %>% dplyr::filter(
              Metric == "TimeJaccard",
              is.na(Subset)
            ) %>% mutate(
              Environment1 = as.factor(Environment1)
            ),
            ggplot2::aes(x = Time, y = Value,
                         color = Environment1, linetype = Environment1)
          ) + ggplot2::geom_line(
          ) + ggplot2::ylab(
            "Temporal Jaccard"
          ),
          # # Gamma Richness
          # ggplot2::ggplot(
          #   Diversity$Diversities$gamma %>% dplyr::filter(
          #     Aggregation == "Gamma"
          #   ),
          #   ggplot2::aes(x = Time, y = Richness)
          # ) + ggplot2::geom_line(
          # ),
          # # Beta Jaccard
          # ggplot2::ggplot(
          #   Diversity$Diversities$beta %>% dplyr::bind_rows() ,
          #   ggplot2::aes(x = Time, y = Jaccard,
          #                group = interaction(Env1, Env2))
          # ) + ggplot2::geom_line(
          # ),
          # Alpha exp(Evenness)
          ggplot2::ggplot(
            Diversity$Diversity %>% dplyr::filter(
              Metric == "Alpha Hill:1",
              is.na(Subset)
            ) %>% mutate(
              Environment1 = as.factor(Environment1)
            ),
            ggplot2::aes(x = Time, y = Value,
                         color = Environment1, linetype = Environment1)
          ) + ggplot2::geom_line(
          ) + ggplot2::ylab("Exp(Evenness)"),
          ggplot2::ggplot(
            Diversity$Diversity %>% dplyr::filter(
              Metric == "TimeBrayCurtis",
              is.na(Subset)
            ) %>% mutate(
              Environment1 = as.factor(Environment1)
            ),
            ggplot2::aes(x = Time, y = Value,
                         color = Environment1, linetype = Environment1)
          ) + ggplot2::geom_line(
          ) + ggplot2::ylab(
            "Temporal Bray Curtis"
          )#,
          # # Gamma exp(Evenness)
          # ggplot2::ggplot(
          #   DivAbund$DivAbund$gamma,
          #   ggplot2::aes(x = Time, y = `1, All`)
          # ) + ggplot2::geom_line(
          # ) + ggplot2::ylab("Exp(Evenness)"),
          # # Beta Bray-Curtis
          # ggplot2::ggplot(
          #   DivAbund$DivAbund$beta ,
          #   ggplot2::aes(x = Time, y = BrayCurtis,
          #                group = interaction(Env1, Env2))
          # ) + ggplot2::geom_line(
          # )
        ),
        plotEventRug
      )

      dashboard <- lapply(dashboard, function(gplot) {
        gplot + ggplot2::geom_vline(
          color = "black", xintercept = times[timestep], linetype = "dashed"
        )
      })

      tempplot <- ggpubr::ggarrange(
        plotlist = combinedplots, # plots
        # ncol = 2, nrow = 1,
        common.legend = TRUE,
        legend = "bottom"
      )

      tempdashboard <- ggpubr::ggarrange(
        plotlist = dashboard,
        # ncol = 3, nrow = 2,
        common.legend = TRUE,
        legend = "bottom"
      )

      tempfull <- ggpubr::ggarrange(
        plotlist = list(tempplot, tempdashboard),
        ncol = 1, nrow = 2, heights = c(3, 1)
      )

      plot(tempfull)
      print(largestTotalEffect)
    }
  }
)
