# As in main package, but with the dependency changed. (RMTRCode2::toCheddar)
CalculateTrophicStructure <- function(
  Pool,
  ReproductionRate = NULL, # Probably easiest to make many calculators
  NumEnvironments,                  # rather than functionalise.
  InteractionMatrices,
  EliminationThreshold
) {
  if (is.null(ReproductionRate)) {
    EffectiveReproductionRate <- function(y) Pool$ReproductionRate
  } else if (!is.function(ReproductionRate)) {
    EffectiveReproductionRate <- function(y) ReproductionRate
  } else { # !NULL & is.function
    EffectiveReproductionRate <- ReproductionRate
  }

  # Borrowing from LM1996-NumPoolCom-FoodWebs-2021-07.Rmd
  nrowPool <- nrow(Pool)
  `%>%` <- magrittr::`%>%`

  # This function should be appliable row-wise to the results.
  # One does need to remove the time column, as usual.
  function(y) {
    # Clean up anything not present.
    y <- ifelse(y <= EliminationThreshold, 0, y)

    # Break y up into its environments.
    EnvsY <- lapply(1:NumEnvironments, function(i, ys) {
      ys[(i - 1) * nrowPool + 1:nrowPool]
    }, ys = y)

    EnvsEdgeVertexLists <- lapply(
      seq_along(EnvsY),
      function(i, envY, mats) {
        # "Red"uced "Com"munity; who's present.
        redCom <- which(envY[[i]] > 0)

        if (length(redCom) == 0) {
          return(list(
            Edges = NA,
            Vertices = NA
          ))
        }

        redPool <- Pool[redCom, ]
        redMat <- matrix(mats[[i]][redCom, redCom],
                         nrow = length(redCom),
                         ncol = length(redCom))

        colnames(redMat) <- paste0('s',as.character(redCom))
        rownames(redMat) <- colnames(redMat)

        names(redPool)[1] <- "node"
        redPool$node <- colnames(redMat)

        Graph <- igraph::graph_from_adjacency_matrix(
          redMat, weighted = TRUE
        )

        Graph <- igraph::set.vertex.attribute(
          Graph, "name", value = colnames(redMat)
        )

        redPool$N <- envY[[i]][redCom]

        # For later analysis, take the matrix diagonal.

        redPool$Intraspecific <- diag(redMat)

        GraphAsDataFrame <- igraph::as_data_frame(Graph)

        # Add in abundances for calculating abundance * (gain or loss)
        GraphAsDataFrame <- dplyr::left_join(
          GraphAsDataFrame,
          dplyr::select(redPool, node, N),
          by = c("to" = "node")
        )

        if (("weight" %in% colnames(GraphAsDataFrame))) {
          # We're in a case where there are edges.
          # In the opposite case, we cannot do this part of the calculation.
          # Split data frame.
          ResCon <- GraphAsDataFrame[GraphAsDataFrame$weight > 0,]
          ConRes <- GraphAsDataFrame[GraphAsDataFrame$weight < 0,]

          # Reorder and rename variables.
          ResCon <- dplyr::select(ResCon,
                                  to, from, # resource = to, consumer = from,
                                  effectPerUnit = weight, resourceAbund = N)
          ConRes <- dplyr::select(ConRes,
                                  to, from, # resource = from, consumer = to,
                                  effectPerUnit = weight, consumerAbund = N)
          ResCon <- dplyr::mutate(dplyr::group_by(ResCon, from),
                                  effectActual = effectPerUnit * resourceAbund,
                                  Type = "Exploit+")
          ConRes <- dplyr::mutate(dplyr::group_by(ConRes, from),
                                  effectActual = effectPerUnit * consumerAbund,
                                  Type = ifelse(from == to,
                                                "SelfReg-",
                                                "Exploit-"))
        }

        # the withs here should now be mostly deprecated. N, node inside.
        reprate <- EffectiveReproductionRate(y)[redCom]
        IntriG <- with(redPool, data.frame(
          from = node, #resource = node,
          to = node, #consumer = node,
          effectPerUnit = ifelse(reprate > 0,
                                 reprate, 0),
          effectActual = ifelse(reprate > 0,
                                N * reprate, 0),
          Type = "Intrisc+",
          stringsAsFactors = FALSE))
        IntriL <- with(redPool, data.frame(
          from = node, #resource = node,
          to = node, #consumer = node,
          effectPerUnit = ifelse(reprate < 0,
                                 reprate, 0),
          effectActual = ifelse(reprate < 0,
                                N * reprate, 0),
          Type = "Intrisc-",
          stringsAsFactors = FALSE))

        if (exists("ResCon")) {
          ResCon <- dplyr::select(ResCon, -resourceAbund)
        } else {
          ResCon <- data.frame()
        }
        if (exists("ConRes")) {
          ConRes <- dplyr::select(ConRes, -consumerAbund)
        } else {
          ConRes <- data.frame()
        }

        EdgeDataFrame <- dplyr::bind_rows(
          ResCon, ConRes,
          IntriG, IntriL
        )

        EdgeDataFrame <- EdgeDataFrame %>% dplyr::rename(
          # Empirically speaking, to and from appear reversed.
          # A consumer (from) should have a negative effect on resource (to),
          # but the organisation so far marks it as positive. We fix this.
          tempname = to,
          to = from
        ) %>% dplyr::rename(
          from = tempname
        ) %>% dplyr::filter(
          # Remove placeholder entries
          effectPerUnit != 0
        ) %>% dplyr::mutate(
          # Useful to keep effects separate
          effectSign = sign(effectPerUnit)
        ) %>% dplyr::group_by(
          to, effectSign
        ) %>% dplyr::mutate(
          # Perform the post mortem of the most influential from's
          effectEfficiency = effectPerUnit / sum(effectPerUnit),
          effectNormalised = effectActual / sum(effectActual)
        ) %>% dplyr::arrange(to)

        list(
          Edges = EdgeDataFrame,
          Vertices = redPool
        )

      },
      envY = EnvsY,
      mats = InteractionMatrices$Mats)

    EnvsCheddar <- lapply(EnvsEdgeVertexLists, toCheddar)

    EnvsTrophic <- lapply(EnvsCheddar,
                          function(x, weight.by) {
                            if (all(!is.na(x)))
                              cheddar::TrophicLevels(x, weight.by = weight.by)
                            else NA
                          },
                          weight.by = "effectNormalised")

    # In principle, I think these are the two return values.
    # In practice, it seems more useful to return the EdgeVertexLists and
    # the Trophic Levels, given the importance of intraspecific interactions.
    # These are what Cheddar does not capture.

    return(list(
      EdgeVertexLists = EnvsEdgeVertexLists,
      TrophicLevels = EnvsTrophic
    ))
  }

}
