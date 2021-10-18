# For playing with Edge-Vertex Lists, especially plotting.
# See LM1996-NumPoolCom-FoodWebs-2021-07.Rmd

toCheddar <- function(EVList, name = " ") {# Edges Vertices List
  if (is.na(EVList$Edges) && is.na(EVList$Vertices)) {
    return(NA)
  }

  links <- EVList$Edges

  # cheddar does not like "cannibalism".
  links <- links[
    links$to != links$from,
  ]

  # "[C]olumns called ‘resource’ and ‘consumer’ must be given."
  links <- dplyr::bind_rows(
    links %>% dplyr::filter(effectSign == 1) %>% dplyr::rename(
      resource = from, consumer = to),
    links %>% dplyr::filter(effectSign == -1) %>% dplyr::rename(
      resource = to, consumer = from),
  ) %>% dplyr::select(-Type) # Cheddar confuses node Type and edge Type.

  cheddar::Community(
    nodes = EVList$Vertices,
    properties = list(
      title = name,
      M.units = "masses",
      N.units = "abund"
    ),
    trophic.links = links
  )
}

toIGraph <- function(EVList, sign = 0) {
  igraph::graph_from_data_frame(
    d = if(sign == 0) {
      EVList$Edges
    } else {
      EVList$Edges[EVList$Edges$effectSign == sign, ]
    },
    directed = TRUE,
    vertices = EVList$Vertices
  )
}


toPostMortem <- function(EVList,
                         threshold = 0, # sets to minimal size edges below
                         nodeSize = c("None", "Abundance", "Size"),
                         edgeScale = 10,
                         reducedTrophic = TRUE) {
  if (tolower(threshold) == "adaptive") {
    threshold = EVList$Edges %>% group_by(
      to, effectSign
    ) %>% summarise(
      max = max(effectNormalised), .groups = "drop"
    ) %>% ungroup %>% pull(max) %>% min
  }

  theGc <- toCheddar(EVList, name = "Trophic Levels")
  theGi <- toIGraph(EVList)

  theGiGain <- toIGraph(EVList, sign = 1)
  theGiLoss <- toIGraph(EVList, sign = -1)

  theLayout <- igraph::layout.circle(theGi)

  theSize <- match.arg(nodeSize, c("Abundance", "Size", "None"))
  if (theSize == "Abundance")
    theVs <- sqrt(igraph::vertex_attr(theGi)$N) * 10
  else if (theSize == "Size") {
    theVs <- igraph::vertex_attr(theGi)$M
    theVs <- sqrt(theVs / min(theVs)) * 10
  } else if (theSize == "None") {
    theVs <- 15
  }

  theColors <- ifelse(
    igraph::vertex_attr(theGi)$Type == "Basal", "skyblue", "red"
  )
  if ("Core" %in% names(igraph::vertex_attr(theGi))) {
    theShapes <- ifelse(igraph::vertex_attr(theGi)$Core,
                        0,
                        1)
  } else {
    theShapes <- 1
    # Note Igraph uses "circle" then "rectangle",
    # but R and cheddar use "rectangle" then "circle", so we will use a !.
  }

  theBoth <- igraph::edge_attr(theGi)$effectNormalised
  theGain <- igraph::edge_attr(theGiGain)$effectNormalised
  theLoss <- igraph::edge_attr(theGiLoss)$effectNormalised

  theBoth[theBoth < threshold] <- 0
  theGain[theGain < threshold] <- 0
  theLoss[theLoss < threshold] <- 0

  # Inform the graphs of which edges are not needed.
  theGi <- igraph::delete_edges(theGi, which(theBoth == 0))
  theGiGain <- igraph::delete_edges(theGiGain, which(theGain == 0))
  theGiLoss <- igraph::delete_edges(theGiLoss, which(theLoss == 0))

  # Remove the same entries so that lengths match.
  theGain <- theGain[theGain > 0]
  theLoss <- theLoss[theLoss > 0]

  theGain <- theGain * edgeScale
  theLoss <- theLoss * edgeScale

  parold <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2), # Two Rows, Two Columns
      mar = c(0, 1.5, 1, 0), # Margins, bottom, left, top, right
      oma = c(0.1, 0.1, 0.1, 0.1) # Outer margins.
  )

  cheddar::PlotWebByLevel(
    theGc,
    show.level.lines = TRUE,
    # Had been using LongWeighted, but that seems to give the upside down T.
    # Flow based seems to be more what we are expecting, given the usage of
    # thresholding and what that shows. The flows here are expected to be
    # flows of energy through the food web.
    level = cheddar::FlowBasedTrophicLevel(theGc, weight.by = "effectNormalised"),
    col = theColors,
    pch = theShapes
  )

  if (!reducedTrophic) {
    plot(
      theGi,
      layout = theLayout,
      vertex.size = theVs,
      edge.width = 1,
      edge.arrow.size = 0.3,
      edge.arrow.width = 1,
      vertex.color = theColors,
      vertex.shape = igraph::shapes()[as.numeric(!theShapes) + 1],
      edge.lty = 2,
      edge.color = "grey",
      edge.arrow.mode = ">",
      main = "Consumption"
    )
  } else {
    EVListRed <- EVList
    EVListRed$Edges <- EVListRed$Edges %>% dplyr::filter(
      effectNormalised >= threshold
    )
    theGc2 <- toCheddar(EVListRed, name = "Strongest Trophic Levels")
    cheddar::PlotWebByLevel(
      theGc2,
      show.level.lines = TRUE,
      level = cheddar::FlowBasedTrophicLevel(theGc2, weight.by = "effectNormalised"),
      col = theColors,
      pch = theShapes
    )
  }

  plot(
    theGiGain,
    layout = theLayout,
    vertex.size = theVs,
    edge.width = theGain,
    edge.arrow.size = 0.3,
    edge.arrow.width = 1,
    vertex.color = theColors,
    vertex.shape = igraph::shapes()[as.numeric(!theShapes) + 1],
    edge.lty = 2,
    edge.color = "blue",
    edge.arrow.mode = ">",
    main = "Consumer's Gains"
  )

  plot(
    theGiLoss,
    layout = theLayout,
    vertex.size = theVs,
    edge.width = theLoss,
    edge.arrow.size = 0.3,
    edge.arrow.width = 2,
    vertex.color = theColors,
    vertex.shape = igraph::shapes()[as.numeric(!theShapes) + 1],
    edge.lty = 3,
    edge.color = "darkred",
    edge.arrow.mode = "<",
    main = "Resource's Losses"
  )

  par(parold)

  EVList$Edges %>% dplyr::ungroup() %>% dplyr::filter(
    effectNormalised >= threshold
  ) %>% dplyr::select(
    -effectSign
  ) %>% dplyr::arrange(
    to, -effectNormalised
  )
}
