# From FoodWebPostMortem.R, 2023-02-24
toIGraph <- function(EVList, sign = 0) {
  `%>%` <- magrittr::`%>%`

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