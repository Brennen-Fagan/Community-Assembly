# From FoodWebPostMortem.R, 2023-02-24
toCheddar <- function(EVList, name = " ") {# Edges Vertices List
  if (all(is.na(EVList$Edges)) && all(is.na(EVList$Vertices))) {
    return(NA)
  }

  links <- EVList$Edges

  # cheddar does not like "cannibalism".
  links <- links[
    links$to != links$from,
    ]

  # "[C]olumns called ‘resource’ and ‘consumer’ must be given."
  links <- dplyr::bind_rows(
    links |> dplyr::filter(effectSign == 1) |> dplyr::rename(
      resource = from, consumer = to),
    links |> dplyr::filter(effectSign == -1) |> dplyr::rename(
      resource = to, consumer = from)
  ) |> dplyr::select(-Type) # Cheddar confuses node Type and edge Type.

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
