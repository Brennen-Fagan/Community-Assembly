unifyAffinityBins <- function(., n = 5) {
  tidytable::separate_wider_delim(
    .,
    col = "AffinityBins", names = c("Left", "Right"), delim = ","
  ) |> tidytable::mutate(
    Left =
      round(as.numeric(gsub(pattern = "^[(]", replacement = "", x = Left))*n)/n,
    Right =
      round(as.numeric(gsub(pattern = "\\]$", replacement = "", x = Right))*n)/n,
    AffinityBins = ifelse(
      is.na(Right), as.character(Left),
      paste0("(", Left, ", ", Right, "]")
    )
  )
}