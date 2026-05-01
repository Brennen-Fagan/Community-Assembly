# Depends on "assignAttributes".
extractBetas <- function(Diversity, Attributes) {
  dplyr::bind_rows(
    lapply(
      seq_along(Diversity$Diversities),
      function(i, d, a) assignAttributes(i, d[[i]]$beta, a),
      d = Diversity$Diversities,
      a = Attributes
    )
  ) |> dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) |> dplyr::group_by(
    Time, Env1, Env2, Set, Number, History, Pool, Noise, Neutral, Space
  ) |> dplyr::summarise(
    Jaccard = median(Jaccard),
    .groups = "drop"
  )
}