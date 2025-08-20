extractAlphas <- function(Diversity, Attributes) {
  dplyr::bind_rows(
    lapply(
      seq_along(Diversity$Diversities),
      function(i, d, a) assignAttributes(i, d[[i]]$alpha, a),
      d = Diversity$Diversities,
      a = Attributes
    )
  ) %>% dplyr::select(-Species) %>% dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) %>% dplyr::group_by(
    Time, Environment, Set, Number, History, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    Richness = floor(median(Richness)),
    Richness_Basal = floor(median(Richness_Basal)),
    Richness_Consumer = floor(median(Richness_Consumer)),
    .groups = "drop"
  )
}