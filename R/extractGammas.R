extractGammas <- function(Diversity, Attributes) {
  dplyr::bind_rows(lapply(
    seq_along(Diversity$Diversities),
    function(i, d, a) assignAttributes(i, d[[i]]$gamma, a),
    d = Diversity$Diversities,
    a = Attributes)) %>% dplyr::filter(
      Aggregation == "Gamma"
    ) %>% dplyr::mutate(
      Time = floor(Time * time_grouping_size)/time_grouping_size
    ) %>% dplyr::group_by(
      Time, Set, Number, History, Pool, Noise, Neutral, Space
    ) %>% summarise(
      Richness = floor(median(Richness)),
      Richness_Basal = floor(median(Basals)),
      Richness_Consumer = floor(median(Consumers)),
      .groups = "drop"
    )
}