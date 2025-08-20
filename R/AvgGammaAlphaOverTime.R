AvgGammaAlphaOverTime <- function(DAG) {
  DAG |> dplyr::group_by(
    Time, Set, Number, Pool, Noise, Neutral, Space
  ) |> dplyr::summarise(
    `Richness, Alpha` = mean(`Richness, Alpha`),
    `Richness_Basal, Alpha` = mean(`Richness_Basal, Alpha`),
    `Richness_Consumer, Alpha` = mean(`Richness_Consumer, Alpha`),
    `Richness, Gamma` = mean(`Richness, Gamma`),
    `Richness_Basal, Gamma` = mean(`Richness_Basal, Gamma`),
    `Richness_Consumer, Gamma` = mean(`Richness_Consumer, Gamma`),
    .groups = "drop"
  ) |> tidyr::pivot_longer(
    cols = c("Richness, Alpha", "Richness, Gamma",
             "Richness_Basal, Alpha", "Richness_Consumer, Alpha",
             "Richness_Basal, Gamma", "Richness_Consumer, Gamma"),
    names_to = "Measurement", values_to = "Value"
  ) |> ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = Time, y= Value,
      group = Measurement, color = Measurement
    )
  ) + ggplot2::geom_line(
  )
}