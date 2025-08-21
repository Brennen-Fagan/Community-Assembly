changeAffinityLevels <- function(df) {
  df |> tidytable::mutate(
    SpeciesAffinity = tidytable::case_when(
      SpeciesAffinity == "rep_0" ~ "100% 0",
      SpeciesAffinity == "evensplit_01" ~ "50% 0, 50% 1",
      SpeciesAffinity == "runif" ~ "Uniform(0, 1)",
      TRUE ~ SpeciesAffinity
    ),
    SpeciesAffinity = factor(SpeciesAffinity, levels = c(
      "100% 0", "50% 0, 50% 1", "Uniform(0, 1)"
    ), ordered = TRUE)
  )
}