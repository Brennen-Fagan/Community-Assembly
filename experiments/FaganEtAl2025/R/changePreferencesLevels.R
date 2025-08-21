changePreferencesLevels <- function(df) {
  df |> tidytable::mutate(
    SpeciesPreferences = tidytable::case_when(
      SpeciesPreferences == "rep_0" ~ "100% 0",
      SpeciesPreferences == "evensplit_01" ~ "50% 0, 50% 1",
      SpeciesPreferences == "runif" ~ "Uniform(0, 1)",
      TRUE ~ SpeciesPreferences
    ),
    SpeciesPreferences = factor(SpeciesPreferences, levels = c(
      "100% 0", "50% 0, 50% 1", "Uniform(0, 1)"
    ), ordered = TRUE)
  )
}
