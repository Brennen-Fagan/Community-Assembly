
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesRichness")) {
  load("diversitiesFlattened10a1_subsetRichness.RData")

  diversitiesRichness <- diversitiesRichness |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}

# if (!exists("InterventionTimes")) {
#   load("TSTS_Interventions_10a1.RData")
# }



