
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesRichness")) {
  paste0("diversitiesFlattened10a1_", dateSwitch, "_Richness.RData"))
  # load("diversitiesFlattened10a1_2025-09-04_Richness.RData")

  diversitiesRichness <- diversitiesRichness |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}

# if (!exists("InterventionTimes")) {
#   load("TSTS_Interventions_10a1.RData")
# }



