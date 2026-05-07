
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesTimeBC")) {
  load(paste0("diversitiesFlattened10a1_", dateSwitch, "_TimeBC.RData"))
  # load("diversitiesFlattened10a1_2025-09-04_TimeBC.RData")

  diversitiesTimeBC <- diversitiesTimeBC |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}

# if (!exists("InterventionTimes")) {
#   load("TSTS_Interventions_10a1.RData")
# }



