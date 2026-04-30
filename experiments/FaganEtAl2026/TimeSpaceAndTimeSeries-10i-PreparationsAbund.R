
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesAbund")) {
  load(paste0("diversitiesFlattened10a1_", dateSwitch, "_Abund.RData"))
  # load("diversitiesFlattened10a1_2025-09-04_Abund.RData")

  diversitiesAbund <- diversitiesAbund |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}

# if (!exists("InterventionTimes")) {
#   load("TSTS_Interventions_10a1.RData")
# }

