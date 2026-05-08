
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesTimeBC")) {
  load(paste0("diversitiesFlattened10a1_", dateSwitch, "_TimeBC.RData"))

  diversitiesTimeBC <- diversitiesTimeBC |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}



