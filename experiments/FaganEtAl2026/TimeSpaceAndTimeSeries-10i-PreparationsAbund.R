
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesAbund")) {
  load(paste0("diversitiesFlattened10a1_", dateSwitch, "_Abund.RData"))

  diversitiesAbund <- diversitiesAbund |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}

