
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesRichness")) {
  load(paste0("diversitiesFlattened10a1_", dateSwitch, "_Richness.RData"))

  diversitiesRichness <- diversitiesRichness |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}



