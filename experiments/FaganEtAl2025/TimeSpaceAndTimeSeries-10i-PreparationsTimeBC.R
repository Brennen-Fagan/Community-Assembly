
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

if (!exists("diversitiesTimeBC")) {
  load("diversitiesFlattened10a1_subsetTimeBC.RData")
  
  diversitiesTimeBC <- diversitiesTimeBC |> changePreferencesLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::left_join(
    endTimes |> dplyr::select(-Times)
  )
}

# if (!exists("InterventionTimes")) {
#   load("TSTS_Interventions_10a1.RData")
# }



