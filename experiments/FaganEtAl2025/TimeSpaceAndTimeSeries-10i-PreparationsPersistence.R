source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

# If not load, then run.

loadPers <- TRUE
fileNamePers <- "TSTS_Pers_10a1.RData"

# ~ 10 Minutes.
if (!exists("Pers")) {
  if (loadPers) {
      load(fileNamePers)
    } else {
      # Persistences: ################################################################
      Pers <- ColExt |> tidytable::rename(
        DispersalParam = Dispersal
      ) |> tidytable::filter(
        EventType != "Present", Success # False Arrivals might mess this up.
      ) |> tidytable::group_by(
        Species, Environment, SpeciesType, Size, ReproductionRate, Speed, 
        Affinity, AffinityBins, 
        PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, 
        EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
        NicheDistance, 
        SpeciesAffinity, SpeciesAffinitySeed, PatchAffinity, PatchAffinitySeed,
        InterventionPatchType, InterventionPatchSeed,
        InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
        InterventionNicheDistance
      ) |> tidytable::mutate(
        InNumber = ifelse(EventType == "Arrival" | EventType == "Dispersal", 1, 0),
        InNumber = cumsum(InNumber)
      ) |> tidytable::ungroup(
      ) |> tidytable::pivot_wider(
        values_from = "Times",
        names_from = EventType,
        id_cols = c(
          # THERE HAS TO BE SOMETHING NICER HERE???
          # Species/Environment/Simulation Identifiers
          Species, Environment, SpeciesType, Size, ReproductionRate, Speed, 
          Affinity, AffinityBins, 
          PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, 
          EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
          NicheDistance, 
          SpeciesAffinity, SpeciesAffinitySeed, PatchAffinity, PatchAffinitySeed,
          InterventionPatchType, InterventionPatchSeed,
          InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
          InterventionNicheDistance,
          # And the true ID variable!
          InNumber
        ),
        values_fill = NA
      ) |> tidytable::mutate(
        In = ifelse(is.na(Dispersal), Arrival, Dispersal),
        Out = ifelse(is.na(Extinct),
                     ifelse(is.na(`Dynamic Loss`),
                            EndOfSimulation,
                            `Dynamic Loss`),
                     Extinct),
        InType = externalNames[ifelse(is.na(Dispersal), "Arrival", "Dispersal")],
        OutType = externalNames[ifelse(is.na(Extinct),
                                       ifelse(is.na(`Dynamic Loss`),
                                              "EndOfSimulation",
                                              "Dynamic Loss"),
                                       "Extinct")],
        Persistence = Out - In,
        # Enhance Readability:
        SpeciesPreferences =
          speciesAffinityDictionaryOrigin$SpeciesAffinities[
            as.numeric(SpeciesAffinity)
            ]
      ) |> changePreferencesLevels(
      ) |> tidytable::select(
        -Dispersal, -Arrival, -Extinct, -`Dynamic Loss`, -EndOfSimulation
      ) |> tidytable::left_join(
        interventionStrings,
        by = c("PatchAffinity", "PoolPatch", "InterventionPatchType"),
        multiple = "all"
      ) |> changeInterventionLevels(
      ) |> tidytable::left_join(
        endTimes
      ) |> unifyAffinityBins()
      
      save(Pers, file = fileNamePers)
    }
}
