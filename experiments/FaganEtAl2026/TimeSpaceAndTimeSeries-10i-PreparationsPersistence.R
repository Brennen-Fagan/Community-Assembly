# Comments: ###################################################################
# Note: there are 187 false dispersal events and 16 double dynamic loss events.
# This is out 9,721,730 pairs of events.

# Of the double dynamic losses, 4 of them had time gaps > 1e-11:
#   0.14, 0.15, 0.28, and 143.44
# so we set the value function to be the maximum of these times.
# Details recorded:
# Species Environment SpeciesType Size ReproductionRate Speed Affinity
# AffinityBins PoolPatch PoolPatchSeed Interactions InteractionsSeed Events
# EventsSeed InitialConditions InitialConditionsSeed DispersalParam
# NicheDistance SpeciesAffinity SpeciesAffinitySeed PatchAffinity
# PatchAffinitySeed InterventionPatchType InterventionPatchSeed
# InterventionTimeType InterventionTimeSeed InterventionDispersal
# InterventionNicheDistance InNumber Arrival Extinct Dynamic Loss
# 117 1 Consumer 0.20918166 -0.09185240 1 0 0 142486 10 4929 10 28 10 1 10 NA
# 3 7 30 4 49 113 1 1 1 p p 2 val 16721.4202298989 NA diff 0.139955025184463
# 50 1 Basal 0.06903233 0.20094599 1 0 0 142486 1 4929 1 28 1 1 1 NA 5 1 1 2 2
# 111 1 1 1 p p 2 val 16269.1906969154 NA diff 0.15305802912917
# 117 1 Consumer 0.20918166 -0.09185240 1 0 0 142486 10 4929 10 28 10 1 10 NA
# 3 7 30 4 49 112 1 1 1 p p 2 val 16721.4202298989 NA diff 0.27991005036165
# 117 1 Consumer 0.20918166 -0.09185240 1 0 0 142486 10 4929 10 28 10 1 10 NA
# 3 7 30 4 49 111 1 1 1 p p 2 val 16721.4202298989 NA diff 143.435277663946

# Of the false dispersal events, this number is inflated as the pre-intervention
# figures are copied from the non-intervention cases. The length of time these
# species remain in the simulation nonetheless vary, so we must acknowledge that
# they do have an effect. On the other hand, they are all in the strong land-use
# interaction case meaning they are not our primary dataset and previous
# analysis of this problem indicated they were not detectable as having
# statistically differed in terms of arrival times. We choose to leave them in
# the data as a result, but exclude them from analyses.

# Preparation: ################################################################
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")

# If not load, then run.

loadPers <- TRUE
fileNamePers <- "TSTS_Pers_10a1.RData"

if (!exists("Pers")) {
  if (loadPers && file.exists(fileNamePers)) {
      load(fileNamePers)
    } else {
      # ~ 10 Minutes.
      # Persistences: #########################################################
      Pers <- ColExt |> tidytable::rename(
        DispersalParam = Dispersal
      ) |> tidytable::filter(
        EventType != "Present", Success # False Arrivals might mess this up.
      ) |> tidytable::group_by(
        Species, Environment, SpeciesType, Size, ReproductionRate, Speed,
        Affinity, AffinityBins,
        PoolPatch:InterventionNicheDistance
      ) |> tidytable::mutate(
        InNumber = ifelse(EventType == "Arrival" |
                            EventType == "Dispersal", 1, 0),
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
          PoolPatch:InterventionNicheDistance,
          # And the true ID variable!
          InNumber
        ),
        values_fill = NA,
        values_fn = max # see comments above.
      ) |> tidytable::mutate(
        In = ifelse(is.na(Dispersal), Arrival, Dispersal),
        Out = ifelse(is.na(Extinct),
                     ifelse(is.na(`Dynamic Loss`),
                            EndOfSimulation,
                            `Dynamic Loss`),
                     Extinct),
        InType = externalNames[ifelse(is.na(Dispersal),
                                      "Arrival", "Dispersal")],
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
