# Statistics used in the manuscript.
# library(fitdistrplus)

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsTimeBC.R")

supplementStatistics <- list()

# Extirpations and Turnover: ##################################################
# "These weak interactions are reflected in the way the system changes through
# time; losses from the community are dominated by stochastic extirpations (no
# adaptation: 93.6%) and turnover is low overall."

# "In contrast... strong interactions... when adaptation is high... losses are
# the result of species interactions (85.2%)... As such, turnover is higher
# overall."

# -- supplementStatistics$STAT$extirpations
# -- supplementStatistics$STAT$turnover

supplementStatistics$inout <- tidytable::bind_rows(
  Pers |> tidytable::filter(
    # SpeciesPreferences == "100% 0",
    NicheDistance == defaultNicheDistance,
    # Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Persistence > 0,
    InType != externalNames["Dispersal"],
    In < Stop, Out > Start
  ) |> tidytable::group_by(
    SpeciesPreferences, InType, OutType, Intervention
  ) |> tidytable::summarise(
    Average = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
  ),
  ColExt |> tidytable::filter(
    # Filter what we can:
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    !Success | EventType == "Present"
  ) |> tidytable::left_join(
    # Start and Stop aren't already present in this version
    endTimes |> tidytable::select(-Times)
  ) |> tidytable::filter(
    Times > Start, Times < Stop
  ) |> tidytable::mutate(
    SpeciesPreferences =
      speciesAffinityDictionaryOrigin$SpeciesAffinities[as.numeric(SpeciesAffinity)]
  ) |> changePreferencesLevels(
    # ) |> tidytable::filter(
    #   SpeciesPreferences == "100% 0"
  ) |> tidytable::left_join(
    interventionStrings,
    by = c("PatchAffinity", "PoolPatch", "InterventionPatchType")
    # ) |> tidytable::filter(
    #   Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)")
  ) |> tidytable::mutate(
    InType = externalNames[
      ifelse(EventType == "Arrival", "Failed Arrival", "Present")
      ],
    OutType = externalNames["NA"]
  ) |> tidytable::group_by(
    SpeciesPreferences, InType, OutType, Intervention
  ) |> tidytable::summarise(
    Average = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
  )
) |> changeInterventionLevels(
)

supplementStatistics$STAT$extirpations <-
  supplementStatistics$inout |> tidytable::filter(
  OutType != "NA", OutType != "Persistent"
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal
) |> tidytable::summarise(
  Fraction = Average / sum(Average) # /44's cancel out.
)

# Note switch to dplyr to group_modify
supplementStatistics$STAT$turnover <-
  diversitiesTimeBC |> tidytable::filter(
    PoolPatchSeed %in% basePoolPatchSeeds,
    NicheDistance == defaultNicheDistance,
    grepl(x = Metric, pattern = "TimeBrayCurtis:", fixed = TRUE),
    is.na(Subset),
    Time > Start, Time < Stop # Not things outside of [Start, Stop]
  ) |> dplyr::group_by(
    # Metric groups
    Metric, Subset,
    # Initial groups
    PoolPatch, InitialConditions, Dispersal, DispersalParam, NicheDistance,
    SpeciesAffinity, PatchAffinity,
    # Intervention groups
    InterventionPatchType, InterventionTimeType, InterventionDispersal,
    InterventionNicheDistance,
    # Human-readable groups
    SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal
  ) |> dplyr::group_modify(
    .f = function(value, key) {
      cbind(#key,
              data.frame(rbind((summary(value$Value)))),
            skewness = moments::skewness(value$Value),
            kurtosis = moments::kurtosis(value$Value)
           # sum(value$Value == 0),
           # sum(value$Value < 0)
           # fitdistrplus::fitdist(value$Value + 1e-12, distr = "beta",
           #                       start = list(shape1 = 1, shape2 = 5))
           # ) # As nice as a beta would be, the 0's cause numerical problems.
      )
    }
  )

# Abundance (All, Basal, and Consumer) Changes: ###############################
# Going from no adaptation to intermediate, average basal abundance decreases by
# a factor of approximately 2, but comparing no adaptation to full adaptation
# shows an increase by a factor of approximately 4.

#TODO!!!

# Long term richness changes between different adaptation scenarios: ##########
# "The result is an overall reduction in richness, even in comparison to the
# extreme adaptation case (Figure 2), with an expected loss of 3.8 species
# (standard error 0.046), but this depends on the strength of the
# species-habitat interaction."

# "While there is some evidence of an edge effect in the land-use preferences --
# intermediate land-use had ~0.6 more species (differences of 0.621 and 0.629
# with standard errors of 0.0250 and 0.0248) -- the differences are otherwise
# minor in other traits."

# "While a positive increase in richness across simulations [with land-use
# change] can be consistently detected, the effect is small, approximately 0.5
# species compared to states that naturally vary on average by between ... and
# ... species."

# Short term richness changes during land-use change: #########################
# "Across the full set of parameter combinations discussed here, we observe
# positive richness changes in 3% of land-use change scenarios, compared to no
# change in 26% and declines in 71% of land-use change scenarios."


# #############################################################################

# diversitiesRichness |> tidytable::filter(
#   NicheDistance == defaultNicheDistance,
#   (PoolPatchSeed %in% as.character(343:386)),
#   Metric == "Alpha Hill:0",
#   is.na(Subset)
# ) |> tidytable::left_join(
#   endTimes |> dplyr::select(-Times)
# ) |> tidytable::group_by(
#   SpeciesAffinity, Intervention, PoolPatchSeed
# ) |> tidytable::arrange(
#   Time
# ) |> tidytable::filter(
#   InterventionInitial != InterventionFinal,
#   Time == Time[1] | Time == Time[2]
# ) |> tidytable::summarise(
#   InterventionChange = abs(
#     as.numeric(gsub(InterventionInitial, pattern = "[(]|[)]", replacement = ""))
#     - as.numeric(gsub(InterventionFinal, pattern = "[(]|[)]", replacement = ""))
#   ),
#   Time = Time[2] - Time[1],
#   Value = Value[2] - Value[1],
#   Method = "Temporal",
#   .groups = "drop"
# ) |> with(table(Intervention, sign(Value), SpeciesAffinity))
