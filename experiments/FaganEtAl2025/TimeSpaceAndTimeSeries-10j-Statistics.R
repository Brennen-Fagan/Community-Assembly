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
# extreme adaptation case (Figure 2), with an expected loss of 3.96 species
# (standard error 0.317), but this depends on the strength of the
# species-habitat interaction."

# This is comparing between Uniform(0, 1) and 100% 0 in (0).

supplementStatistics$diff_1000_5050 <-
  diversitiesRichness |> tidytable::filter(
    PoolPatchSeed %in% basePoolPatchSeeds,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)"),
    Time > Start, Time < Stop, # Not things outside of [Start, Stop]
    is.na(Subset),
    SpeciesPreferences %in% c("100% 0", "50% 0, 50% 1")
  ) |> tidytable::group_by(
    Time, PoolPatchSeed
  ) |> tidytable::mutate(
    Value = ifelse(SpeciesPreferences == "100% 0", # Reference
                   -Value, Value)
  ) |> tidytable::summarise(
    Value = sum(Value)
  )

supplementStatistics$STAT$diff_1000_5050 <- list(
  summary(supplementStatistics$diff_1000_5050$Value),
  sd(supplementStatistics$diff_1000_5050$Value),
  nlme::lme(
    data = supplementStatistics$diff_1000_5050,
    fixed = Value ~ 1, # Time, # Time ~ little to no effect.
    random = ~ 1 | PoolPatchSeed, # Nat. Variation in Intercept
    correlation = nlme::corAR1(form = ~ Time | PoolPatchSeed)
  ) |> summary()
)

# "While there is some evidence of an edge effect in the land-use preferences --
# intermediate land-use had ~0.6-0.7 more species (differences of 0.59 and 0.73
# with standard errors of 0.17 and 0.13) -- the differences are otherwise
# minor in other traits."

# This is comparing amongst Uniform(0, 1): (0.5) - (0) and (0.5) - (1).

supplementStatistics$diff_unif_5v0 <- tidytable::left_join(
  diversitiesRichness |> tidytable::filter(
    PoolPatchSeed %in% basePoolPatchSeeds,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(1)"),
    Time > Start, Time < Stop, # Not things outside of [Start, Stop]
    is.na(Subset),
    SpeciesPreferences %in% c("Uniform(0, 1)")
  ) |> tidytable::rename(
    Value1 = Value,
    Preferences1 = SpeciesPreferences,
    Intervention1 = Intervention
  ),
  diversitiesRichness |> tidytable::filter(
    PoolPatchSeed %in% basePoolPatchSeeds,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0.5)"),
    Time > Start, Time < Stop, # Not things outside of [Start, Stop]
    is.na(Subset),
    SpeciesPreferences %in% c("Uniform(0, 1)")
  ) |> tidytable::rename(
    Value2 = Value,
    Preferences2 = SpeciesPreferences
  ) |> tidytable::select(
    -tidytable::starts_with("Intervention"),
    -PatchAffinity, -PatchAffinitySeed, -Start, -Stop
  )
) |> tidytable::mutate(
  Value = round(Value2 - Value1, 4), # Prevent numerical errors
  Group = paste(Intervention1, PoolPatchSeed)
)

supplementStatistics$STAT$diff_unif_5v0 <-
  supplementStatistics$diff_unif_5v0 |> dplyr::group_by(
    Intervention1
  ) |> dplyr::group_map(
    .f = function(values, key) {
      list(
        key,
        summary(values$Value),
        sd(values$Value),
        nlme::lme(
          data = values,
          fixed = Value ~ 1, # Time, # Time ~ little to no effect.
          random = ~ 1 | Group, # Nat. Variation in Intercept
          correlation = nlme::corAR1(form = ~ Time | Group)
        ) |> summary()
      )
    }
  )


# "While a positive increase in richness across simulations [with land-use
# change] can be consistently detected, the effect is small, approximately 0.5
# species compared to states that naturally vary on average by between ... and
# ... species."

#TODO It's even lower than this. Probably should just delete.

supplementStatistics$baseStrings <- c(
  # Everything we need to match up the simulations with their natural contrast.
  "Time", "Environment1", "Environment2", "Metric", "Subset",
  "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
  "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
  "Dispersal", "NicheDistance", "SpeciesAffinity", "SpeciesAffinitySeed",
  #"PatchAffinity", "PatchAffinitySeed", # These are *starting* states.
  "SpeciesPreferences", "DispersalParam", "InterventionFinal"
)

# We pair by the corresponding case without an intervention, which otherwise
# has all parameters identical to the intervention case (naturally paired).
# We then take the difference (contrast) and try to "model" the contrasts over
# the set of interventions and species preferences.
# (Starting with summary statistics, then linear regression, etc. as needed.)
supplementStatistics$interventionContrast <- tidytable::left_join(
  # Cases with intervention
  diversitiesRichness |> tidytable::filter(
    PoolPatchSeed %in% basePoolPatchSeeds,
    NicheDistance == defaultNicheDistance,
    InterventionInitial != InterventionFinal,
    Time > Start, Time < Stop, # Not things outside of [Start, Stop]
    is.na(Subset)
  ) |> tidytable::rename(
    ValueIntervention = Value
  ),
  # Join to cases without intervention by ENDING state (NOT INITIAL)
  diversitiesRichness |> tidytable::filter(
    PoolPatchSeed %in% basePoolPatchSeeds,
    NicheDistance == defaultNicheDistance,
    InterventionInitial == InterventionFinal,
    Time > Start, Time < Stop, # Not things outside of [Start, Stop]
    is.na(Subset)
  ) |> tidytable::select(
    Time:SpeciesAffinitySeed, SpeciesPreferences, DispersalParam,
    InterventionFinal
  ) |> tidytable::rename(
    ValueNoIntervention = Value
  ),
  by = supplementStatistics$baseStrings
) |> tidytable::mutate(
  Value = ValueIntervention - ValueNoIntervention,
  Time = Time - 20000 # Scale so that it's time since start of comparison.
  # Switch to dplyr for group_map
) |> dplyr::group_by(
  SpeciesPreferences, InterventionFinal
) |> dplyr::group_map(
  .f = function(values, key) {
    list(
      key,
      summary(values$Value),
      nlme::lme(
        data = values,
        fixed = Value ~ Time, # Time ~ little to no effect.
        random = ~ 1 | Intervention/PoolPatchSeed, # Nat. Variation in Intercept
        correlation = nlme::corAR1(form = ~ Time | Intervention/PoolPatchSeed)
      )
    )
  }
)

# Short term richness changes during land-use change: #########################
# "Across the full set of parameter combinations discussed here, we observe
# positive richness changes in 3% of land-use change scenarios, compared to no
# change in 26% and declines in 71% of land-use change scenarios."

# ABOVE NEEDS TO BE ADJUSTED IN THE TEXT, e.g. we're not reaching the 26%.
# supplementStatistics$STAT$shortTermLoss

supplementStatistics$shortTermLoss <-
  diversitiesRichness |> tidytable::filter(
    PoolPatchSeed %in% basePoolPatchSeeds,
    NicheDistance == defaultNicheDistance,
    Metric == "Alpha Hill:0",
    InterventionInitial != InterventionFinal,
    is.na(Subset)
  ) |> tidytable::group_by(
    SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal,
    PoolPatchSeed
  ) |> tidytable::arrange(
    Time
  ) |> tidytable::summarise(
    # InterventionChange = abs(
    #   as.numeric(gsub(InterventionInitial, pattern = "[(]|[)]", replacement = ""))
    #   - as.numeric(gsub(InterventionFinal, pattern = "[(]|[)]", replacement = ""))
    # ),
    Time = round(Time - Time[2], digits = 4), # Make numerically safe.
    # Note the different conventions to make analysis easier
    PostIntervention = Time != Time[1],
    ValueDiffPreIntervention = round(Value - Value[1]), # Make numerically safe.
    .groups = "drop"
  ) |> tidytable::rename(
    TimeSinceIntervention = Time
  ) |> tidytable::filter(
    TimeSinceIntervention <= 51,
    TimeSinceIntervention == floor(TimeSinceIntervention)
  ) |> tidytable::group_by(
    SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal,
    TimeSinceIntervention
  ) |> tidytable::summarise(
    # Across PoolPatchSeeds
    Total = tidytable::n(),
    Neg = sum(ValueDiffPreIntervention < 0),
    Zero = sum(ValueDiffPreIntervention == 0),
    Pos = sum(ValueDiffPreIntervention > 0)
  )

supplementStatistics$PLOT$shortTermLoss <-
  supplementStatistics$shortTermLoss |> tidytable::pivot_longer(
    cols = Neg:Pos, names_to = "Type", values_to = "Counts"
  ) |> tidytable::mutate(
    Percentage = Counts / Total * 100,
    Text = paste0(formatC(Percentage, digits = 1, format = "f"), "%")
  ) |> ggplot2::ggplot(
    ggplot2::aes(x = TimeSinceIntervention, color = Intervention,
                 group = interaction(Type, Intervention, SpeciesPreferences))
  ) + ggplot2::geom_line(
    ggplot2::aes(y = Percentage, linetype = Type),
    show.legend = TRUE
  ) + ggplot2::geom_text(
    ggplot2::aes(y = Percentage, label = Text),
    show.legend = FALSE
  ) + ggplot2::facet_grid(
    SpeciesPreferences + InterventionInitial ~ InterventionFinal
  ) + ggplot2::scale_color_manual(
    values = colorPalette
  )

# Too fine to examine inside of R. Need to magnify further.
ggplot2::ggsave(
  plot = supplementStatistics$PLOT$shortTermLoss,
  filename = "Figure_supplementStatistics_ShortTermLoss.png",
  height = 200, width = 200, units = "cm", limitsize = FALSE
)

# Take some slices so we end up with something more writable in a paper.
supplementStatistics$STAT$shortTermLoss <-
  # 3 times for each of the 3 preferences with each of the 3 statistics.
  supplementStatistics$shortTermLoss |> tidytable::filter(
    TimeSinceIntervention %in% c(10, 20, 50)
  ) |> tidytable::group_by(
    SpeciesPreferences, TimeSinceIntervention
  ) |> tidytable::summarise(
    Total = sum(Total),
    Neg = sum(Neg),
    Zero = sum(Zero),
    Pos = sum(Pos)
  )

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
