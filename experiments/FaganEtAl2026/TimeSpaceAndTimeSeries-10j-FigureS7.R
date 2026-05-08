### Networks oriented: ######################################################
#### Setup: #################################################################
##### Resources: ############################################################
# Data:
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

load("TSTS_Interventions_10a1.RData")
InterventionTimes <- InterventionTimes |> tidytable::select(
  -tidytable::starts_with("Intervention")
) # Only for the base cases, so Intervention information all NA.

# Functions:
source(file.path("R", "flattenDiversity.R")) # Req'd by below
source(file.path("R", "generateNetworks.R")) # To create inset graphs.
library(patchwork) # Plot assembly
library(ggExtra) # N6 marginal distributions

##### Data Management: ######################################################
figureS7 <- list(
  graph = list(
    # Note that, unusually, we need all max(time), but only seed for others.
    seed = "2", # "11", "17", "2"!,                 # for examples
    time = c(100, 2000, 25000),                     # for examples
    timeInterventions = c(0, 10, 100, 1000, 10000), # for examples
    pref = c("100% 0", "Uniform(0, 1)"),            # for KDEs
    interventions = c("(0)", "(0.5)", "(1)",        # for KDEs
                      "(0)->(0.5)", "(0.5)->(0)")   # for examples
  ),
  interventions = c("(0)", "(0)->(0.5)", "(0.5)->(0)",
                    "(0.5)", "(0.5)->(1)", "(1)"),
  ci = 0.75
)

# Apply initial specification set-up, but then identify practical times
# post-intervention rather than the approximates supplied above for the
# examples.
figureS7$graph$specification <-
  diversitiesRichness |> tidytable::select(c(
    # Which network:
    "Time", "Environment1",
    # Which File (Base):
    "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
    "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
    "Dispersal", "NicheDistance",
    # Oops, there was a collision causing human readable to replace machine.
    # Will be replaced SpeciesAffinity#2 with -> SpeciesPreferences.
    "SpeciesAffinity",
    "SpeciesAffinitySeed", "PatchAffinity", "PatchAffinitySeed",
    # Which File (Intervention):
    "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
    "InterventionTimeSeed", "InterventionDispersal",
    "InterventionNicheDistance",
    # Ease of Use
    "SpeciesPreferences", "Intervention"
  )) |> tidytable::left_join(
    InterventionTimes |> tidytable::select(
      TimeIntervention, PoolPatch:PatchAffinitySeed
    ),
    by = c("PoolPatch", "PoolPatchSeed", "Interactions",
           "InteractionsSeed", "Events",
           "EventsSeed", "InitialConditions",
           "InitialConditionsSeed", "Dispersal",
           "NicheDistance", "SpeciesAffinity",
           "SpeciesAffinitySeed", "PatchAffinity",
           "PatchAffinitySeed")
  ) |> tidytable::mutate(
    TimeSinceIntervention = Time - TimeIntervention
  ) |> tidytable::filter(
    SpeciesPreferences %in% figureS7$graph$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% figureS7$graph$interventions | # for KDEs
      (PoolPatchSeed %in% figureS7$graph$seed &      # for examples
         Intervention %in% figureS7$interventions),  # for examples
    Time == max(figureS7$graph$time) |               # for KDEs
      PoolPatchSeed %in% figureS7$graph$seed         # for examples
  ) |> tidytable::distinct(
  )

figureS7$graph$timeInterventions <-
  figureS7$graph$specification |> tidytable::filter(
    PoolPatchSeed %in% figureS7$graph$seed
  ) |> with(
    TimeSinceIntervention[
      outer(
        TimeSinceIntervention,
        figureS7$graph$timeInterventions,
        function(x, y) abs(x - y)
      ) |> apply(2, which.min)
      ]
  )


# Why to the level of summary? Because the PlotMeanAndInner function
# isn't built to handle the multiple resolutions that we have in the
# actual data, which makes it harder to portray the data accurately.
figureS7$dataSummary <- tidytable::bind_rows(
  diversitiesRichness
) |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0)->(0.5)", # (0)->(0.5) used for N6
                      "(0.5)", "(0.5)->(0)", # all others used in N3.
                      "(0.5)->(1)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric %in% c("Alpha Hill:0", "Alpha Abundance", "TimeBrayCurtis: 10"),
  is.na(Subset) # Overall values
) |> tidytable::left_join(
  InterventionTimes
) |> tidytable::mutate(
  Metric = factor(
    Metric,
    levels = c("Alpha Hill:0", "Alpha Abundance", "TimeBrayCurtis: 10"),
    labels = c("Richness", "Abundance", "Bray-Curtis Dissimilarity"),
    ordered = TRUE
  ),
  Time = round(Time - TimeIntervention, 6) # remove false differences
) |> tidytable::filter(
  Time < 16050, # Need the start for the inset.
  # Avoid singletons.
  abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0,
  Time > -1000
) |> tidytable::mutate(
  Time = tidytable::case_when( # Create groupings for times.
    Time < -5 ~ round(Time, -1),
    # In between sampling regimes vs Initial Conditions
    Time < 0 & (InterventionInitial == InterventionFinal) ~ -5,
    Time < 0 & (InterventionInitial != InterventionFinal) ~ -0.1,
    Time <= 50 ~ round(Time, 0), # After 50, jumps to 10s
    Time < 1115 ~ round(Time, -1), # After 1115, slowly expands until ~1863
    Time < 1117 ~ 1110, # Marginal
    Time < 1298 ~ (round(Time/11)*11), # ~1300
    Time < 1395 ~ (round(Time/15)*15), # ~1400
    Time < 1500 ~ (round(Time/30)*30),
    Time < 1600 ~ (round(Time/50)*50),
    Time < 1831 ~ (round(Time/60)*60),
    Time < 16350 ~ round(Time, -2), # 100 gap hols until the end.
    TRUE ~ Time
  )
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, Time
) |> tidytable::summarise(
  Value = median(Value), .groups = "drop"
) |> tidytable::group_by(
  # Average across simulations
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, Time
) |> tidytable::summarise(
  Lower = quantile(
    Value,
    p = (1 - figureS7$ci) - (1 - figureS7$ci)/2,
    na.rm = TRUE
  ),
  Average = mean(Value),
  Upper = quantile(
    Value,
    p = figureS7$ci + (1 - figureS7$ci)/2,
    na.rm = TRUE
  )
)

##### Functions: ############################################################
figureS7$renameSpeciesPreferences <- function(dat) {
  dat |> tidytable::mutate(
    SpeciesPreferences = tidytable::case_when(
      SpeciesPreferences == "100% 0" ~ "Single Adaptation Type",
      SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
      SpeciesPreferences == "Uniform(0, 1)" ~ "Many Adaptation Types",
      TRUE ~ SpeciesPreferences
    )
  )
}

### Figure: ##################################################################
# Two variants:
#   Time series of the RAT (no C, would need to massively refactor code).
#     Unfortunately, due to holes in my sampling strategy, BrayCurtis
#     doesn't have all the points necessary to look convincing.
#       (1, 2, & 3 end early, 10 has holes, and these are substitutable.)
#   Progression of the networks with R only.
#
###### Prep KDE Data: #####################################################
# Goal: x = "Affinity", y = "Size", scatter plot at a time, with marginal
# KDE plots, so that we can show how the distribution of species shifts
# after the intervention for a system with uniform habitat adaptations.

ColExtInterventionStrings <- ColExt %>% dplyr::select(
  PatchAffinity, PoolPatch, InterventionPatchType
) %>% dplyr::distinct(
) %>% dplyr::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    PatchAffinity, PoolPatch, InterventionPatchType
  ))
)

figureS7$ColExt <- ColExt |> tidytable::filter(
  EventType != "Present", Success, # Only Events that Did Something.
  NicheDistance == defaultNicheDistance # Base case.
) |> tidytable::mutate( # Make human readable formats.
  SpeciesPreferences =
    speciesAffinityDictionaryOrigin$SpeciesAffinities[
      as.numeric(SpeciesAffinity)
      ]
) |> tidytable::left_join(
  ColExtInterventionStrings,
  by = c("PatchAffinity", "PoolPatch", "InterventionPatchType"),
  multiple = "all"
) |> changePreferencesLevels(
) |> changeInterventionLevels(
) |> tidytable::filter( # Reduce to the levels we're interested in.
  SpeciesPreferences == "Uniform(0, 1)",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)",  # To display.
                      "(0)", "(0.5)") # To sense check.
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::left_join( # So we can calculate our post-int. slices
  InterventionTimes
) |> tidytable::group_by(
  # Species Characteristics
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed,
  Affinity,
  # Simulation Characteristics
  PoolPatch:TimeIntervention
) |> tidytable::mutate( # Create the groupings for when present in system.
  In = ifelse(EventType == "Arrival", "In", "Out"),
  InOutNumber = cumsum(In == "In"),
  Time = Times - TimeIntervention
) |> tidytable::pivot_wider( # Rearrange so we have time ranges
  names_from = In, values_from  = Time,
  id_cols = c(
    # Species Characteristics
    Species, Environment, SpeciesType, Size, ReproductionRate, Speed,
    Affinity,
    # Simulation Characteristics
    PoolPatch:TimeIntervention,
    InOutNumber
  )
) |> tidytable::filter( # Remove ranges that aren't included in any slice.
  outer( # Check all combinations intervention and range, then summarize.
    FUN = function(X, Y, low, high) {low[X] < Y & Y < high[X]},
    X = seq_along(In),
    Y = figureS7$graph$timeInterventions,
    low = In, high = Out
  ) |> apply(MARGIN = 1, FUN = any)
)


###### Main Plot: #########################################################
figureS7$plotA <- lapply(
  unique(figureS7$dataSummary$Metric), function(met, dat) {
    plt <- ggplot2::ggplot(
      dat |> tidytable::filter(
        Time >= -100, Metric == met,
        Intervention %in% c("(0)", "(0.5)->(0)", "(0.5)", "(0)->(0.5)"),
        SpeciesPreferences == "Uniform(0, 1)"
      ) |> figureS7$renameSpeciesPreferences(
      ),
      aes(
        x = Time, y = Average,
        color = Intervention,
        fill = Intervention,
        group = interaction(Intervention, Metric),
        alpha = Intervention
      )
    ) + ggplot2::geom_vline(
      xintercept = 0, color = "black", linetype = "dashed"
    ) + ggplot2::geom_line(
    ) + ggplot2::scale_color_manual(
      values = colorPalette, aesthetics = c("color", "fill"),
      name = ""
    ) + ggplot2::guides(
      linetype = "none",
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      legend.position = if(met == "Richness") {c(0.05, 0.35)} else {"none"},
      panel.spacing = ggplot2::unit(1, "lines"),
      strip.text = ggplot2::element_text(size = 12)
    ) + ggplot2::labs(
      x = "Time Since Intervention",
      y = met
    ) + ggplot2::coord_cartesian(
      xlim = c(-20, NA)
    ) + ggplot2::scale_x_continuous(
      breaks = c(0, 1, 10, 100, 1000, 10000, 15000),
      transform = scales::transform_pseudo_log(sigma = 10)
    ) + ggplot2::scale_alpha_manual(
      values = c("(0)" = 1, "(0.5)->(0)" = 1,
                 "(0.5)" = 1, "(0)->(0.5)" = 1),
      guide = "none"
    )
    return(plt)
  },
  dat = figureS7$dataSummary
)

###### Network Plots: ####################################################
figureS7$plotB <- with(
  list(
    plotfun = function(inter, time) {
      list(
        (
          figureS7$ColExt |> tidytable::filter(
            In < time, time < Out, Intervention == inter
          ) |> ggplot2::ggplot(
            ggplot2::aes(
              x = Affinity, y = Size, color = SpeciesType
            )
          ) + ggplot2::geom_point(
            show.legend = FALSE
          ) + ggplot2::scale_y_log10(
          ) + ggplot2::scale_x_continuous(
            breaks = c(0, 0.5, 1)
          ) + ggplot2::theme_minimal(
          ) + ggplot2::labs(
            x = if(time == 0) "Habitat Adaptation",
            y = if(time == 0) "Species Size"
          ) + scale_color_manual(
            values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
          ) + ggplot2::coord_cartesian(
            xlim = c(0, 1), ylim = c(10^-2, 10^0.5), expand = FALSE
          )
        ) |> ggExtra::ggMarginal(
        )
      )
    }
  ),
  expand_grid(
    Intervention = sort(unique(figureS7$ColExt$Intervention)),
    Time = figureS7$graph$timeInterventions
  ) |> tidytable::mutate_rowwise(
    Plot = plotfun(Intervention, Time)
  )
)

###### Save: #############################################################
for (i in 1+5*(0:3)) {
  ggplot2::ggsave(
    # Use Patchwork to Combine
    with(figureS7,
         plotA[[1]] + ggplot2::geom_vline(
           xintercept = plotB$Time[i+0:4], linetype = "dashed"
         )  + ggplot2::scale_alpha_manual(
           values = (table(plotB$Intervention[i+0:4])+1)/6,
           guide = "none"
         ) +
           plotB$Plot[[i]] + plotB$Plot[[i+1]] + plotB$Plot[[i+2]] +
           plotB$Plot[[i+3]] + plotB$Plot[[i+4]] +
           patchwork::plot_layout(
             design =
               "BBCCDDEEFF
                BBCCDDEEFF
                AAAAAAAAAA
                AAAAAAAAAA"
           )
    ),
    path = dirImages,
    filename = paste0("FigureS7_RichnessAndKDEs_", i, ".png"),
    units = "cm", width = 20, height = 11
  )
}