# Plotting with both ColExt and Divs
load("ColExt9a9_flat.RData")
load("diversitiesFlattened9a9_subset2.RData")

# Problems with X11
options(bitmapType = "cairo")

# Grey interval that we compute over, usually after intervention (~50%)
# If second number is less than 1, we lose persistent species.
end <- c(0.6, 1)

# Libraries: ##################################################################
source("TimeSpaceAndTimeSeries-9-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')

library(ggplot2)
library(ggpubr)
library(tidytable)

# colors: #####################################################################
#                                (0),   (0.5),   (1)
colorPalette <- c(#              Cyan, Magenta, Yellow, Black
  "(0)" =        DescTools::CmykToRgb(1,   0,   0,   0.25),
  "(0)->(0.5)" = DescTools::CmykToRgb(1,   0.5, 0,   0.25),
  "(0)->(1)" =   DescTools::CmykToRgb(1,   0,   0.75, 0.25),
  "(0.5)" =      DescTools::CmykToRgb(0,   1,   0,   0.25),
  "(0.5)->(0)" = DescTools::CmykToRgb(0.5, 1,   0,   0.25),
  "(0.5)->(1)" = DescTools::CmykToRgb(0,   1,   0.75, 0.25),
  "(1)" =        DescTools::CmykToRgb(0,   0,   1,   0.25),
  "(1)->(0)" =   DescTools::CmykToRgb(0.5, 0,   1,   0.25),
  "(1)->(0.5)" = DescTools::CmykToRgb(0,   0.5, 1,   0.25)
)

linetypePalette <- c(
  "100% 0" = "solid",
  "50% 0, 50% 1" = "longdash",
  "Uniform(0, 1)" = "dotted"
)

# renames: ####################################################################
# For presentation (i.e., "Arrival" is a working term, but not 100% accurate.)
externalNames <- c(
  "Arrival"         = "Colonisation",
  "Failed Arrival"  = "Failure",
  "Present"         = "Present",
  "Dispersal"       = "Adjacent",
  "Extinct"         = "Neutral Ext.",
  "Dynamic Loss"    = "Dynamic Ext.",
  "EndOfSimulation" = "Persistent",
  "NA"              = "NA"
)

# Functions: ##################################################################
### Manipulation: #############################################################
changeAffinityLevels <- function(df) {
  df %>% tidytable::mutate(
    SpeciesAffinity = tidytable::case_when(
      SpeciesAffinity == "rep_0" ~ "100% 0",
      SpeciesAffinity == "evensplit_01" ~ "50% 0, 50% 1",
      SpeciesAffinity == "runif" ~ "Uniform(0, 1)",
      TRUE ~ SpeciesAffinity
    ),
    SpeciesAffinity = factor(SpeciesAffinity, levels = c(
      "100% 0", "50% 0, 50% 1", "Uniform(0, 1)"
    ), ordered = TRUE)
  )
}

changeInterventionLevels <- function(df) {
  df %>% tidytable::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0.5)", "(0)->(1)",
        "(0.5)->(0)", "(0.5)", "(0.5)->(1)",
        "(1)->(0)", "(1)->(0.5)", "(1)"
      ), ordered = TRUE
    )
  )
}

### Plotting: #################################################################
plotMeanAndInner <- function(
  data, CIs = c(0.5, 0.95),
  facets = as.formula(Intervention ~ SpeciesAffinity)
) {
  data$Subset <- ifelse(is.na(data$Subset), "NA", data$Subset)
  baseplot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Time, y = Value,
      group = interaction(
        Subset,
        Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
      ),
      color = Intervention, fill = Intervention, linetype = SpeciesAffinity
    )
  )
  for (CI in CIs) {
    baseplot <- baseplot + ggplot2::geom_ribbon(
      data = data %>% tidytable::mutate(
        Time = round(Time, digits = -2)
      ) %>% tidytable::group_by(
        Time, Subset,
        Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
      ) %>% tidytable::summarise(
        top = quantile(Value, probs = CI+(1-CI)/2, na.rm = TRUE),
        bot = quantile(Value, probs = (1-CI)-(1-CI)/2, na.rm = TRUE)
      ), mapping = ggplot2::aes(
        x = Time, ymin = bot, ymax = top,
        group = interaction(
          Subset,
          Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
        ),
        fill = Intervention,
        color = Intervention
      ), inherit.aes = FALSE,
      alpha = 0.25, linewidth = 0.25 #, linetype = "dotted"
    )
  }
  baseplot <- (baseplot
               # ) + ggplot2::geom_smooth(
               # color = "black"
  ) + ggplot2::geom_line(
    data = data %>% tidytable::mutate(
      Time = round(Time, digits = -2)
    ) %>% tidytable::group_by(
      Time, Subset,
      Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
    ) %>% tidytable::summarise(
      Value = median(Value, na.rm = TRUE)
    )
  ) + ggplot2::facet_grid(
    facets
  ) + ggplot2::scale_color_manual(
    values = colorPalette, aesthetics = c("color", "fill"),
    name = "Island Land-use"
  ) + ggplot2::scale_linetype_manual(
    name = "Species Preferences",
    values = linetypePalette
  ) + ggplot2::theme_minimal(
  ) + ggplot2::guides(
    color = "none",
    fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
  )
  return(baseplot)
}

plotValueChart <- function(
  data, facets = as.formula(Intervention ~ SpeciesAffinity)
) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = interaction(InType,
                                 OutType,
                                 sep = "\n"),#"",
                 y = ChartValue,
                 fill = interaction(InType,
                                    OutType,
                                    sep = "/"))
  ) + ggplot2::geom_col(
    show.legend = FALSE
  ) + ggplot2::facet_grid(
    facets
  ) + ggplot2::theme_minimal(
  ) + ggplot2::labs(
    x = "Enter/Exit Type", y = "Count"#, fill = "In/Out"
  )
}

# Primary figure generator for Figures 2 and 3 currently.
createOverviewFigure <- function(
  .divs, # diversitiesAll
  .ps, # Pers
  .ces, # ColExt
  .ets, # endTimes
  ... # commonfilters
) {
  ggpubr::ggarrange(
    # Richness
    plotMeanAndInner(
      .divs %>% tidytable::filter(
        ...,
        Metric == "Alpha Hill:0",
        is.na(Subset)
      ),
      CIs = 0.75,
      facets = as.formula(. ~ .)
    ) + ggplot2::labs(
      y = "Richness"
    ) + ggplot2::theme(
      legend.position = "bottom", legend.direction = "vertical"
    ) + ggplot2::geom_rect(
      data = data.frame(
        1 # 1 rectangle per row, so dummy df to prevent overplotting
      ),
      xmin = min(.ets$Start),
      xmax = max(.ets$Stop),
      ymin = 0, ymax = 45,
      fill = "grey",
      alpha = 0.2,
      inherit.aes = FALSE
    ),
    plotValueChart(
      rbind(
        .ps %>% tidytable::filter(
          ...,
          Persistence > 0,
          InType != externalNames["Dispersal"],
          In < Stop, Out > Start
        ) %>% tidytable::group_by(
          SpeciesAffinity, InType, OutType, Intervention
        ) %>% tidytable::summarise(
          ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
        ) %>% dplyr::mutate( # Tidytable renders as character again!
          Intervention =
            factor(Intervention,
                   levels = rev(c("(0)", "(0)->(0.5)", "(0)->(1)",
                                  "(0.5)",
                                  "(1)")),
                   ordered = TRUE)
        ),
        .ces %>% tidytable::filter(
          ...,
          !Success | EventType == "Present",
          Times > Start, Times < Stop
        ) %>% tidytable::mutate(
          InType = externalNames[
            ifelse(EventType == "Arrival", "Failed Arrival", "Present")
            ],
          OutType = externalNames["NA"]
        ) %>% tidytable::group_by(
          SpeciesAffinity, InType, OutType, Intervention
        ) %>% tidytable::summarise(
          ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
        ) %>% dplyr::mutate( # Tidytable renders as character again!
          Intervention =
            factor(Intervention,
                   levels = rev(c("(0)", "(0)->(0.5)", "(0)->(1)",
                                  "(0.5)",
                                  "(1)")),
                   ordered = TRUE)
        )
      ),
      facets = as.formula(Intervention ~ .)
    ) + ggplot2::theme(
      legend.position = "bottom"
    ) + ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 2)
    ) + ggplot2::ylab(
      "Average Count in a Simulation"
    ),
    ggplot2::ggplot(
      .ps %>% tidytable::filter(
        ...,
        Persistence > 0,
        In < Stop, Out > Start
      ),
      ggplot2::aes(
        x = Size,
        # color = AffinityBins,
        y = Persistence
      )
    ) + ggplot2::geom_hex(
    ) + ggplot2::theme_minimal(
    ) + ggplot2::scale_x_log10(
    ) + ggplot2::scale_y_log10(
    ) + ggplot2::scale_fill_viridis_c(
      trans = "log10", direction = -1
    ) + ggplot2::theme(
      legend.position = "bottom"
    ) + ggplot2::facet_grid(
      factor(Intervention, levels = rev(levels(Intervention))) ~ .
    ) + ggplot2::geom_vline(
      xintercept = 0.10, color = "red"
    ),
    nrow = 1, widths = c(1, 0.6, 0.6)#, common.legend = TRUE
  )
}

# START HERE: ##################################################################

# Enhance readability, from 9g TablePlots
diversitiesInterventionStrings <- ColExt %>% tidytable::select(
  Affinity, PoolPatch, InterventionPatchType
) %>% tidytable::distinct(
) %>% tidytable::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    Affinity, PoolPatch, InterventionPatchType
  ))
)

# Work out the end times so we can truncate the simulations
# so that we are making sure our comparisons are equivalent.
endTimes <- ColExt %>% tidytable::rename(
  DispersalParam = Dispersal
) %>% tidytable::filter(
  EventType == "EndOfSimulation"
) %>% tidytable::select(
  Times, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, Affinity, AffinitySeed
) %>% tidytable::distinct(
) %>% tidytable::group_by(
  # One of these had an early stop. We "fix" it by going to its descendants.
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, Affinity, AffinitySeed
) %>% tidytable::summarise(
  Times = max(Times),
  .groups = "drop"
) %>% tidytable::mutate( # In the plots:
  Start = end[1] * Times, # Neglect anything with an out time before this.
  Stop = end[2] * Times # Neglect anything with an in time after this.
)

Pers <- ColExt %>% tidytable::rename(
  DispersalParam = Dispersal
) %>% tidytable::filter(
  EventType != "Present", Success # False Arrivals might mess this up.
) %>% tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed, Affinity,
  AffinityBins, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed,
  Events, EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance
) %>% tidytable::mutate(
  InNumber = ifelse(EventType == "Arrival" | EventType == "Dispersal", 1, 0),
  InNumber = cumsum(InNumber)
) %>% tidytable::ungroup(
) %>% tidytable::pivot_wider(
  values_from = "Times",
  names_from = EventType,
  id_cols = c(
    # Species/Environment/Simulation Identifiers
    Species, Environment, SpeciesType, Size, ReproductionRate, Speed, Affinity,
    AffinityBins, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed,
    Events, EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
    NicheDistance, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
    InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
    InterventionNicheDistance,
    # And the true ID variable!
    InNumber
  ),
  values_fill = NA
) %>% tidytable::mutate(
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
  SpeciesAffinity =
    affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
) %>% changeAffinityLevels(
) %>% tidytable::select(
  -Dispersal, -Arrival, -Extinct, -`Dynamic Loss`, -EndOfSimulation
) %>% tidytable::left_join(
  diversitiesInterventionStrings,
  by = c("Affinity", "PoolPatch", "InterventionPatchType"),
  multiple = "all"
) %>% changeInterventionLevels(
) %>% tidytable::left_join(
  endTimes
)

diversitiesAll <- diversitiesAll %>% changeAffinityLevels(
) %>% changeInterventionLevels(
) %>% tidytable::mutate(
  #TODO consider moving this into changeInterventionLevels.
  InterventionInitial = tidytable::case_when(
    Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)") ~ "(0)",
    Intervention %in% c("(0.5)->(0)", "(0.5)", "(0.5)->(1)") ~ "(0.5)",
    Intervention %in% c("(1)->(0)", "(1)->(0.5)", "(1)") ~ "(1)",
    TRUE ~ NA_character_
  ),
  InterventionInitial = factor(
    InterventionInitial,
    levels = c(
      "(0)", "(0.5)",  "(1)"
    ), ordered = TRUE
  ),
  InterventionFinal = tidytable::case_when(
    Intervention %in% c("(0)", "(0.5)->(0)", "(1)->(0)") ~ "(0)",
    Intervention %in% c("(0)->(0.5)", "(0.5)",  "(1)->(0.5)") ~ "(0.5)",
    Intervention %in% c("(0)->(1)", "(0.5)->(1)", "(1)") ~ "(1)",
    TRUE ~ NA_character_
  ),
  InterventionFinal = factor(
    InterventionFinal,
    levels = c(
      "(0)", "(0.5)",  "(1)"
    ), ordered = TRUE
  )
)

ColExt <- ColExt %>% tidytable::rename(
  DispersalParam = Dispersal
) %>% tidytable::mutate(
  SpeciesAffinity =
    affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
) %>% changeAffinityLevels(
) %>% tidytable::left_join(
  diversitiesInterventionStrings,
  by = c("Affinity", "PoolPatch", "InterventionPatchType"),
  multiple = "all"
) %>% changeInterventionLevels(
) %>% tidytable::left_join(
  endTimes %>% tidytable::select(-Times),
  by = c(
    "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed", "Events",
    "EventsSeed", "InitialConditions", "InitialConditionsSeed",
    "DispersalParam", "NicheDistance", "Affinity", "AffinitySeed"
  )
)

# Define Species Color Scale from AffinityBins: ##############################
# Can't be done before hand because the affinity bins are calculated differently
# for each of the runs.
SpeciesPaletteDF <- data.frame(
  Bins = sort(levels(Pers$AffinityBins))
) %>% tidytable::separate_wider_delim(
  Bins, delim = c(","),
  cols_remove = FALSE
) %>% tidytable::mutate(
  Bins1 = gsub(pattern = "(", replacement = "", Bins1, fixed = TRUE),
  Bins2 = gsub(pattern = "]", replacement = "", Bins2, fixed = TRUE),
  Bins2 = ifelse(is.na(Bins2), Bins1, Bins2),
  Bins1 = as.numeric(Bins1),
  Bins2 = as.numeric(Bins2),
  ColorNumber = (Bins1 + Bins2)/2,
  color =
    # Define the Color Palette to reasonable resolution
    colorRampPalette(colorPalette[c("(0)", "(0.5)", "(1)")])(length(Bins))[
      # Apply Color Palette judiciously, stackoverflow.com/a/35656396
      findInterval(ColorNumber, seq(0, 1, length.out = length(Bins)))
      ]
) %>% tidytable::arrange(
  Bins
)

SpeciesPalette <- c(SpeciesPaletteDF$color)
names(SpeciesPalette) <- SpeciesPaletteDF$Bins

# Plots: #####################################################################
plot1 <- createOverviewFigure(
  .divs = diversitiesAll,
  .ps = Pers,
  .ces = ColExt,
  .ets = endTimes,
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  !(PoolPatchSeed %in% c("341", "342"))
)

ggplot2::ggsave(
  plot1,
  filename = "Figure2_Prototype5.png",
  units = "px", height = 3200, width = 4800 # Words are too small
  # height = 1600, width = 2400 # Words over-run
)

plot2 <- createOverviewFigure(
  .divs = diversitiesAll,
  .ps = Pers,
  .ces = ColExt,
  .ets = endTimes,
  SpeciesAffinity == "100% 0",
  NicheDistance == "5",
  Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)"),
  !(PoolPatchSeed %in% c("341", "342"))
)

ggplot2::ggsave(
  plot2,
  filename = "Figure3_Prototype3.png",
  units = "px", height = 3200, width = 4800 # Words are too small
  # height = 1600, width = 2400 # Words over-run
)
