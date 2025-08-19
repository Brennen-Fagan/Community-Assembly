library(ggplot2)
library(dplyr)

# dataframes: #################################################################
load("div_2024-11-30-341.RData")
load("pres_2024-11-30-341.RData")

divTemp <- divTemp %>% dplyr::mutate(
  SpeciesAffinity = dplyr::case_when(
    SpeciesAffinity == "rep_0" ~ "100% 0",
    SpeciesAffinity == "evensplit_01" ~ "50% 0, 50% 1",
    SpeciesAffinity == "runif" ~ "Uniform(0, 1)",
    TRUE ~ SpeciesAffinity
  ),
  SpeciesAffinity = factor(SpeciesAffinity, levels = c(
    "100% 0", "50% 0, 50% 1", "Uniform(0, 1)"
  ), ordered = TRUE),
  Intervention = factor(
    Intervention,
    levels = c(
      "(0)", "(0)->(0.5)", "(0)->(1)",
      "(0.5)->(0)", "(0.5)", "(0.5)->(1)",
      "(1)->(0)", "(1)->(0.5)", "(1)"
    ), ordered = TRUE
  ),
  InterventionInitial = dplyr::case_when(
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
  InterventionFinal = dplyr::case_when(
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
presTemp <- presTemp %>% dplyr::mutate(
  SpeciesAffinity = dplyr::case_when(
    SpeciesAffinity == "rep_0" ~ "100% 0",
    SpeciesAffinity == "evensplit_01" ~ "50% 0, 50% 1",
    SpeciesAffinity == "runif" ~ "Uniform(0, 1)",
    TRUE ~ SpeciesAffinity
  ),
  SpeciesAffinity = factor(SpeciesAffinity, levels = c(
    "100% 0", "50% 0, 50% 1", "Uniform(0, 1)"
  ), ordered = TRUE),
  Intervention = factor(
    Intervention,
    levels = c(
      "(0)", "(0)->(0.5)", "(0)->(1)",
      "(0.5)->(0)", "(0.5)", "(0.5)->(1)",
      "(1)->(0)", "(1)->(0.5)", "(1)"
    ), ordered = TRUE
  ),
  InterventionInitial = dplyr::case_when(
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
  InterventionFinal = dplyr::case_when(
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

## Columns: ###################################################################
### Both: #####################################################################
# "Time": The time of the record of the diversity or presence.
# "PoolPatch": Internal code controlling the pool and patch properties.
# "PoolPatchSeed": As above, but randomisation. Two independent runs (341, 342).
# "Interactions": As above, but interactions properties.
# "InteractionsSeed": I suspect you get the idea.
# "Events"
# "EventsSeed"
# "InitialConditions"
# "InitialConditionsSeed"
# "Dispersal": No dispersal allowed.
# "NicheDistance": I've sent the most intense (7 => 10x (dis)advantage).
#                  I have two weaker levels (2x and 5x).
# "Affinity"
# "AffinitySeed"
# "InterventionPatchType"
# "InterventionPatchSeed"
# "InterventionTimeType"
# "InterventionTimeSeed"
# "InterventionDispersal"
# "InterventionNicheDistance"
# "Intervention": The human readable version of patch type and intervention.
# "SpeciesAffinity": The more human readable version of pool type.

### divTemp (Diversity Data <- diversitiesFlattened): #########################
# "Environment1": Only 1 Environment at the moment.
# "Environment2": NA as above. (In case we do beta comparisons.)
# "Metric": I've only sent through Richness (Alpha Hill:0).
# "Value": Richness values.
# "Subset": I've only sent through overall richness.

### presTemp (Presence/Species Data <- presencesFlattened): ###################
# "Species": Species ID.
# "Abundance": Density of individuals in the environment.
# "Environment": Only 1 Environment at the moment.
# "Size": Effectively the mass of the species (IIRC!).
# "Type": Basal or Consumer archetype.
# "Affinity": The preferred habitat type for which the species receives a boost
# "TimeGroup": Easy way of selecting three time slices to pull out.

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

# plots: ######################################################################
## Diversity: #################################################################
### Richness: #################################################################
# Line plot of richness.
ImageRichness <- ggplot2::ggplot(
  divTemp,
  ggplot2::aes(x = Time, y = Value, color = Intervention,
               linetype = InterventionPatchType)
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::labs(
  y = "Richness"
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::guides(
  linetype = "none"
) + ggplot2::theme_minimal(
)

## Presence: ##################################################################
# Dot plots.
### Abundance: ################################################################
ImageAbundance <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    !is.na(TimeGroup),
    PoolPatchSeed == PoolPatchSeed[1]
  ),
  ggplot2::aes(x = Time, y = Abundance, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_point(
  position = ggplot2::position_dodge(width = 6000), alpha = 0.5, size = 1
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

### Size: #####################################################################
ImageSize <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    !is.na(TimeGroup),
    PoolPatchSeed == PoolPatchSeed[1]
  ),
  ggplot2::aes(x = Time, y = Size, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_point(
  position = ggplot2::position_dodge(width = 6000), alpha = 0.5, size = 1
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

### Affinity: #################################################################
# Note the use of count rather than point to deal with discrete overplotting.
ImageAffinity <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    !is.na(TimeGroup),
    PoolPatchSeed == PoolPatchSeed[1]
  ),
  ggplot2::aes(x = Time, y = Affinity, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_count(
  position = ggplot2::position_dodge(width = 6000), alpha = 0.6
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::guides(
  linetype = "none"
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)


### Abundance Line: ###########################################################
ImageAbundance2_Geometric <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Time, Intervention, SpeciesAffinity
  ) %>% dplyr::summarise(
    Abundance = exp(mean(log(Abundance))),
    .groups = "drop"
  ),
  ggplot2::aes(x = Time, y = Abundance, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)
ImageAbundance2 <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Time, Intervention, SpeciesAffinity
  ) %>% dplyr::summarise(
    Abundance = mean(Abundance),
    .groups = "drop"
  ),
  ggplot2::aes(x = Time, y = Abundance, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

ImageAbundance3 <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Intervention, Species, SpeciesAffinity
  ) %>% dplyr::mutate(
    Run = cumsum(c(0, na.omit(Time > lag(Time)+10)))
  ),
  ggplot2::aes(x = Time, y = Abundance, color = Intervention,
               group = interaction(Intervention, Species, Run))
) + ggplot2::geom_line(
  alpha = 0.5
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  Intervention ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

### Size Line: ################################################################
ImageSize2_Geometric <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Time, Intervention, SpeciesAffinity
  ) %>% dplyr::summarise(
    Size = exp(mean(log(Size))),
    .groups = "drop"
  ),
  ggplot2::aes(x = Time, y = Size, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)
ImageSize2 <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Time, Intervention, SpeciesAffinity
  ) %>% dplyr::summarise(
    Size = mean(Size),
    .groups = "drop"
  ),
  ggplot2::aes(x = Time, y = Size, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

ImageSize3 <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Intervention, Species, SpeciesAffinity
  ) %>% dplyr::mutate(
    Run = cumsum(c(0, na.omit(Time > lag(Time)+10)))
  ),
  ggplot2::aes(x = Time, y = Size, color = Affinity,
               group = interaction(Intervention, Species, Run))
) + ggplot2::geom_line(
  alpha = 0.5
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  Intervention ~ SpeciesAffinity,
  scales = "free"
# ) + ggplot2::scale_color_manual(
  # values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

### Affinity Line: ############################################################

ImageAffinity2 <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Time, Intervention, SpeciesAffinity
  ) %>% dplyr::summarise(
    Affinity = mean(Affinity),
    .groups = "drop"
  ),
  ggplot2::aes(x = Time, y = Affinity, color = Intervention,
               group = interaction(Intervention))
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

ImageAffinity3 <- ggplot2::ggplot(
  presTemp %>% dplyr::filter(
    PoolPatchSeed == PoolPatchSeed[1]
  ) %>% dplyr::group_by(
    Intervention, Species, SpeciesAffinity
  ) %>% dplyr::mutate(
    Run = cumsum(c(0, na.omit(Time > lag(Time)+10)))
  ),
  ggplot2::aes(x = Time, y = Affinity, color = Intervention,
               group = interaction(Intervention, Species, Run))
) + ggplot2::geom_line(
  alpha = 0.5
) + ggplot2::facet_grid(
  # Row:
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  # Column:
  Intervention ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)


## Combine: ###################################################################
# Vertical arrangement
ggpubr::ggarrange(
  # Richness, no x description.
  ImageRichness + ggplot2::theme(axis.title.x = ggplot2::element_blank()),
  # Abundance, no x description and no facet labels.
  ImageAbundance + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                  strip.background = ggplot2::element_blank(),
                                  strip.text = ggplot2::element_blank()),
  # Size, no x description and no facet labels.
  ImageSize + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                             strip.background = ggplot2::element_blank(),
                             strip.text = ggplot2::element_blank()),
  # Affinity, no facet labels.
  ImageAffinity + ggplot2::theme(strip.background = ggplot2::element_blank(),
                                 strip.text = ggplot2::element_blank()),
  # Common legend, but note that counts are not carried over.
  # Align forces the plotting areas to match (across the column).
  common.legend = TRUE, legend = "right", ncol = 1, align = "v")

ggpubr::ggarrange(
  # Richness, no x description.
  ImageRichness + ggplot2::theme(axis.title.x = ggplot2::element_blank()),
  # Abundance, no x description and no facet labels.
  ImageAbundance2_Geometric + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                  strip.background = ggplot2::element_blank(),
                                  strip.text = ggplot2::element_blank()),
  # Size, no x description and no facet labels.
  ImageSize2_Geometric + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                             strip.background = ggplot2::element_blank(),
                             strip.text = ggplot2::element_blank()),
  # Affinity, no facet labels.
  ImageAffinity2 + ggplot2::theme(strip.background = ggplot2::element_blank(),
                                 strip.text = ggplot2::element_blank()),
  # Common legend, but note that counts are not carried over.
  # Align forces the plotting areas to match (across the column).
  common.legend = TRUE, legend = "right", ncol = 1, align = "v")

