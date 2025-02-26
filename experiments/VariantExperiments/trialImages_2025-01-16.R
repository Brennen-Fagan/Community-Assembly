# requires presences to compute size/trait distributions through time.


# Problems with X11
options(bitmapType = "cairo")

# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)

source("TimeSpaceAndTimeSeries-9-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')

target <- "diversitiesFlattened9a8_plottable.RData"
# or
# target <- "diversitiesFlattened9a8-logistic_plottable.RData"
load(target)

if (target == "diversitiesFlattened9a8_plottable.RData") {
  date <- "2024-11-30"
  path <- file.path(".", "Deprecated")
  label <- "Base"
} else {
  date <- "2025-01-26" # and 27
  path <- file.path(".")
  label <- "Logistic"
}

diversitiesFlattenedSubset <- diversitiesFlattened %>% dplyr::group_by(
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention
) %>% dplyr::filter(
  Metric == "Alpha Hill:0", is.na(Subset), NicheDistance == 7,
  Intervention %in% c(
    "(0)", "(0.5)", "(1)", "(0)->(0.5)", "(0)->(1)",
    "(0.5)->(0)", "(0.5)->(1)", "(1)->(0)", "(1)->(0.5)"
  ),
  SpeciesAffinity %in% c("rep_0", "evensplit_01", "runif")
)
# Leaves pool variation (x4),
# the independent runs (x2),
# and then the base cases and interventions (x3 and x3) = 72 runs

# We have the presence data, but there's too much to load and analyse all at
# once, so we need to be more selective.
tempIDs <- dplyr::bind_rows(
  diversitiesFlattenedSubset %>% dplyr::ungroup() %>% dplyr::select(
    # PoolPatch:AffinitySeed # Affinity can get out of order...?
    PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
    InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
    Affinity, AffinitySeed
  ) %>% dplyr::distinct() %>% tidyr::unite(
    col = "Col2", ends_with("Seed"), sep = "-"
  ) %>% tidyr::unite(
    col = "Col1", !starts_with("Col"), sep = "-"
  ) %>% tidyr::unite(col = "ID"),
  diversitiesFlattenedSubset %>% dplyr::ungroup() %>% dplyr::select(
    PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events, EventsSeed,
    InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
    Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
    InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
    InterventionNicheDistance
  ) %>% dplyr::distinct() %>% tidyr::unite(
    col = "Col4", ends_with("Seed") & starts_with("Intervention"), sep = "-"
  ) %>% tidyr::unite(
    col = "Col3", starts_with("Intervention"), sep = "-"
  ) %>% tidyr::unite(
    col = "Col2", ends_with("Seed"), sep = "-"
  ) %>% tidyr::unite(
    col = "Col1", !starts_with("Col"), sep = "-"
  ) %>% tidyr::unite(col = "ID")
) %>% dplyr::mutate(ID = paste0("TSTS_Diversity_", ID, ".RData"))

tempDirs <- diversitiesFlattenedSubset %>% dplyr::ungroup() %>% dplyr::select(
  PoolPatch, Interactions, PoolPatchSeed, InteractionsSeed
) %>% dplyr::distinct() %>% tidyr::unite(
  col = "Col2", ends_with("Seed"), sep = "-"
) %>% tidyr::unite(
  col = "Col1", !starts_with("Col"), sep = "-"
) %>% tidyr::unite(
  col = "ID"
  ) %>% dplyr::mutate(
    ID = file.path(path,
                   paste0("TSTS_Simulations_", ID, "_", "2024-11-30"#"2025-01-27"
                          ))
  )

tempLoad <- dir(path = tempDirs$ID, full.names = TRUE)
tempLoad <- tempLoad[basename(tempLoad) %in% tempIDs$ID]
stopifnot(length(tempLoad)%%27 == 0,
          length(tempLoad) > 0) # was length == 54 for two sets.

diversities <- lapply(
  tempLoad, function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = dirname(x), "File" = basename(x)))
  })

presencesFlattened <- do.call(rbind, lapply(diversities, function(d) {
  if ("FullID" %in% names(d$Ellipsis)) {
    id <- strsplit(
      strsplit(d$Ellipsis$FullID, "_", fixed = TRUE)[[1]], # Split seeds off.
      "-", fixed = TRUE)
  } else if ("ID" %in% names(d$Ellipsis)) {
    id <- strsplit(
      strsplit(d$Ellipsis$ID, "_", fixed = TRUE)[[1]], # Split seeds off.
      "-", fixed = TRUE)
  } else {
    id <- strsplit(
      strsplit(
        strsplit(d$File, ".", fixed = TRUE)[[1]][1], # Remove .RData.
        "_", fixed = TRUE)[[1]][3:4], # Remove TSTS_Type and split seeds off.
      "-", fixed = TRUE # Separate out the id values.
    )
  }

  if (length(id) < 3) {
    # I.e., no intervention.
    id[[3]] <- rep(NA, 4)
    id[[4]] <- rep(NA, 2)
  }

  # pres <- d$Presence %>% dplyr::mutate(
  #   Environment1 = Environment,
  #   Environment2 = NA
  # ) %>% dplyr::group_by(
  #   Time, Environment1, Environment2
  # ) %>% dplyr::mutate(
  #   `Average Size:0` = mean(Size),
  #   `Average Size:1` = mean(Size*Abundance)/sum(Abundance),
  #   `St.Dev. Size:0` = sqrt(var(Size)),
  #   `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance))),
  #   `Ratio Con/Bas:0` = sum(Type == "Consumer")/sum(Type == "Basal"),
  #   `Ratio Con/Bas:1` = sum((Type == "Consumer") * Abundance) /
  #     sum((Type == "Basal") * Abundance),
  #   `Average Aff.:0` = mean(Affinity),
  #   `Average Aff.:1` = mean(Affinity*Abundance)/sum(Abundance)
  # ) %>% tidyr::pivot_longer(
  #   cols = `Average Size:0`:`Average Aff.:1`,
  #   names_to = "Metric", values_to = "Value"
  # ) %>% dplyr::mutate(
  #   Subset = NA
  # )
  #
  # pres_subset <- d$Presence %>% dplyr::mutate(
  #   Environment1 = Environment,
  #   Environment2 = NA,
  #   Subset = paste0(Type, "_", Affinity)
  # ) %>% dplyr::group_by(
  #   Time, Environment1, Environment2, Subset
  # ) %>% dplyr::mutate(
  #   `Average Size:0` = mean(Size),
  #   `Average Size:1` = mean(Size*Abundance)/sum(Abundance),
  #   `St.Dev. Size:0` = sqrt(var(Size)),
  #   `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance)))
  # ) %>% tidyr::pivot_longer(
  #   cols = `Average Size:0`:`St.Dev. Size:1`,
  #   names_to = "Metric", values_to = "Value"
  # )

  d$Presence %>%
    # dplyr::bind_rows(pres, pres_subset) %>%
    dplyr::mutate(
    PoolPatch = id[[1]][1],
    PoolPatchSeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    InitialConditions = id[[1]][4],
    InitialConditionsSeed = id[[2]][4],
    Dispersal = id[[1]][5],
    NicheDistance = id[[1]][6],
    AffinitySetting = id[[1]][7],
    AffinitySeed = id[[2]][5],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
  )
}))

rm(diversities)

diversitiesInterventionStrings <- presencesFlattened %>% dplyr::select(
  AffinitySetting, PoolPatch, InterventionPatchType
) %>% dplyr::distinct(
) %>% dplyr::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    AffinitySetting, PoolPatch, InterventionPatchType
  ))
)

presencesFlattened <- presencesFlattened %>% dplyr::left_join(
  diversitiesInterventionStrings,
  by = c("AffinitySetting", "PoolPatch", "InterventionPatchType")
)

presencesFlattened <- presencesFlattened %>% dplyr::mutate(
  SpeciesAffinity =
    affinityDictionaryOrigin$SpeciesAffinities[as.numeric(AffinitySetting)]
)

# colors:
colorPalette <- c(
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

# plots:
# ImageRichness <- ggplot2::ggplot(
#   diversitiesFlattenedSubset %>% dplyr::ungroup(
#   ) %>% dplyr::mutate(
#     InterventionPatchType = ifelse(is.na(InterventionPatchType),
#                                    "0", InterventionPatchType),
#   #   TimeSubset = dplyr::case_when(
#   #     Time < 5000 ~ "1: Early",
#   #     Time > 30000 ~ "3: Late",
#   #     Time > 15000 & Time < 22000 ~ "2: Middle",
#   #     TRUE ~ NA_character_
#   #   )
#   Intervention = factor(Intervention, levels = c(
#     "(0)", "(0)->(0.5)", "(0)->(1)",
#     "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
#     "(1)", "(1)->(0)", "(1)->(0.5)"
#   ), ordered = TRUE)
#   ) %>% dplyr::filter(
#   #   !is.na(TimeSubset)
#     PoolPatchSeed == PoolPatchSeed[1]
#   ),
#   ggplot2::aes(x = Time, y = Value, color = Intervention,
#                linetype = InterventionPatchType)
# ) + ggplot2::geom_line(
# ) + ggplot2::facet_grid(
#   # SpeciesAffinity ~ .,# PoolPatchSeed, # + TimeSubset,
#   . ~ SpeciesAffinity,
#   scales = "free"
# ) + ggplot2::labs(
#   linetype = "", y = "Richness"
# ) + ggplot2::scale_color_manual(
#   values = colorPalette
# ) + ggplot2::guides(
#   linetype = "none"
# )
#
# ImageAbundance <- ggplot2::ggplot(
#   presencesFlattened %>% dplyr::ungroup(
#   ) %>% dplyr::mutate(
#     InterventionPatchType = ifelse(is.na(InterventionPatchType),
#                                    "0", InterventionPatchType),
#     TimeGroup = dplyr::case_when(
#       round(Time, -1) == 100 ~ 1,
#       round(Time, -1) == 10000 ~ 2,
#       round(Time, -1) == 20000 ~ 3,
#       round(Time, -1) == 30000 ~ 4,
#       TRUE ~ NA_real_
#     ),
#     Intervention = factor(Intervention, levels = c(
#       "(0)", "(0)->(0.5)", "(0)->(1)",
#       "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
#       "(1)", "(1)->(0)", "(1)->(0.5)"
#     ), ordered = TRUE)
#   ) %>% dplyr::filter(!is.na(TimeGroup), PoolPatchSeed == PoolPatchSeed[1]),
#   ggplot2::aes(x = Time, y = Abundance, fill = Intervention,
#                # linetype = InterventionPatchType,
#                group = interaction(Time, Intervention))
# ) + ggplot2::geom_violin(
#   adjust = 1/4, trim = TRUE, position = ggplot2::position_dodge(),
#   width = 8000, scale = "width", alpha = 0.5
# ) + ggplot2::facet_grid(
#   # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
#   . ~ SpeciesAffinity,
#   scales = "free"
# ) + ggplot2::labs(
#   linetype = ""
# ) + ggplot2::scale_color_manual(
#   values = colorPalette, aesthetics = c("color", "fill")
# ) + ggplot2::guides(
#   linetype = "none"
# ) + ggplot2::scale_y_log10(
# )
#
# ImageSize <- ggplot2::ggplot(
#   presencesFlattened %>% dplyr::ungroup(
#   ) %>% dplyr::mutate(
#     InterventionPatchType = ifelse(is.na(InterventionPatchType),
#                                    "0", InterventionPatchType),
#     TimeGroup = dplyr::case_when(
#       round(Time, -1) == 100 ~ 1,
#       round(Time, -1) == 10000 ~ 2,
#       round(Time, -1) == 20000 ~ 3,
#       round(Time, -1) == 30000 ~ 4,
#       TRUE ~ NA_real_
#     ),
#     Intervention = factor(Intervention, levels = c(
#       "(0)", "(0)->(0.5)", "(0)->(1)",
#       "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
#       "(1)", "(1)->(0)", "(1)->(0.5)"
#     ), ordered = TRUE)
#   ) %>% dplyr::filter(!is.na(TimeGroup), PoolPatchSeed == PoolPatchSeed[1]),
#   ggplot2::aes(x = Time, y = Size, fill = Intervention,
#                # linetype = InterventionPatchType,
#                group = interaction(Time, Intervention))
# ) + ggplot2::geom_violin(
#   adjust = 1/4, trim = TRUE, position = ggplot2::position_dodge(),
#   width = 8000, scale = "width", alpha = 0.5
# ) + ggplot2::facet_grid(
#   # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
#   . ~ SpeciesAffinity,
#   scales = "free"
# ) + ggplot2::labs(
#   linetype = ""
# ) + ggplot2::scale_color_manual(
#   values = colorPalette, aesthetics = c("color", "fill")
# ) + ggplot2::guides(
#   linetype = "none"
# ) + ggplot2::scale_y_log10(
# )
#
# ImageAffinity <- ggplot2::ggplot(
#   presencesFlattened %>% dplyr::ungroup(
#   ) %>% dplyr::mutate(
#     InterventionPatchType = ifelse(is.na(InterventionPatchType),
#                                    "0", InterventionPatchType),
#     TimeGroup = dplyr::case_when(
#       round(Time, -1) == 100 ~ 1,
#       round(Time, -1) == 10000 ~ 2,
#       round(Time, -1) == 20000 ~ 3,
#       round(Time, -1) == 30000 ~ 4,
#       TRUE ~ NA_real_
#     ),
#     Intervention = factor(Intervention, levels = c(
#       "(0)", "(0)->(0.5)", "(0)->(1)",
#       "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
#       "(1)", "(1)->(0)", "(1)->(0.5)"
#     ), ordered = TRUE)
#   ) %>% dplyr::filter(!is.na(TimeGroup), PoolPatchSeed == PoolPatchSeed[1]),
#   ggplot2::aes(x = Time, y = Affinity, fill = Intervention,
#                # linetype = InterventionPatchType,
#                group = interaction(Time, Intervention))
# ) + ggplot2::geom_violin(
#   adjust = 1/4, trim = TRUE, position = ggplot2::position_dodge(),
#   width = 8000, scale = "width", alpha = 0.5
# ) + ggplot2::facet_grid(
#   # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
#   . ~ SpeciesAffinity,
#   scales = "free"
# ) + ggplot2::labs(
#   linetype = ""
# ) + ggplot2::scale_color_manual(
#   values = colorPalette, aesthetics = c("color", "fill")
# ) + ggplot2::guides(
#   linetype = "none"
# # ) + ggplot2::scale_y_log10(
# )
#
# # ggpubr::ggarrange(ImageRichness, ImageAbundance, ImageSize, ImageAffinity,
# #                   common.legend = TRUE, legend = "bottom", nrow = 1)
# ggpubr::ggarrange(
#   ImageRichness + ggplot2::theme(axis.title.x = ggplot2::element_blank()),
#   ImageAbundance + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
#                                   strip.background = ggplot2::element_blank(),
#                                   strip.text = ggplot2::element_blank()),
#   ImageSize + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
#                              strip.background = ggplot2::element_blank(),
#                              strip.text = ggplot2::element_blank()),
#   ImageAffinity + ggplot2::theme(strip.background = ggplot2::element_blank(),
#                                  strip.text = ggplot2::element_blank()),
#                   common.legend = TRUE, legend = "right", ncol = 1, align = "v")

ImageRichness <- ggplot2::ggplot(
  diversitiesFlattenedSubset %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    InterventionPatchType = ifelse(is.na(InterventionPatchType),
                                   "0", InterventionPatchType),
    #   TimeSubset = dplyr::case_when(
    #     Time < 5000 ~ "1: Early",
    #     Time > 30000 ~ "3: Late",
    #     Time > 15000 & Time < 22000 ~ "2: Middle",
    #     TRUE ~ NA_character_
    #   )
    Intervention = factor(Intervention, levels = c(
      "(0)", "(0)->(0.5)", "(0)->(1)",
      "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
      "(1)", "(1)->(0)", "(1)->(0.5)"
    ), ordered = TRUE)
  ) %>% dplyr::filter(
    #   !is.na(TimeSubset)
    PoolPatchSeed == unique(PoolPatchSeed)[2]
  ),
  ggplot2::aes(x = Time, y = Value, color = Intervention,
               linetype = InterventionPatchType)
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # SpeciesAffinity ~ .,# PoolPatchSeed, # + TimeSubset,
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::labs(
  linetype = "", y = "Richness"
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::guides(
  linetype = "none"
) + ggplot2::theme_minimal(
)

ImageAbundance <- ggplot2::ggplot(
  presencesFlattened %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    InterventionPatchType = ifelse(is.na(InterventionPatchType),
                                   "0", InterventionPatchType),
    TimeGroup = dplyr::case_when(
      # round(Time, -1) == 100 ~ 1,
      round(Time, -1) == 10000 ~ 2,
      round(Time, -1) == 20000 ~ 3,
      round(Time, -1) == 30000 ~ 4,
      TRUE ~ NA_real_
    ),
    Intervention = factor(Intervention, levels = c(
      "(0)", "(0)->(0.5)", "(0)->(1)",
      "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
      "(1)", "(1)->(0)", "(1)->(0.5)"
    ), ordered = TRUE)
  ) %>% dplyr::filter(
    !is.na(TimeGroup), PoolPatchSeed == unique(PoolPatchSeed)[2]
    ),
  ggplot2::aes(x = Time, y = Abundance, color = Intervention,
               # linetype = InterventionPatchType,
               group = interaction(Intervention))
) + ggplot2::geom_point(
  position = ggplot2::position_dodge(width = 6000), alpha = 0.5, size = 1
) + ggplot2::facet_grid(
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::labs(
  linetype = ""
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::guides(
  linetype = "none"
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

ImageSize <- ggplot2::ggplot(
  presencesFlattened %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    InterventionPatchType = ifelse(is.na(InterventionPatchType),
                                   "0", InterventionPatchType),
    TimeGroup = dplyr::case_when(
      # round(Time, -1) == 100 ~ 1,
      round(Time, -1) == 10000 ~ 2,
      round(Time, -1) == 20000 ~ 3,
      round(Time, -1) == 30000 ~ 4,
      TRUE ~ NA_real_
    ),
    Intervention = factor(Intervention, levels = c(
      "(0)", "(0)->(0.5)", "(0)->(1)",
      "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
      "(1)", "(1)->(0)", "(1)->(0.5)"
    ), ordered = TRUE)
  ) %>% dplyr::filter(
    !is.na(TimeGroup), PoolPatchSeed == unique(PoolPatchSeed)[2]
    ),
  ggplot2::aes(x = Time, y = Size, color = Intervention,
               # linetype = InterventionPatchType,
               group = interaction(Intervention))
) + ggplot2::geom_point(
  position = ggplot2::position_dodge(width = 6000), alpha = 0.5, size = 1
) + ggplot2::facet_grid(
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::labs(
  linetype = ""
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::guides(
  linetype = "none"
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

ImageAffinity <- ggplot2::ggplot(
  presencesFlattened %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    InterventionPatchType = ifelse(is.na(InterventionPatchType),
                                   "0", InterventionPatchType),
    TimeGroup = dplyr::case_when(
      # round(Time, -1) == 100 ~ 1,
      round(Time, -1) == 10000 ~ 2,
      round(Time, -1) == 20000 ~ 3,
      round(Time, -1) == 30000 ~ 4,
      TRUE ~ NA_real_
    ),
    Intervention = factor(Intervention, levels = c(
      "(0)", "(0)->(0.5)", "(0)->(1)",
      "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
      "(1)", "(1)->(0)", "(1)->(0.5)"
    ), ordered = TRUE)
  ) %>% dplyr::filter(
    !is.na(TimeGroup), PoolPatchSeed == unique(PoolPatchSeed)[2]
    ),
  ggplot2::aes(x = Time, y = Affinity, color = Intervention,
               # linetype = InterventionPatchType,
               group = interaction(Intervention))
) + ggplot2::geom_count(
  position = ggplot2::position_dodge(width = 6000), alpha = 0.6
) + ggplot2::facet_grid(
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::labs(
  linetype = ""
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::guides(
  linetype = "none"
  # ) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_continuous(
  breaks = c(10000, 20000, 30000),
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)

# ggpubr::ggarrange(ImageRichness, ImageAbundance, ImageSize, ImageAffinity,
#                   common.legend = TRUE, legend = "bottom", nrow = 1)
ImageComposite <- ggpubr::ggarrange(
  ImageRichness + ggplot2::theme(axis.title.x = ggplot2::element_blank()),
  ImageAbundance + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                  strip.background = ggplot2::element_blank(),
                                  strip.text = ggplot2::element_blank()),
  ImageSize + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                             strip.background = ggplot2::element_blank(),
                             strip.text = ggplot2::element_blank()),
  ImageAffinity + ggplot2::theme(strip.background = ggplot2::element_blank(),
                                 strip.text = ggplot2::element_blank()),
  common.legend = TRUE, legend = "right", ncol = 1, align = "v")

ggplot2::ggsave(
  ImageComposite, filename = paste0("Image_Composite342_", label, ".png"),
  width = 3200, height = 2400, units = "px"
)

# Biomasses?
ImageBiomass <- ggplot2::ggplot(
  presencesFlattened %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    InterventionPatchType = ifelse(is.na(InterventionPatchType),
                                   "0", InterventionPatchType),
    TimeGroup = dplyr::case_when(
      # round(Time, -1) == 100 ~ 1,
      round(Time, -1) == 10000 ~ 2,
      round(Time, -1) == 20000 ~ 3,
      round(Time, -1) == 30000 ~ 4,
      TRUE ~ NA_real_
    ),
    Intervention = factor(Intervention, levels = c(
      "(0)", "(0)->(0.5)", "(0)->(1)",
      "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
      "(1)", "(1)->(0)", "(1)->(0.5)"
    ), ordered = TRUE)
  ) %>% dplyr::filter(PoolPatchSeed == unique(PoolPatchSeed)[2]) %>% dplyr::group_by(
    Time, Intervention, SpeciesAffinity
  ) %>% dplyr::summarise(
    Biomass = sum(Size*Abundance)
  ),
  ggplot2::aes(x = Time, y = Biomass, color = Intervention,
               # linetype = InterventionPatchType,
               group = interaction(Intervention))
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  # SpeciesAffinity ~ .,#Intervention,# PoolPatchSeed, # + TimeSubset,
  . ~ SpeciesAffinity,
  scales = "free"
) + ggplot2::labs(
  linetype = ""
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill")
) + ggplot2::guides(
  linetype = "none"
) + ggplot2::scale_y_log10(
  minor_breaks = 10^(-4:4)
) + ggplot2::scale_x_continuous(
  minor_breaks = NULL
) + ggplot2::theme_minimal(
)
ggplot2::ggsave(
  ImageBiomass, filename = paste0("Image_Biomass342_", label, ".png"),
  width = 3200, height = 2400, units = "px"
  )
