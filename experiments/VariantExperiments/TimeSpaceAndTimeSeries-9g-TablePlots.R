# Try to do 8g, but in place.

# datfolders <- dir(pattern = "TSTS_Simulations_")
# datfolders <- dir(pattern = "TSTS_Simulations_.+2024-11-30$") # Regex
datfolders <- dir(path = ".",
                  pattern = "TSTS_Simulations_.+2025-01-2[1-4]$",# Regex
                  full.names = TRUE)

overwrite <- TRUE

# Problems with X11
options(bitmapType = "cairo")

# Libraries: ##################################################################
source("TimeSpaceAndTimeSeries-9-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')

library(ggplot2)
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
  "Uniform(0, 1)" = "dotdash"
)

# functions: ##################################################################

source("flattenDiversity.R")

# Load Data: ##################################################################

for (datfolder in datfolders) {
  datfolderID <- paste0(
    strsplit(datfolder, split = "_")[[1]][-c(1:2)],
    collapse = "_")

  filestring <- paste0("diversitiesFlattened9_",datfolderID,".RData")

  if (file.exists(filestring) && !overwrite) {next()}

  diversities <- lapply(
    dir(datfolder, full.names = TRUE, pattern = "Diversity"), function(x) {
      names <- load(x)
      stopifnot(length(names) == 1)
      obj <- get(names)
      # obj <- tidytable::data.table(obj) # Cannot Table yet.
      return(c(obj, "Dir" = dirname(x), "File" = basename(x)))
    })

  diversitiesFlattened <- vector(mode = "list", length = length(diversities))
  for(i in 1:length(diversitiesFlattened)) {
    # Pop front of diversities, process, put in flattened, remove.
    # Hence, len(diversities) changes, but pre-alloc => len(flattened) is not.
    diversitiesFlattened[[i]] <- flattenDiversity(diversities[[1]])
    diversities[[1]] <- NULL
    gc()
  }

  rm(diversities)

  diversitiesFlattened <- do.call(rbind, diversitiesFlattened)

  # Human readable patch affinities
  diversitiesInterventionStrings <- diversitiesFlattened %>% dplyr::select(
    Affinity, PoolPatch, InterventionPatchType
  ) %>% dplyr::distinct(
  ) %>% dplyr::mutate(
    Intervention = unlist(mapply(
      FUN = interventionNamingScheme,
      Affinity, PoolPatch, InterventionPatchType
    ))
  )

  # Col-wise append
  diversitiesFlattened <- diversitiesFlattened %>% dplyr::left_join(
    diversitiesInterventionStrings,
    by = c("Affinity", "PoolPatch", "InterventionPatchType"),
    multiple = "all"
  )

  # Human readable species affinities
  diversitiesFlattened <- diversitiesFlattened %>% dplyr::mutate(
    SpeciesAffinity =
      affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
  )

  # Correct the NA for richness values
  diversitiesFlattened <- diversitiesFlattened %>% dplyr::mutate(
    Value = dplyr::case_when(
      Metric == "Alpha Hill:0" & is.na(Value) ~ 0,
      TRUE ~ Value
    )
  )

  save(diversitiesFlattened,
       file = filestring)

  rm(diversitiesFlattened)
  gc()
}

diversitiesAll <- NULL
for (datfolder in datfolders) {
  datfolderID <- paste0(
    strsplit(datfolder, split = "_")[[1]][-c(1:2)],
    collapse = "_")

  filestring <- paste0("diversitiesFlattened9_", datfolderID,".RData")

  stopifnot(file.exists(filestring))

  objnames <- load(filestring)
  obj <- get(objnames)

  # We just can't get past the memory barrier, so we'll need to reduce the
  # amount of data we are looking at.
  obj <- obj %>% tidytable::filter(
    (round(Time, digits = -1) %% 100) == 0,
    Intervention %in% c(
      "(0)",   "(0)->(0.5)", "(0)->(1)",
      "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
      "(1)",   "(1)->(0)",   "(1)->(0.5)")
    )

  if(!is.null(diversitiesAll)) {
    diversitiesAll <- tidytable::bind_rows(diversitiesAll, obj)
  } else {
    diversitiesAll <- obj
  }

  rm(obj)
  rm(list = objnames)
  gc()
}

save(diversitiesAll, file = "diversitiesFlattened9a9_subset4.RData")

# START HERE: ##################################################################

diversitiesAll <- diversitiesAll %>% tidytable::mutate(
  SpeciesAffinity = tidytable::case_when(
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

plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  PoolPatchSeed %in% as.character(343:386)
))

lapply(unique(diversitiesAll$Metric), function(metric) {
  print(metric)
  thePlot <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
    is.na(Subset), Metric == metric, PoolPatchSeed %in% as.character(343:386)
  )) + ggplot2::labs(y = metric)
  if (! grepl(pattern = "Alpha", x = metric, fixed = TRUE) &&
      ! grepl(pattern = "Size", x = metric, fixed = TRUE) &&
      ! grepl(pattern = "Ratio", x = metric, fixed = TRUE)) {
    # The alphas routinely escape [0, 1], but the Betas and Aff. don't
    thePlot <- thePlot + ggplot2::coord_cartesian(ylim = c(0, 1))
  }
  if (grepl(pattern = "Ratio", x = metric, fixed = TRUE) ||
      grepl(pattern = "Average Size", x = metric, fixed = TRUE) # ||
      # grepl(pattern = "Average LSize", x = metric, fixed = TRUE)
      ) {
    thePlot <- thePlot + ggplot2::scale_y_log10()
  }
  ggplot2::ggsave(
    thePlot,
    filename = paste0("FigureMetric_",
                      gsub(metric, pattern = "[ .:/]",
                           replacement = ""),
                      ".png"),
    units = "px", height = 2400, width = 3200
  )
})

# plot1 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
#   Intervention %in% c("(0)", "(0.5)", "(1)")#, "(0)->(0.5)", "(0.5)->(0)")
# ), facets = as.formula(. ~ SpeciesAffinity), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are species preferences, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot1,
#   filename = "Figure2_Prototype1.png",
#   units = "px", height = 1600, width = 3200
#   )
#
# plot2 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
#   Intervention %in% c("(0)", "(0.5)", "(1)")#, "(0)->(0.5)", "(0.5)->(0)")
# ), facets = as.formula(.~Intervention), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are patch type, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot2,
#   filename = "Figure2_Prototype2.png",
#   units = "px", height = 1600, width = 3200
# )

plot1 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("100% 0"),
  Intervention %in% c("(0)", "(0.5)", "(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~.), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, inner 75% intervals"
)
ggplot2::ggsave(
  plot1,
  filename = "Figure2_Prototype3.png",
  units = "px", height = 1600, width = 2000
)

plot1s1 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Ratio Con/Bas:0",
  SpeciesAffinity %in% c("100% 0"),
  Intervention %in% c("(0)", "(0.5)", "(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~.), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab(
  "No. of Consumer Species / No. of Basal Species"
  ) + ggplot2::labs(
  subtitle = "200 species, inner 75% intervals"
)
plot1s2 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Average Size:0",
  SpeciesAffinity %in% c("100% 0"),
  Intervention %in% c("(0)", "(0.5)", "(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~.), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab(
  "Average size of species present"
) + ggplot2::labs(
  subtitle = "200 species, inner 75% intervals"
)
ggplot2::ggsave(
  ggpubr::ggarrange(plot1s1, plot1s2,
                    common.legend = TRUE, legend = "right", ncol = 2),
  filename = "Figure2s1_Prototype1.png",
  units = "px", height = 1600, width = 3200
)

plot2 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("100% 0"),
  Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~.), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, inner 75% intervals"
)
plot2s1 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Ratio Con/Bas:0",
  SpeciesAffinity %in% c("100% 0"),
  Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~.), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab(
  "No. of Consumer Species / No. of Basal Species"
) + ggplot2::labs(
  subtitle = "200 species, inner 75% intervals"
)
plot2s2 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Average Size:0",
  SpeciesAffinity %in% c("100% 0"),
  Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~.), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab(
  "Average size of species present"
) + ggplot2::labs(
  subtitle = "200 species, inner 75% intervals"
)
ggplot2::ggsave(
  ggpubr::ggarrange(plot2, plot2s1, plot2s2,
                    common.legend = TRUE, legend = "right", ncol = 3),
  filename = "Figure3_Prototype2.png",
  units = "px", height = 1600, width = 3200
)


# plot3 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
#   Intervention %in% c("(0)", "(0.5)", "(0)->(0.5)", "(0.5)->(0)")
# ), facets = as.formula(.~InterventionInitial), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are starting patch type, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot3,
#   filename = "Figure3_Prototype1.png",
#   units = "px", height = 1600, width = 3200
# )
#
# plot3s1 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
#   Intervention %in% c("(0)", "(1)", "(0)->(1)", "(1)->(0)")
# ), facets = as.formula(.~InterventionInitial), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are starting patch type, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot3s1,
#   filename = "Figure3s1_Prototype1.png",
#   units = "px", height = 1600, width = 3200
# )
#
# plot3s2 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("50% 0, 50% 1")#,
#   # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
# ), facets = as.formula(.~InterventionInitial), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are starting patch type, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot3s2,
#   filename = "Figure3s2_Prototype1.png",
#   units = "px", height = 1600, width = 3200
# )
#
# plot3s4_1000 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("100% 0")#,
#   # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
# ), facets = as.formula(.~InterventionFinal), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are ending patch type, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot3s4_1000,
#   filename = "Figure3s4_Prototype1_1000.png",
#   units = "px", height = 1600, width = 3200
# )
#
# plot3s4_5050 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("50% 0, 50% 1")#,
#   # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
# ), facets = as.formula(.~InterventionFinal), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are ending patch type, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot3s4_5050,
#   filename = "Figure3s4_Prototype1_5050.png",
#   units = "px", height = 1600, width = 3200
# )
#
# plot3s4_unif <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
#   is.na(Subset), Metric == "Alpha Hill:0",
#   SpeciesAffinity %in% c("Uniform(0, 1)")#,
#   # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
# ), facets = as.formula(.~InterventionFinal), CIs = c(0.75)
# ) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
#   subtitle = "200 species, labels are ending patch type, inner 75% intervals"
# )
# ggplot2::ggsave(
#   plot3s4_unif,
#   filename = "Figure3s4_Prototype1_Uniform.png",
#   units = "px", height = 1600, width = 3200
# )

scatterBase <- ggplot2::ggplot(
  diversitiesAll %>% dplyr::filter(
    is.na(Subset),
    Metric %in% c("Alpha Hill:0", "Average Aff.:0", "Average Size:0",
                  "Ratio Con/Bas:0", "St.Dev. Size:0", "TimeJaccard")
  ) %>% dplyr::distinct(
  ) %>% dplyr::group_by(
    Environment1, Environment2, Metric, Subset, PoolPatch, PoolPatchSeed,
    Interactions, InteractionsSeed, Events, EventsSeed, InitialConditions,
    InitialConditionsSeed, Dispersal, NicheDistance, AffinitySeed,
    InterventionPatchType, InterventionTimeType, InterventionTimeSeed,
    InterventionDispersal, InterventionNicheDistance,
    Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
  ) %>% dplyr::summarise(
    Value = mean(Value),
    .groups = "drop"
  ) %>% tidyr::pivot_wider(names_from = Metric, values_from = Value),
  ggplot2::aes(
    color = Intervention, fill = Intervention, shape = SpeciesAffinity
  )
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Island Land-use"
)

# "Alpha Hill:0", "Average Aff.:0", "Average Size:0",
# "Ratio Con/Bas:0", "St.Dev. Size:0", "TimeJaccard"
scatterBase + ggplot2::geom_point(
  ggplot2::aes(x = `Alpha Hill:0`, y = `TimeJaccard`)
) + ggplot2::facet_grid(Intervention ~ SpeciesAffinity)

