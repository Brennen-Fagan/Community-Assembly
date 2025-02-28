# Try to do 8g, but in place.

# datfolders <- dir(pattern = "TSTS_Simulations_")
# datfolders <- dir(pattern = "TSTS_Simulations_.+2024-11-30$") # Regex
datfolders <- dir(path = ".",
                  pattern = "TSTS_Simulations_.+2025-01-2[1-4]$",# Regex
                  full.names = TRUE)

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

# functions: ##################################################################

# Not a finished function!
interventionNamingScheme <- function(aff, ppa, ipt) {
  aDO <- affinityDictionaryOrigin[aff, ]
  ppDO <- poolpatchDictionaryOrigin[ppa, ]

  if (explicit <- grepl(pattern = "rep", aDO$PatchAffinities)) {
    initState <-
      paste0("(", paste( # NOT PRETTY FOR 10, MAY WANT TO JUST REPORT FUNC CALL
        vals <- retrieveFunction(aDO$PatchAffinities)(ppDO$NumberEnvironments),
        collapse = ", "), ")")
  } else {
    initState <-
      paste0(aDO$PatchAffinities, "(", ppDO$NumberEnvironments, ")")
  }

  if(is.na(ipt)) {return(initState)}

  ipDO <- interventionPatchDictionaryOrigin[ipt, ]

  if (ppDO$NumberEnvironments == 1) {
    finState <- paste0(
      "(", retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments), ")"
    )
  } else if (is.na(ipDO$InterventionLocation) ||
             !explicit ||
             !grepl(pattern = "rep", ipDO$PatchAffinities)) {

    finState <- # InterventionPercentage is a bit of a misnomer!
      paste0(ipDO$InterventionPercentage * 100, "%", ipDO$PatchAffinities)

  } else if (ipDO$InterventionLocation == 0) {# Left

    valsnew <- retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments)

    finState <- paste0("(", paste(
      valsnew[1:(ppDO$NumberEnvironments*ipDO$InterventionPercentage)],
      vals[
        (ppDO$NumberEnvironments*ipDO$InterventionPercentage + 1):
          ppDO$NumberEnvironments],
      collapse = ", ", sep = ", "), ")")

  } else if (ipDO$InterventionLocation == 1) {# Right

    valsnew <- retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments)

    finState <- paste0("(", paste(
      vals[1:(ppDO$NumberEnvironments*(1 - ipDO$InterventionPercentage))],
      valsnew[
        (ppDO$NumberEnvironments*(1 - ipDO$InterventionPercentage) + 1):
          ppDO$NumberEnvironments],
      collapse = ", ", sep = ", "), ")")
  }

  return(paste0(initState, "->", finState))
}

flattenDiversity <- function(d) {
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

  pres <- tidytable::data.table(d$Presence) %>% tidytable::mutate(
    Environment1 = Environment,
    Environment2 = NA
  ) %>% tidytable::group_by(
    Time, Environment1, Environment2
  ) %>% tidytable::mutate(
    `Average Size:0` = mean(Size),
    `Average Size:1` = sum(Size*Abundance)/sum(Abundance),
    `St.Dev. Size:0` = sqrt(var(Size)),
    `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance))),
    `Ratio Con/Bas:0` = sum(Type == "Consumer")/sum(Type == "Basal"),
    `Ratio Con/Bas:1` = sum((Type == "Consumer") * Abundance) /
      sum((Type == "Basal") * Abundance),
    `Average Aff.:0` = mean(Affinity),
    `Average Aff.:1` = sum(Affinity*Abundance)/sum(Abundance)
  ) %>% tidytable::pivot_longer(
    cols = `Average Size:0`:`Average Aff.:1`,
    names_to = "Metric", values_to = "Value"
  ) %>% tidytable::mutate(
    Subset = NA
  )

  pres_subset <- tidytable::data.table(d$Presence) %>% tidytable::mutate(
    Environment1 = Environment,
    Environment2 = NA,
    Subset = paste0(Type, "_", Affinity)
  ) %>% tidytable::group_by(
    Time, Environment1, Environment2, Subset
  ) %>% tidytable::mutate(
    `Average Size:0` = mean(Size),
    `Average Size:1` = sum(Size*Abundance)/sum(Abundance),
    `St.Dev. Size:0` = sqrt(var(Size)),
    `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance)))
  ) %>% tidytable::pivot_longer(
    cols = `Average Size:0`:`St.Dev. Size:1`,
    names_to = "Metric", values_to = "Value"
  )

  tidytable::data.table(d$Diversity) %>% tidytable::bind_rows(
    pres, pres_subset
  ) %>% tidytable::mutate(
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
    Affinity = id[[1]][7],
    AffinitySeed = id[[2]][5],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
  )
}

# Load Data: ##################################################################

for (datfolder in datfolders) {
  datfolderID <- paste0(
    strsplit(datfolder, split = "_")[[1]][-c(1:2)],
    collapse = "_")

  filestring <- paste0("diversitiesFlattened9_",datfolderID,".RData")

  if (file.exists(filestring)) {next()}

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

save(diversitiesAll, file = "diversitiesFlattened9a9_subset2.RData")

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
  baseplot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Time, y = Value,
      group = interaction(
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
        Time,
        Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
      ) %>% tidytable::summarise(
        top = quantile(Value, probs = CI+(1-CI)/2, na.rm = TRUE),
        bot = quantile(Value, probs = (1-CI)-(1-CI)/2, na.rm = TRUE)
      ), mapping = ggplot2::aes(
        x = Time, ymin = bot, ymax = top,
        group = interaction(
          Intervention, InterventionInitial, InterventionFinal, SpeciesAffinity
          ),
        fill = Intervention,
        color = Intervention
      ), inherit.aes = FALSE,
      alpha = 0.25, linewidth = 0.25 #, linetype = "dotted"
    )
  }
  baseplot <- baseplot + ggplot2::geom_smooth(
    # color = "black"
  ) + ggplot2::facet_grid(
    facets
  ) + ggplot2::scale_color_manual(
    values = colorPalette, aesthetics = c("color", "fill"),
    name = "Island Land-use"
  ) + ggplot2::scale_linetype(
    name = "Species Preferences"
  )
  return(baseplot)
}

plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0"
))

lapply(unique(diversitiesAll$Metric), function(metric) {
  thePlot <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
    is.na(Subset), Metric == metric
  )) + ggplot2::labs(y = metric)
  if (! grepl(pattern = "Alpha", x = metric, fixed = TRUE)) {
    # The alphas routinely escape [0, 1], but the Betas, Ratios, and Sizes don't
    thePlot <- thePlot + ggplot2::coord_cartesian(ylim = c(0, 1))
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

plot1 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
  Intervention %in% c("(0)", "(0.5)", "(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(. ~ SpeciesAffinity), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are species preferences, inner 75% intervals"
)
ggplot2::ggsave(
  plot1,
  filename = "Figure2_Prototype1.png",
  units = "px", height = 1600, width = 3200
  )

plot2 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
  Intervention %in% c("(0)", "(0.5)", "(1)")#, "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~Intervention), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are patch type, inner 75% intervals"
)
ggplot2::ggsave(
  plot2,
  filename = "Figure2_Prototype2.png",
  units = "px", height = 1600, width = 3200
)

plot3 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
  Intervention %in% c("(0)", "(0.5)", "(0)->(0.5)", "(0.5)->(0)")
), facets = as.formula(.~InterventionInitial), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are starting patch type, inner 75% intervals"
)
ggplot2::ggsave(
  plot3,
  filename = "Figure3_Prototype1.png",
  units = "px", height = 1600, width = 3200
)

plot3s1 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)"),
  Intervention %in% c("(0)", "(1)", "(0)->(1)", "(1)->(0)")
), facets = as.formula(.~InterventionInitial), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are starting patch type, inner 75% intervals"
)
ggplot2::ggsave(
  plot3s1,
  filename = "Figure3s1_Prototype1.png",
  units = "px", height = 1600, width = 3200
)

plot3s2 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("50% 0, 50% 1")#,
  # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
), facets = as.formula(.~InterventionInitial), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are starting patch type, inner 75% intervals"
)
ggplot2::ggsave(
  plot3s2,
  filename = "Figure3s2_Prototype1.png",
  units = "px", height = 1600, width = 3200
)

plot3s4_1000 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("100% 0")#,
  # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
), facets = as.formula(.~InterventionFinal), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are ending patch type, inner 75% intervals"
)
ggplot2::ggsave(
  plot3s4_1000,
  filename = "Figure3s4_Prototype1_1000.png",
  units = "px", height = 1600, width = 3200
)

plot3s4_5050 <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("50% 0, 50% 1")#,
  # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
), facets = as.formula(.~InterventionFinal), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are ending patch type, inner 75% intervals"
)
ggplot2::ggsave(
  plot3s4_5050,
  filename = "Figure3s4_Prototype1_5050.png",
  units = "px", height = 1600, width = 3200
)

plot3s4_unif <- plotMeanAndInner(diversitiesAll %>% tidytable::filter(
  is.na(Subset), Metric == "Alpha Hill:0",
  SpeciesAffinity %in% c("Uniform(0, 1)")#,
  # Intervention %in% c("(0)", "(0.5)", "(1)", "(0)->(1)", "(1)->(0)")
), facets = as.formula(.~InterventionFinal), CIs = c(0.75)
) + ggplot2::theme_bw() + ggplot2::ylab("Richness") + ggplot2::labs(
  subtitle = "200 species, labels are ending patch type, inner 75% intervals"
)
ggplot2::ggsave(
  plot3s4_unif,
  filename = "Figure3s4_Prototype1_Uniform.png",
  units = "px", height = 1600, width = 3200
)

