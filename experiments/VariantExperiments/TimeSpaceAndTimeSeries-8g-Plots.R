# datfolders <- dir(pattern = "TSTS_Simulations_")
# datfolders <- dir(pattern = "TSTS_Simulations_.+2024-11-30$") # Regex
datfolders <- dir(path = ".",
                  pattern = "TSTS_Simulations_.+2025-01-2[67]$",# Regex
                  full.names = TRUE)

# Problems with X11
options(bitmapType = "cairo")

# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)

source("TimeSpaceAndTimeSeries-9-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')

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

# Load Data: ##################################################################
diversities <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Diversity"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = dirname(x), "File" = basename(x)))
  })

# Append Information of Parameters etc. to data.frames before consolidation.
diversitiesFlattened <- do.call(rbind, lapply(diversities, function(d) {
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

  pres <- d$Presence %>% dplyr::mutate(
    Environment1 = Environment,
    Environment2 = NA
  ) %>% dplyr::group_by(
    Time, Environment1, Environment2
  ) %>% dplyr::mutate(
    `Average Size:0` = mean(Size),
    `Average Size:1` = mean(Size*Abundance)/sum(Abundance),
    `St.Dev. Size:0` = sqrt(var(Size)),
    `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance))),
    `Ratio Con/Bas:0` = sum(Type == "Consumer")/sum(Type == "Basal"),
    `Ratio Con/Bas:1` = sum((Type == "Consumer") * Abundance) /
      sum((Type == "Basal") * Abundance),
    `Average Aff.:0` = mean(Affinity),
    `Average Aff.:1` = mean(Affinity*Abundance)/sum(Abundance)
  ) %>% tidyr::pivot_longer(
    cols = `Average Size:0`:`Average Aff.:1`,
    names_to = "Metric", values_to = "Value"
  ) %>% dplyr::mutate(
    Subset = NA
  )

  pres_subset <- d$Presence %>% dplyr::mutate(
    Environment1 = Environment,
    Environment2 = NA,
    Subset = paste0(Type, "_", Affinity)
  ) %>% dplyr::group_by(
    Time, Environment1, Environment2, Subset
  ) %>% dplyr::mutate(
    `Average Size:0` = mean(Size),
    `Average Size:1` = mean(Size*Abundance)/sum(Abundance),
    `St.Dev. Size:0` = sqrt(var(Size)),
    `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance)))
  ) %>% tidyr::pivot_longer(
    cols = `Average Size:0`:`St.Dev. Size:1`,
    names_to = "Metric", values_to = "Value"
  )

  d$Diversity %>% dplyr::bind_rows(pres, pres_subset) %>% dplyr::mutate(
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
}))

rm(diversities)

diversitiesInterventionStrings <- diversitiesFlattened %>% dplyr::select(
  Affinity, PoolPatch, InterventionPatchType
) %>% dplyr::distinct(
) %>% dplyr::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    Affinity, PoolPatch, InterventionPatchType
  ))
)

diversitiesFlattened <- diversitiesFlattened %>% dplyr::left_join(
  diversitiesInterventionStrings,
  by = c("Affinity", "PoolPatch", "InterventionPatchType")
)

diversitiesFlattened <- diversitiesFlattened %>% dplyr::mutate(
  SpeciesAffinity =
    affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)]
)

diversitiesFlattened <- diversitiesFlattened %>% dplyr::mutate(
  Value = dplyr::case_when(
    Metric == "Alpha Hill:0" & is.na(Value) ~ 0,
    TRUE ~ Value
  )
)

# diversitiesFlattened <- diversitiesFlattened %>% dplyr::filter(
#   Metric %in% c(
#     "Alpha Hill:0", "Alpha Hill:1", "Alpha Hill:Inf",
#     "TimeBrayCurtis", "TimeBrayCurtisBalance", "TimeBrayCurtisGradient",
#     "TimeJaccard", "TimeJaccardNestedness", "TimeJaccardTurnover"
#   )
# )

diversitiesFlattenedAveragedBySeed <- diversitiesFlattened %>% dplyr::group_by(
  Environment1, Environment2, Metric, Subset, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  Window = round(Time/1000)*1000
) %>% dplyr::arrange(Time) %>% dplyr::summarise(
  Mean = mean(Value),
  StDev = sqrt(var(Value)),
  Slope = if (sum(!is.na(Value)) > 1) coef(lm(Value ~ Time))[2] else {NA},
  Difference = dplyr::last(Value) - dplyr::first(Value), # Number gained.
  .groups = "drop"
) %>% dplyr::mutate(
  # Many extremely small values.
  Slope = ifelse(abs(Slope) < sqrt(.Machine$double.eps), 0, Slope)
)

save(diversitiesFlattened,
     diversitiesFlattenedAveragedBySeed,
     file = "diversitiesFlattened9_plottable.RData")

basecase <- c(
  "BasalConsumerRatio" = 0.5,
  "NSpecies" = 200,
  "PoolK1InteractionEffectiveness" = 0.01,
  "PoolK2ConsumerSizeAdvantage" = 10,
  "PoolK3ConsumerPredationRange" = 0.1,
  "PoolK4ConsumerEfficiency" = 0.2,
  # "PoolK5BasalBiomass" =  100,
  "PoolK5BasalBiomass" =  30,
  "PoolK6CoefOfVariation" = 0.1,
  "PoolBasalLogBodySize" = "c(-2, -1)",
  "PoolConsumerLogBodySize" = "c(-1, 1)",
  "InteractionK1InteractionEffectiveness" = 0.01,
  "InteractionK2ConsumerSizeAdvantage" = 10,
  "InteractionK3ConsumerPredationRange" = 0.1,
  "InteractionK4ConsumerEfficiency" = 0.2,
  # "InteractionK5BasalBiomass" =  100,
  "InteractionK5BasalBiomass" =  30,
  "InteractionK6CoefOfVariation" = 0.1,
  # "InteractionEliminationThreshold" = 1e-04,
  # "ColonizationPropaguleSize" = 0.4
  "InteractionEliminationThreshold" = 1e-05,
  "ColonizationPropaguleSize" = 0.04
)

# Add columns for whatever regression we perform,
# but remove those that aren't varying or whose variations aren't useful.
diversitiesFlattenedAveragedBySeed <-
  diversitiesFlattenedAveragedBySeed %>% dplyr::left_join(
    poolpatchDictionaryOrigin %>% dplyr::mutate(ID = as.character(ID)),
    by = c("PoolPatch" = "ID")
  ) %>% dplyr::left_join(
    dynamicsDictionaryOrigin %>% dplyr::mutate(ID = as.character(ID)),
    by = c("Interactions" = "ID")
  ) %>% dplyr::left_join(
    eventsDictionaryOrigin %>% dplyr::mutate(ID = as.character(ID)),
    by = c("Events" = "ID")
  ) %>% dplyr::select_if( # Not Varying
    .predicate = function(x) length(unique(x)) > 1
  ) %>% dplyr::select(
    -dplyr::any_of(c("PoolPatch", "Interactions", "Events")), # Not Useful After Join.
    -dplyr::ends_with("Seed")
  ) %>% dplyr::mutate(
    dplyr::across(dplyr::any_of(names(basecase)),
                  .fns = function(column) {
                    base <- basecase[[dplyr::cur_column()]]
                    base <- as(base, typeof(column[1]))
                    column <- format(c(column, base), scientific = FALSE)
                    base <- column[length(column)]
                    column <- column[-length(column)]
                    levs <- sort(unique(column))
                    relevel(
                      factor(column, levels = levs),
                      ref = base
                    )
                  })
  )

# Inside my code, it looks like PoolK1-K5 aren't called. What happens when
# we eliminate those columns? (They show up significant and large mind you,
# but they might be driven by the changes in pools and interaction matrices.)

dfabs2 <- diversitiesFlattenedAveragedBySeed %>% dplyr::select(
  -dplyr::any_of(c("PoolK1InteractionEffectiveness",
                   "PoolK2ConsumerSizeAdvantage",
                   "PoolK3ConsumerPredationRange",
                   "PoolK4ConsumerEfficiency",
                   "PoolK5BasalBiomass"))
)

metrics <- c("Mean", "StDev", "Slope", "Difference")
checks <- lapply(
  paste0(metrics, "~",
         paste0(names(diversitiesFlattenedAveragedBySeed)[
           !names(diversitiesFlattenedAveragedBySeed) %in%
             c(metrics, "Metric", "Subset")
           ],
           collapse = "+")),
  as.formula)
checks2 <- lapply(
  paste0(metrics, "~",
         paste0(names(dfabs2)[
           !names(dfabs2) %in%
             c(metrics, "Metric", "Subset")
           ],
           collapse = "+")),
  as.formula)

straightlines <- diversitiesFlattenedAveragedBySeed %>% dplyr::filter(
  Metric == "Alpha Hill:0"
)  %>% dplyr::group_by(
  Metric, Subset
) %>% dplyr::group_map(.f = function(rows, key) {
  # Regression, controlling for
  retval <- lm(checks[[1]], data = rows)
  list(key, retval)
})
straightlines2 <- dfabs2 %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::group_by(
  Metric, Subset
) %>% dplyr::group_map(.f = function(rows, key) {
  # Regression, controlling for
  retval <- lm(checks2[[1]], data = rows)
  list(key, retval)
})

stopifnot(unlist(lapply(straightlines, function(sl) {
  all(names(straightlines[[1]][[2]]$coefficients) ==
        names(sl[[2]]$coefficients))
})))

coefficientmatrix <- cbind(
  do.call(
    rbind,
    lapply(straightlines, function(x) x[[1]])
  ),
  do.call(
    rbind, lapply(straightlines, function(x)
      data.frame(t(data.frame(x[[2]]$coefficients))))
  )
)
coefficientmatrix2 <- cbind(
  do.call(
    rbind,
    lapply(straightlines2, function(x) x[[1]])
  ),
  do.call(
    rbind, lapply(straightlines2, function(x)
      data.frame(t(data.frame(x[[2]]$coefficients))))
  )
)

coefficientsdf <- tidyr::pivot_longer(
  coefficientmatrix, cols = -c(Metric, Subset),
  names_to = "CoefficientLevel",
  values_to = "Change"
) %>% dplyr::mutate(
  CoefficientLevel = ifelse(grepl("Intercept", CoefficientLevel),
                            "Intercept", CoefficientLevel)
) %>% tidyr::separate(
  col = CoefficientLevel, into = c("Coefficient", "Level"),
  fill = "right",
  sep = "(?<=[a-z])(?=[.0-9]+$)"
  # Split between letters left and ending number right
) %>% dplyr::group_by(
  Coefficient
) %>% dplyr::group_modify(
  # Use the basecase as a 0 color level. Then +1 is a single step increase.
  .f = function(value, key) {
    if (! key$Coefficient[1] %in% names(basecase)) {
      return(value %>% dplyr::mutate(StepsFromBase = 0))
    }
    basevalue <- basecase[names(basecase) == key$Coefficient[1]]
    value$Level <- gsub(pattern = "^[.]", replacement = "",
                        x = value$Level)
    uniquevalues <- unique(value$Level)
    allvalues <- c(gsub(pattern = "^[.]", replacement = "",
                        x = basevalue), uniquevalues)
    numvalues <- tryCatch({as.numeric(allvalues)},
                          warning = function(w) {
                            return(NA)
                          })
    if (!any(is.na(numvalues))) {
      allvalues <- numvalues
      value$Level <- as.numeric(value$Level)
    }
    allvalues <- sort(allvalues)
    numbers <- 1:length(allvalues) - which(basevalue == allvalues)
    value %>% dplyr::left_join(
      data.frame(Level = allvalues, StepsFromBase = numbers),
      by = c("Level")
    ) %>% dplyr::mutate(Level = as.character(Level))
  }
)
coefficientsdf2 <- tidyr::pivot_longer(
  coefficientmatrix2, cols = -c(Metric, Subset),
  names_to = "CoefficientLevel",
  values_to = "Change"
) %>% dplyr::mutate(
  CoefficientLevel = ifelse(grepl("Intercept", CoefficientLevel),
                            "Intercept", CoefficientLevel)
) %>% tidyr::separate(
  col = CoefficientLevel, into = c("Coefficient", "Level"),
  fill = "right",
  sep = "(?<=[a-z])(?=[.0-9]+$)"
  # Split between letters left and ending number right
) %>% dplyr::group_by(
  Coefficient
) %>% dplyr::group_modify(
  # Use the basecase as a 0 color level. Then +1 is a single step increase.
  .f = function(value, key) {
    if (! key$Coefficient[1] %in% names(basecase)) {
      return(value %>% dplyr::mutate(StepsFromBase = 0))
    }
    basevalue <- basecase[names(basecase) == key$Coefficient[1]]
    value$Level <- gsub(pattern = "^[.]", replacement = "",
                        x = value$Level)
    uniquevalues <- unique(value$Level)
    allvalues <- c(gsub(pattern = "^[.]", replacement = "",
                        x = basevalue), uniquevalues)
    numvalues <- tryCatch({as.numeric(allvalues)},
                          warning = function(w) {
                            return(NA)
                          })
    if (!any(is.na(numvalues))) {
      allvalues <- numvalues
      value$Level <- as.numeric(value$Level)
    }
    allvalues <- sort(allvalues)
    numbers <- 1:length(allvalues) - which(basevalue == allvalues)
    value %>% dplyr::left_join(
      data.frame(Level = allvalues, StepsFromBase = numbers),
      by = c("Level")
    ) %>% dplyr::mutate(Level = as.character(Level))
  }
)

# Overall Plot:
ggplot2::ggplot(
  coefficientsdf %>% dplyr::filter(
    # ! Metric %in% c("Alpha Hill:1", "Alpha Hill:Inf")
    Coefficient != "Intercept"
  ),
  ggplot2::aes(fill = StepsFromBase,
               x = Coefficient, y = Change)
) + ggplot2::geom_hline(
  yintercept = 0
) + ggplot2::geom_point(
  color = "black", shape = 21, size = 4
) + ggplot2::facet_grid(
  Metric ~  Subset, scales = "free_y"
) + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
) + ggplot2::scale_fill_gradient2(
  low = "orange", high = "blue", mid = "grey"
) + ggplot2::scale_y_continuous(
  minor_breaks = function(y) {seq(floor(y[1]), ceiling(y[2]))}
)
# Good news:
# Since we weren't worried about the effects of the affinity, different
# affinities should and do appear the same. Easy double check with 9a distDO.
# While we are averaging over large time periods, it appears that the window
# does not have effect, perhaps because of the size of the averaging.
# It also looks like the effects are largely the same between different metrics
# for the same alpha diversity.
ggplot2::ggplot(
  coefficientsdf %>% dplyr::filter(
    # ! Metric %in% c("Alpha Hill:1", "Alpha Hill:Inf"),
    grepl(pattern = "BodySize", Coefficient)
  ),
  ggplot2::aes(fill = Level,
               x = Coefficient, y = Change)
) + ggplot2::geom_hline(
  yintercept = 0
) + ggplot2::geom_point(
  color = "black", shape = 21, size = 4
) + ggplot2::facet_grid(
  Metric ~  Subset, scales = "free_y"
) + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
  # ) + ggplot2::scale_fill_gradient2(low = "orange", high = "blue", mid = "grey"
) + ggplot2::scale_y_continuous(
  minor_breaks = function(y) {seq(floor(y[1]), ceiling(y[2]))}
)

ggplot2::ggplot(
  diversitiesFlattenedAveragedBySeed %>% tidyr::separate(
    Subset, into = c("Type", "Affinity"), sep = "[_]"
  ) %>% dplyr::mutate(
    Type = ifelse(is.na(Type), "All", Type)
  ) %>% dplyr::filter(
    Metric == "Alpha Hill:0"
  ),
  ggplot2::aes(
    y = Mean,
    fill = InteractionEliminationThreshold,
    x = ColonizationPropaguleSize
  )
) + ggplot2::geom_boxplot(
  notch = TRUE
) + ggplot2::facet_grid(
  Metric ~  Type, scales = "free_y"
) + ggplot2::scale_y_continuous(
  minor_breaks = function(y) {seq(floor(y[1]), ceiling(y[2]))}
)

# Proper Plotting: ############################################################
plotDiversityOverview <- function(d, measures) {
  d <- d %>% dplyr::filter(
    Metric %in% measures
  )
  hits0 <- d %>% dplyr::filter(
    Time != 0
  ) %>% dplyr::group_by(
    Environment1, Environment2, Metric, Subset,
    PoolPatch, PoolPatchSeed,
    Interactions, InteractionsSeed,
    Events, EventsSeed, Dispersal, NicheDistance, Affinity,
    InterventionPatchType, InterventionPatchSeed,
    InterventionTimeType, InterventionTimeSeed,
    InterventionDispersal, InterventionNicheDistance,
    Intervention
  ) %>% dplyr::summarise(
    Any0 = any(Value == 0),
    position = max(Value),
    .groups = "drop"
  ) %>% dplyr::filter(
    Any0
  )
  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = Time,
      y = Value,
      color = interaction(Intervention)#,
      # alpha = Alpha
    )
  ) + ggplot2::geom_line(
    data = d %>% dplyr::filter(is.na(InterventionPatchType)),
    mapping = ggplot2::aes(
      group = paste(Environment1, Environment2, Metric, Subset,
                    PoolPatch, PoolPatchSeed,
                    Interactions, InteractionsSeed,
                    Events, EventsSeed, Dispersal, NicheDistance, Affinity,
                    InterventionPatchType, InterventionPatchSeed,
                    InterventionTimeType, InterventionTimeSeed,
                    InterventionDispersal, InterventionNicheDistance)
    ),
    color = "black"
  ) + ggplot2::geom_line(
    data = d %>% dplyr::filter(!is.na(InterventionPatchType)),
    mapping = ggplot2::aes(
      group = paste(Environment1, Environment2, Metric, Subset,
                    PoolPatch, PoolPatchSeed,
                    Interactions, InteractionsSeed,
                    Events, EventsSeed, Dispersal, NicheDistance, Affinity,
                    InterventionPatchType, InterventionPatchSeed,
                    InterventionTimeType, InterventionTimeSeed,
                    InterventionDispersal, InterventionNicheDistance)
    )#,
    # alpha = 0.6
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
    y = "Value", # Number of Species",
    x = paste0("Time (Characteristic Scale)"),
    color = "Intervention",
    fill = "Intervention"
  ) + ggplot2::scale_size(
    guide = "none"
  ) + ggplot2::coord_cartesian(
    ylim = c(0, NA)
  ) + ggplot2::facet_wrap(
    . ~ Metric + Subset + SpeciesAffinity
  )
  if (nrow(hits0)!=0) {
    obj + ggplot2::geom_text(
      data = hits0,
      mapping = ggplot2::aes(
        x = max(d$Time)*1.05,
        y = position,
        group = paste(Environment1, Environment2, Metric, Subset,
                      PoolPatch, PoolPatchSeed,
                      Interactions, InteractionsSeed,
                      Events, EventsSeed, Dispersal, NicheDistance, Affinity,
                      InterventionPatchType, InterventionPatchSeed,
                      InterventionTimeType, InterventionTimeSeed,
                      InterventionDispersal, InterventionNicheDistance),
        color = interaction(Intervention)
      ),
      label = "*"
    )
  } else obj
}

# e.g. plotDiversityOverview(diversitiesFlattened %>% dplyr::filter(PoolPatchSeed == 30, Intervention %in% c("(0)", "(0)->(1)"), SpeciesAffinity %in% c("evensplit_01"), NicheDistance == 6), "Alpha Hill:0")

# Or just a raw richness plot with different colors and smoothings.
ggplot(
  diversitiesFlattened %>% filter(
    Metric == "Alpha Hill:0"
  ) %>% mutate(
    Subset = factor(ifelse(is.na(Subset), "All", Subset),
                    levels = c("All",
                               "Consumer_0", "Consumer_1",
                               "Basal_0", "Basal_1"),
                    ordered = TRUE)
  ),
  aes(x = Time,
      y = Value,
      group = interaction(PoolPatchSeed, InteractionsSeed, EventsSeed,
                          InitialConditionsSeed,AffinitySeed),
      color = interaction(PoolPatchSeed, InteractionsSeed, EventsSeed,
                          InitialConditionsSeed,AffinitySeed)
  )
) + geom_line(
  show.legend = FALSE
) + geom_line(
  data = diversitiesFlattenedAveragedBySeed %>% filter(
    Metric == "Alpha Hill:0"
  ) %>% mutate(
    Subset = factor(ifelse(is.na(Subset), "All", Subset),
                    levels = c("All",
                               "Consumer_0", "Consumer_1",
                               "Basal_0", "Basal_1"),
                    ordered = TRUE
    ),
    Time = Window,
    Value = Mean),
  color = "black"
) + facet_grid(
  Subset ~ PoolPatchSeed
)

# Or the intervention overview plot:
plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("evensplit_01"),
    NicheDistance == 7,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
        )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
  )


# Hurst Exponents?

diversitiesFlattened %>% dplyr::filter(
  SpeciesAffinity %in% c("evensplit_01"),
  NicheDistance == 7,
  PoolPatchSeed == "341",
  Metric == "Alpha Hill:0",
  is.na(Subset)
) %>% dplyr::group_by(Intervention) %>% dplyr::arrange(Time) %>% dplyr::mutate(
  Hurst1 = Value - lag(Value),
  Hurst2 = Value - lag(Value,n = 2),
  Hurst4 = Value - lag(Value,n = 4),
  Hurst8 = Value - lag(Value,n = 8),
  Hurst16 = Value - lag(Value,n = 16),
  Hurst32 = Value - lag(Value,n = 32),
  Hurst64 = Value - lag(Value,n = 64),
  Hurst128 = Value - lag(Value,n = 128),
  Hurst256 = Value - lag(Value,n = 256),
  Hurst512 = Value - lag(Value,n = 512),
  Hurst1024 = Value - lag(Value,n = 1024)
) %>% dplyr::summarise(
  dplyr::across(starts_with("Hurst"), .fns = var, na.rm = TRUE)
)
