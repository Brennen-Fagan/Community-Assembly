datfolders <- dir(pattern = "TSTS_Simulations_")

# Problems with X11
options(bitmapType = "cairo")

# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)

source("TimeSpaceAndTimeSeries-0-Dictionaries.R")
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


# Fix issue with the wrong name for "Time" (time).
# with(diversities[[1]]$Diversity %>% dplyr::mutate(
#   Time2 = Time,
#   Time = ifelse(is.na(time), 0, time) +
#     ifelse(is.na(Time), 0, Time)),
#   all((is.na(Time2) & Time == time) | (is.na(time) & Time == Time2))
# )
diversities <- lapply(
  diversities,
  function(d) {
    if ("time" %in% names(d$Diversity) &&
        "Time" %in% names(d$Diversity)) {
      d$Diversity <- d$Diversity %>% dplyr::mutate(
        Time = ifelse(is.na(time), 0, time) +
          ifelse(is.na(Time), 0, Time)
      ) %>% dplyr::select(-time) %>% dplyr::select(
        Time, dplyr::everything()
      )
    }
    return(d)
  }
)

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

  d$Diversity %>% dplyr::mutate(
    PoolPatch = id[[1]][1],
    PoolPatchSeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    Dispersal = id[[1]][4],
    NicheDistance = id[[1]][5],
    Affinity = id[[1]][6],
    AffinitySeed = id[[2]][4],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
  )
}))

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

# save(diversitiesFlattened, file = "diversitiesFlattened_plottable.RData")

diversitiesFlattenedAveragedBySeed <- diversitiesFlattened %>% dplyr::group_by(
  Environment1, Environment2, Metric, Subset, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  Window = round(Time/10)*10
) %>% dplyr::arrange(Time) %>% dplyr::summarise(
  Mean = mean(Value),
  Slope = if (sum(!is.na(Value)) > 1) coef(lm(Value ~ Time))[2] else {NA},
  Difference = dplyr::first(Value) - dplyr::last(Value),
  .groups = "drop"
)

diversitiesFlattenedAveragedAcrossSeed <- diversitiesFlattenedAveragedBySeed %>%
  dplyr::group_by(
    Environment1, Environment2, Metric, Subset, PoolPatch,
    Interactions, Events, Dispersal, NicheDistance,
    Affinity, InterventionPatchType,
    InterventionTimeType, InterventionDispersal,
    InterventionNicheDistance, Intervention, SpeciesAffinity,
    Window
  ) %>% dplyr::summarise(
    Mean = mean(Mean),
    Slope = mean(Slope),
    Difference = mean(Difference)
  )

# Proper Plotting: ############################################################
plotDiversityOverview <- function(d, measures) {
  d <- d %>% dplyr::filter(
    Metric %in% measures
  )
  ggplot2::ggplot(
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
}

# e.g. plotDiversityOverview(diversitiesFlattened %>% dplyr::filter(PoolPatchSeed == 30, Intervention %in% c("(0)", "(0)->(1)"), SpeciesAffinity %in% c("evensplit_01"), NicheDistance == 6), "Alpha Hill:0")
