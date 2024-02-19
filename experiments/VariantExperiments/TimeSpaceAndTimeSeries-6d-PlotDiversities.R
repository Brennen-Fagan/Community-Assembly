# Introduction: ###############################################################
# Follows TimeSpaceAndTimeSeries-6c-Diversities.R
# Composed roughly of two sets of plots.
#   1) Plots corresponding to our old figure 3 of diversities.
#   2) Plots corresponding to our old figure 3 of presences, ordered by niche.

# This time around, plots will be generated within folder effectively.
datfolders <- c(
  "TSTS_Simulations_1-1_1-1_2024-02-13",
  "TSTS_Simulations_1-1_2-2_2024-02-14",
  "TSTS_Simulations_2-1_2-2_2024-02-14",
  "TSTS_Simulations_10-1_2-2_2024-02-15",
  "TSTS_Simulations_6-1_2-2_2024-02-15"
)

# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer) # Shading: stackoverflow.com/a/24436825
library(patchwork)

library(RMTRCode2)

# Load Data: ##################################################################

# results <- lapply(
#   dir(datfolders, full.names = TRUE, pattern = "Simulation"), function(x) {
#     names <- load(x)
#     stopifnot(length(names) == 1)
#     return(get(names))
#   })

# samples <- lapply(
#   dir(datfolders, full.names = TRUE, pattern = "Sampling"), function(x) {
#     names <- load(x)
#     stopifnot(length(names) == 1)
#     return(get(names))
#   })

diversities <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Diversity"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = dirname(x), "File" = basename(x)))
  })

presences <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Presence"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = dirname(x), "File" = basename(x)))
  })

stopifnot(#length(results) >= length(datfolders),
          #length(samples) == length(datfolders),
          #length(diversities) == length(results),
          length(presences) == length(diversities))

# Convert formats for plotting: ###############################################

convertThinnedDiversitiesListToDF <- function(d) {
  rbind(
    d$alpha %>% dplyr::select(
      -Species
    ) %>% tidyr::pivot_longer(
      cols = c("Richness", "Richness_Basal", "Richness_Consumer"),
      names_to = "Measurement", values_to = "Value"
    ) %>% dplyr::mutate(
      Measurement = paste("Alpha", Measurement),
      Aggregation = NA,
      Environment2 = NA,
      Affinity = if("NicheValues" %in% names(d)) lapply(d$NicheValues, as.character) else NA
    ),
    do.call(rbind, lapply(d$beta, function(b) {
      b %>% dplyr::rename(
        Environment = Env1,
        Environment2 = Env2,
        Value = Jaccard
      ) %>% dplyr::mutate(
        Measurement = "Beta Jaccard",
        Aggregation = NA,
        Affinity = if("NicheValues" %in% names(d)) lapply(d$NicheValues, as.character) else NA
      )
    })),
    d$gamma %>% dplyr::rename(
      Richness_Basal = Basals,
      Richness_Consumer = Consumers
    ) %>% tidyr::pivot_longer(
      cols =  c("Richness", "Richness_Basal", "Richness_Consumer"),
      names_to = "Measurement", values_to = "Value"
    ) %>% dplyr::mutate(
      Measurement = paste("Gamma", Measurement),
      Environment = NA,
      Environment2 = NA,
      # Maybe this needs to be a lapply if we have multiple dimensions...
      Affinity = if("NicheValues" %in% names(d)) lapply(d$NicheValues, as.character) else NA
    )
  )
}

### Diversity: ################################################################
diversities <- do.call(rbind, lapply(diversities, function(d) {
  if ("FullID" %in% names(d$Ellipsis)) {
    id <- strsplit(
      strsplit(d$Ellipsis$FullID, "_", fixed = "TRUE"), # Split seeds off.
               "-", fixed = TRUE)
  } else {
    id <- strsplit(
      strsplit(
        strsplit(d$File, ".", fixed = TRUE)[[1]][1], # Remove .RData.
        "_", fixed = TRUE)[[1]][3:4], # Remove TSTS_Type and split seeds off.
      "-", fixed = TRUE # Separate out the id values.
    )
  }

  # thinAndCalculateDiversities creates a list with alpha, beta and gamma,
  # where beta is in turn a list separated by times.
  # But, thinAndCalculateDiversities was *also* used on the different niche/
  # affinities, as well as on the whole system.

  retval <- convertThinnedDiversitiesListToDF(d$Diversities)

  retval <- rbind(
    retval,
    do.call(rbind, lapply(d$Affinity, convertThinnedDiversitiesListToDF))
  )

  retval %>% dplyr::mutate(
    PoolPatchAffinity = id[[1]][1],
    PoolPatchAffinitySeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    Dispersal = id[[1]][4],
    NicheDistance = id[[1]][5]
  )
}))

# Computationally does not seem feasible to run on the entire thing!!
diversitiesRounded <- diversities %>%  dplyr::group_by(
  round(Time), # Time
  Environment, Environment2, # Location
  Affinity, # Effectively: Species set
  Measurement, Aggregation, # Measurement

  # By Run:
  PoolPatchAffinity, PoolPatchAffinitySeed, Interactions, InteractionsSeed,
  Events, EventsSeed, Dispersal, NicheDistance
) %>% dplyr::summarise(
  Value = median(Value)
) %>% dplyr::rename(
  Time = `round(Time)`
)

diversityRibbons <- diversitiesRounded %>% dplyr::filter(
  !(is.na(Environment)), # Not Gamma
  Measurement == "Alpha Richness" |
    Measurement == "Beta Jaccard" # Not PoolType Specific
) %>%  dplyr::group_by(
  Time, # Time

  Affinity, # Effectively: Species set
  Measurement, Aggregation, # Measurement

  PoolPatchAffinity, PoolPatchAffinitySeed, Interactions, InteractionsSeed,
  Events, EventsSeed, Dispersal, NicheDistance
) %>% dplyr::summarise(
  Low = unlist(dplyr::across(dplyr::any_of("Value"),
                             .fns = ~ quantile(.x, p = 0.1, na.rm = TRUE))),
  High = unlist(dplyr::across(dplyr::any_of("Value"),
                              .fns = ~ quantile(.x, p = 0.9, na.rm = TRUE))),
  .groups = "drop"
)

# diversityRibbons_Gamma <- diversitiesRounded %>% dplyr::filter(
#   is.na(Environment),
#   Measurement == "Gamma Richness"
# ) %>%  dplyr::group_by(
#   Time, # Time
#
#   Affinity, # Effectively: Species set
#   Measurement, Aggregation, # Measurement
#
#   PoolPatchAffinity, PoolPatchAffinitySeed, Interactions, InteractionsSeed,
#   Events, EventsSeed, Dispersal, NicheDistance
# ) %>% dplyr::summarise(
#   Low = unlist(dplyr::across(dplyr::any_of("Value"),
#                              .fns = ~ quantile(.x, p = 0.1, na.rm = TRUE))),
#   High = unlist(dplyr::across(dplyr::any_of("Value"),
#                               .fns = ~ quantile(.x, p = 0.9, na.rm = TRUE))),
#   .groups = "drop"
# ) %>% dplyr::mutate(
#   Measurement = "Regional Rich."
# )

diversitiesRounded <- diversitiesRounded %>% dplyr::filter(
  Measurement %in% c(
    "Beta Jaccard", "Gamma Gamma", "Gamma Mean", "Alpha Richness"
  )
) %>% dplyr::mutate(
  Measurement2 = dplyr::case_when(
    Measurement == "Beta Jaccard" ~ "Spatial Diss.",
    Measurement == "Gamma Gamma" ~ "Regional Rich.",
    Measurement == "Gamma Mean" ~ "Local Rich.", # Panel
    Measurement == "Alpha Richness"  ~ "Local Rich.", # Otherwise
    TRUE ~ Measurement
  )
)

### Presence: #################################################################
presences <- do.call(rbind, lapply(presences, function(p) {
  if ("FullID" %in% names(p$Ellipsis)) {
    id <- strsplit(
      strsplit(p$Ellipsis$FullID, "_", fixed = "TRUE"), # Split seeds off.
      "-", fixed = TRUE)
  } else {
    id <- strsplit(
      strsplit(
        strsplit(p$File, ".", fixed = TRUE)[[1]][1], # Remove .RData.
        "_", fixed = TRUE)[[1]][3:4], # Remove TSTS_Type and split seeds off.
      "-", fixed = TRUE # Separate out the id values.
    )
  }

  retval <-

  p$SpeciesPresences %>% dplyr::mutate(
    InIntervention = Environment %in% InterventionPatches
  ) %>% dplyr::group_by(
    Species, round(Time), Environment, InIntervention
  ) %>% dplyr::summarise(
    Abundance = mean(Abundance)
  ) %>% dplyr::rename(
    Time = `round(Time)`
  ) %>% dplyr::ungroup(
  ) %>% dplyr::group_by(
    Species, Time
  ) %>% dplyr::summarise(
    Count = dplyr::n(),
    CountInControl = sum(!InIntervention),
    CountInIntervention = sum(InIntervention),
    Abundance = sum(Abundance),
    AbundanceInControl = sum(Abundance[!InIntervention]),
    AbundanceInIntervention = sum(Abundance[InIntervention])
  ) %>% dplyr::mutate(
    SimulationSet = id[1],
    InterventionType = id[2],
    InterventionParameters = id[3],
    SimulationNumber = id[4],
    ParentFileNumber = id[5]
  )



  retval %>% dplyr::mutate(
    PoolPatchAffinity = id[[1]][1],
    PoolPatchAffinitySeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    Dispersal = id[[1]][4],
    NicheDistance = id[[1]][5]
  )
}))

  retval <- convertThinnedDiversitiesListToDF(d$Diversities)

  retval <- rbind(
    retval,
    do.call(rbind, lapply(d$Affinity, convertThinnedDiversitiesListToDF))
  )

# Plotting: ###################################################################

# Inspiration: stackoverflow.com/a/24436825
# We'll take the 3rd number (which is parameters) for color,
# while we'll take the 2nd number (how we perform the transition) for shading.
categories <- aggregate(InterventionType ~ InterventionParameters,
                        diversitiesRounded, function(x) length(unique(x)))
category.palettes <- if(nrow(categories) == 3) {
  c("Oranges", "Blues", "Purples")
} else if (nrow(categories) == 2) {
  c("Oranges", "Blues")
} else if (nrow(categories) == 1) {
  c("Oranges")
} else {
  stop("More than 3 color palettes needed. We recommend you double-check.")
}
colors <- unlist(lapply(
  1:nrow(categories),
  function(i) {
    colorRampPalette(
      RColorBrewer::brewer.pal(
        9, # inverse intensity of shading
        category.palettes[i]
      )[3:7] # shades chosen
    )(categories[i, 2])
  }))
colors[1] <- "#000000" # set 1.1 to black
colorNameKeys <- sort(unique(diversitiesRounded$InterventionParameters))

### Diversity Plots: ##########################################################
# Note the ribbons are for Intervention or Intervention-Control only.
PLOT_B <- ggplot2::ggplot(
  rbind(diversitiesRounded %>% dplyr::filter(
    Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    Environment == "Mean" | Environment == "Gamma"
  ), diversitiesRounded %>% dplyr::filter(
    InIntervention == 1,
    Measurement2 %in% c("Local Rich.")
  ) %>% dplyr::group_by(
    Time,
    SimulationSet, InterventionType, InterventionParameters, # Run
    SimulationNumber, ParentFileNumber,
    Measurement2
  ) %>% dplyr::summarise(
    Value = mean(Value)
  )),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = interaction(InterventionType, InterventionParameters)
  )
) + ggplot2::geom_line(
  # alpha = 0.4,
  mapping = ggplot2::aes(
    group = interaction(InterventionType, InterventionParameters)#,
    # alpha = ifelse(Measurement2 == "Regional Rich.", 1, 0.4)
    #   )
    # ) + ggplot2::geom_line(
    #   data = diversitiesRounded %>% dplyr::filter(
    #     Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    #     Environment == "Mean"
  ),
  size = 1.5
) + ggplot2::geom_ribbon(
  data = dplyr::bind_rows(
    diversityRibbons,
    diversityRibbons_Gamma
  ) %>% dplyr::mutate(
    Measurement2 = dplyr::case_when(
      Measurement == "Jaccard" ~ "Spatial Diss.",
      Measurement == "Richness" ~ "Local Rich.",
      TRUE ~ Measurement
    )
  ),
  mapping = ggplot2::aes(
    ymin = Low,
    ymax = High,
    x = Time,
    fill = interaction(InterventionType, InterventionParameters)
  ),
  alpha = 0.1,
  inherit.aes = FALSE
) + ggplot2::theme_bw(
) + ggplot2::labs(
  y = "Value", # Number of Species",
  x = "Time (Characteristic Scale)"
  # tag = "b)"
  # x = ""
) + ggplot2::theme(
  # plot.tag.position = c(0.02, 0.98),
  strip.text.x = ggplot2::element_text(size = 8)
) + ggplot2::scale_color_manual(
  name = "Intervention",
  values = colors,
  # breaks = paste0(colorNameKeys, ".1"),
  # labels = colorNameKeys,
  aesthetics = c("colour", "fill")
) + ggplot2::facet_wrap(
  . ~ factor(
    Measurement2, ordered = T,
    levels = c("Local Rich.", "Regional Rich.", "Spatial Diss.")
  ), nrow = 1, scales = "free_y"
) + ggplot2::scale_alpha(guide = "none") + ggplot2::coord_cartesian(
  ylim = c(0, NA)
)

### Presence Plots: ###########################################################

PLOT_T <- ggplot2::ggplot(
  presences,
  ggplot2::aes(x = Time, y = Species, color = CountInIntervention)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1,
  limits = c(1, max(presences$CountInIntervention))
) + ggplot2::facet_grid(
  . ~ interaction(InterventionType, InterventionParameters)
) + ggplot2::geom_hline(
  yintercept = 1/3 * (ncol(results[[1]]$Abundance) - 1) / results[[1]]$NumEnvironments, color = "red"
) + ggplot2::labs(
  y = "Species Number",
  x = "Time (Characteristic Scale)"#,
  # tag = "a)"
) + ggplot2::theme_bw(
) + ggplot2::theme(
  # axis.text.x = ggplot2::element_blank(),
  # plot.tag.position = c(0.02, 0.98)
) + ggplot2::scale_y_continuous(
  limits = c(0, (ncol(results[[1]]$Abundance) - 1) / results[[1]]$NumEnvironments)
)

presences_1 <- presences %>% dplyr::filter(InterventionType == "1")

PLOT_T_diffs <- ggplot2::ggplot(
  presences %>% dplyr::left_join(
    presences_1 %>% dplyr::select(Time, Species, CountInIntervention),
    by = c("Time", "Species"), suffix = c("", ".1")
  ) %>% dplyr::mutate(
    InterventionDifference = CountInIntervention - ifelse(
      !is.na(CountInIntervention.1), CountInIntervention.1, 0)
  ),
  ggplot2::aes(x = Time, y = Species, color = InterventionDifference)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_distiller(
  palette = "RdYlBu",
  # direction = -1,
  # limits = c(-5, 5)
) + ggplot2::facet_grid(
  . ~ interaction(InterventionType, InterventionParameters)
) + ggplot2::geom_hline(
  yintercept = 1/3 * (ncol(results[[1]]$Abundance) - 1) / results[[1]]$NumEnvironments, color = "red"
) + ggplot2::labs(
  y = "Species Number",
  x = "Time (Characteristic Scale)"#,
  # tag = "a)"
) + ggplot2::theme_bw(
) + ggplot2::theme(
  # axis.text.x = ggplot2::element_blank(),
  # plot.tag.position = c(0.02, 0.98)
) + ggplot2::scale_y_continuous(
  limits = c(0, (ncol(results[[1]]$Abundance) - 1) / results[[1]]$NumEnvironments)
)
