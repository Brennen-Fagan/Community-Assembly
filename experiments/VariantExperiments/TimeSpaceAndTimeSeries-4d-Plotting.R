# Introduction: ###############################################################
# Follows TimeSpaceAndTimeSeries-4b-Sampling.R
# Composed roughly of two sets of plots.
#   1) Plots corresponding to our old figure 3 of diversities.
#   2) Plots corresponding to our samplings, e.g. esp. File 1b.

# Warning: care with combining dat folders.
#          At minimum, the first digit should agree (same simulations).
#          The second and third digits correspond to the intervention types.
#          In particular, 2nd = 1 -> no intervention, so 3rd is meaningless.
datfolders <- c(
  "TSTS_Simulations_1-1-1_2024-01-16",
  "TSTS_Simulations_1-2-1_2024-01-10",
  "TSTS_Simulations_1-2-2_2024-01-10",
  "TSTS_Simulations_1-2-3_2024-01-10",
  "TSTS_Simulations_1-3-1_2024-01-12",
  "TSTS_Simulations_1-3-2_2024-01-12",
  "TSTS_Simulations_1-3-3_2024-01-15"#,
  # "TSTS_Simulations_2-1-6_2024-01-19",
  # "TSTS_Simulations_2-2-6_2024-01-19",
  # "TSTS_Simulations_2-3-6_2024-01-19"#,
  # "TSTS_Simulations_3-1-7_2024-01-22",
  # "TSTS_Simulations_3-2-7_2024-01-22",
  # "TSTS_Simulations_3-3-7_2024-01-22"
)

# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer) # Shading: stackoverflow.com/a/24436825
library(patchwork)

library(RMTRCode2)

# Load Data: ##################################################################

results <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Simulation"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

samples <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Sampling"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

diversities <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Diversity"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

presences <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Presence"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

stopifnot(length(results) >= length(datfolders),
          length(samples) == length(datfolders),
          length(diversities) == length(results),
          length(presences) == length(results))

# Convert formats for plotting: ###############################################

### Diversity: ################################################################
diversities <- do.call(rbind, lapply(diversities, function(d) {
  id <- strsplit(d$Ellipsis$FullID, "-", fixed = TRUE)[[1]]
  InterventionPatches <- as.character(d$Ellipsis$Intervention$Patches)
  if (identical(InterventionPatches, character(0))) {
    InterventionPatches <- as.character(d$Ellipsis$PatchesIntervention)
  }
  d$Diversities %>% dplyr::mutate(
    SimulationSet = id[1],
    InterventionType = id[2],
    InterventionParameters = id[3],
    SimulationNumber = id[4],
    ParentFileNumber = id[5]
  ) %>% tidyr::separate(
    # Acknowledge that we have two patch types and that can be important.
    Environment, into = c("Environment1", "Environment2"),
    sep = " ", remove = FALSE, fill = "right"
  ) %>% dplyr::mutate(
    InIntervention = ((Environment1 %in% InterventionPatches) +
                   ifelse(!is.na(Environment2),
                          (Environment2 %in% InterventionPatches),
                          0))
    # I.e. Single Patch: 1 if Experiment, 0 if Control
    #      Double Patch: 2 if both Experiment, 1 if Mixed, 0 if both Control
    #      Region      : 0
  )
}))

# Computationally does not seem feasible to run on the entire thing!!
diversitiesRounded <- diversities %>% dplyr::filter(
  round(Time) <= 10000
) %>%  dplyr::group_by(
  round(Time), # Time
  Environment,
  SimulationSet, InterventionType, InterventionParameters, # Run
  SimulationNumber, ParentFileNumber, InIntervention, Environment1, Environment2,
  Measurement # Measurement
) %>% dplyr::summarise(
  Value = median(Value)
) %>% dplyr::rename(
  Time = `round(Time)`
)

diversityRibbons <- diversitiesRounded %>% dplyr::filter(
  !(Environment %in% c("Mean", "Gamma")), # Not Gamma
  Measurement == "Richness" | Measurement == "Jaccard", # Not PoolType Specific
  (InIntervention == 1 & is.na(Environment2)) |
    (InIntervention == 1) # Mix/Intervention
) %>%  dplyr::group_by(
  Time, # Time
  SimulationSet, InterventionType, InterventionParameters, # Run
  SimulationNumber, ParentFileNumber,
  Measurement # Measurement
) %>% dplyr::summarise(
  Low = unlist(dplyr::across(dplyr::any_of("Value"),
                             .fns = ~ quantile(.x, p = 0.1, na.rm = TRUE))),
  High = unlist(dplyr::across(dplyr::any_of("Value"),
                              .fns = ~ quantile(.x, p = 0.9, na.rm = TRUE))),
  .groups = "drop"
)

diversityRibbons_Gamma <- diversitiesRounded %>% dplyr::filter(
  (Environment %in% c("Gamma")),
  Measurement == "Richness"
) %>%  dplyr::group_by(
  Time, # Time
  SimulationSet, InterventionType, InterventionParameters, # Run
  SimulationNumber, ParentFileNumber,
  Measurement # Measurement
) %>% dplyr::summarise(
  Low = unlist(dplyr::across(dplyr::any_of("Value"),
                             .fns = ~ quantile(.x, p = 0.1, na.rm = TRUE))),
  High = unlist(dplyr::across(dplyr::any_of("Value"),
                              .fns = ~ quantile(.x, p = 0.9, na.rm = TRUE))),
  .groups = "drop"
) %>% dplyr::mutate(
  Measurement = "Regional Rich."
)

diversitiesRounded <- diversitiesRounded %>% dplyr::mutate(
  Measurement2 = dplyr::case_when(
    Measurement == "Jaccard" ~ "Spatial Diss.",
    Measurement == "Richness" & Environment == "Gamma" ~ "Regional Rich.",
    Measurement == "Richness" & Environment == "Mean" ~ "Local Rich.", # Panel
    Measurement == "Richness"  ~ "Local Rich.", # Otherwise
    TRUE ~ Measurement
  )
)

### Presence: #################################################################
presences <- do.call(rbind, lapply(presences, function(p) {
  id <- strsplit(p$Ellipsis$FullID, "-", fixed = TRUE)[[1]]
  InterventionPatches <- as.character(p$Ellipsis$Intervention$Patches)
  if (identical(InterventionPatches, character(0))) {
    InterventionPatches <- as.character(p$Ellipsis$PatchesIntervention)
  }
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
}))

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
) + ggplot2::scale_color_viridis_c(
  # palette = "RdYlBu",
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

### Sampling Plots: ###########################################################
