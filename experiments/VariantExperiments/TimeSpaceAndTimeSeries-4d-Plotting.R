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
  controlPatches <- as.character(d$Ellipsis$Intervention$Patches)
  if (identical(controlPatches, character(0))) {
    controlPatches <- as.character(d$Ellipsis$PatchesIntervention)
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
    InControl = ((Environment1 %in% controlPatches) +
                   ifelse(!is.na(Environment2),
                          (Environment2 %in% controlPatches),
                          0))
    # I.e. Single Patch: 0 if Experiment, 1 if Control
    #      Double Patch: 0 if both Experiment, 1 if Mixed, 2 if both Control
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
  SimulationNumber, ParentFileNumber, InControl,
  Measurement # Measurement
) %>% dplyr::summarise(
  Value = median(Value)
) %>% dplyr::rename(
  Time = `round(Time)`
)

diversityRibbons <- diversitiesRounded %>% dplyr::filter(
  !(Environment %in% c("Mean", "Gamma")), # Not Gamma
  Measurement == "Richness" | Measurement == "Jaccard", # Not PoolType Specific
  (InControl == 0 & is.na(Environment2)) | (InControl == 1) # Mix/Intervention
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

# Plotting: ###################################################################

### Presence Plots: ###########################################################

### Diversity Plots: ##########################################################
PLOT_B <- ggplot2::ggplot(
  diversities %>% dplyr::filter(
    Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    Environment != "Mean"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space)
  )
) + ggplot2::geom_line(
  # alpha = 0.4,
  mapping = ggplot2::aes(
    group = interaction(Dispersal, Environment),
    alpha = ifelse(Measurement2 == "Regional Rich.", 1, 0.4)
  )
) + ggplot2::geom_line(
  data = Diversity %>% dplyr::filter(
    Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    Environment == "Mean"
  ),
  size = 1.5
) + ggplot2::geom_ribbon(
  data = dplyr::bind_rows(
    DiversityRibbons,
    DiversityRibbons_Gamma
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
    fill = pasteCustom(Dispersal, Space)
  ),
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::theme_bw(
) + ggplot2::labs(
  y = "Value", # Number of Species",
  x = paste0("Time, ", divide_time_by, " units"),
  tag = "b)"
  # x = ""
) + ggplot2::theme(
  plot.tag.position = c(0.02, 0.98),
  strip.text.x = ggplot2::element_text(size = 8)
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::scale_fill_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4")
) + ggplot2::facet_wrap(
  . ~ factor(
    Measurement2, ordered = T,
    levels = c("Local Rich.", "Regional Rich.", "Spatial Diss.")
  ), nrow = 1, scales = "free_y"
) + ggplot2::scale_alpha(guide = "none") + ggplot2::coord_cartesian(
  ylim = c(0, NA)
)

### Sampling Plots: ###########################################################
