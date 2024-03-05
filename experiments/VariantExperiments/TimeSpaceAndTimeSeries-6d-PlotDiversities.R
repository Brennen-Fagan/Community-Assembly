# Introduction: ###############################################################
# Follows TimeSpaceAndTimeSeries-6c-Diversities.R
# Composed roughly of two sets of plots.
#   1) Plots corresponding to our old figure 3 of diversities.
#   2) Plots corresponding to our old figure 3 of presences, ordered by niche.

# This time around, plots will be generated within folder effectively.
datfolders <- c(
  # "TSTS_Simulations_1-1_1-1_2024-02-13",
  # "TSTS_Simulations_1-1_2-2_2024-02-14",
  # "TSTS_Simulations_2-1_2-2_2024-02-14",
  # "TSTS_Simulations_10-1_2-2_2024-02-15",
  # "TSTS_Simulations_6-1_2-2_2024-02-15"
  "TSTS_Simulations_11-1_4-4_2024-02-23"
)
# Problems with X11
options(bitmapType = "cairo")

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

# Fixed data: #################################################################
if (dplyr::is_grouped_df(presences[[1]]$SpeciesPresences)) {
  presences <- lapply(presences, function(p) {
    p$SpeciesPresences <- dplyr::ungroup(p$SpeciesPresences)
    p
  })
}

# Desired density:
pointsPerTimeUnit <- 1/100

# Convert formats for plotting: ###############################################

convertThinnedDiversitiesListToDF <- function(
  d, pPTU = pointsPerTimeUnit
) {
  # Shared Format
  retval <- rbind(
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

  # Balanced Thinning

  retval <- retval %>% dplyr::arrange(
    Time
  ) %>% dplyr::mutate( # Thin according to weighted time grouping.
    TimeGroup = floor(Time * pPTU) / pPTU
  ) %>% dplyr::group_by(
    TimeGroup, Measurement, Aggregation, Environment, Environment2, Affinity
  ) %>% dplyr::group_modify(
    .f = function(.x, .y) {
      ## Add beginning and end of time group:
      rbind(
        if(!unname(.y$TimeGroup) %in% .x$Time)
          data.frame(Time = unname(.y$TimeGroup),
                     Value = NA),
        .x,
        if(!any(.x$Time > unname(.y$TimeGroup) + 0.99/pPTU))
          data.frame(Time = unname(.y$TimeGroup) + 0.99/pPTU,
                   Value = .x[nrow(.x),]$Value)
        )
    }
  ) %>% dplyr::ungroup(
  ) %>% dplyr::group_by(
    Measurement, Aggregation, Environment, Environment2, Affinity
  ) %>% dplyr::mutate(
    Value = ifelse(is.na(Value), dplyr::lag(Value), Value), # All but first
    Weights = c(diff(Time), NA)
  ) %>% dplyr::ungroup(
  )

  ## Summarize
  retval %>% dplyr::group_by(
    TimeGroup, Measurement, Aggregation, Environment, Environment2, Affinity
  ) %>% dplyr::summarise(
    Value = Hmisc::wtd.quantile(Value, Weights, normwt = TRUE, probs = 0.5),
    #Value = sum(Weights * Value) / sum(Weights), # Mean
    Time = unique(TimeGroup)[1],
    .groups = "drop"
  ) %>% dplyr::select(-TimeGroup)
}

### Diversity: ################################################################
diversities <- do.call(rbind, lapply(diversities, function(d) {
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
    NicheDistance = id[[1]][5],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
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
  Events, EventsSeed, Dispersal, NicheDistance,
  InterventionPatchType, InterventionPatchSeed, InterventionTimeType,
  InterventionTimeSeed, InterventionDispersal, InterventionNicheDistance
) %>% dplyr::summarise(
  Value = median(Value),
  .groups = "drop"
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
  Events, EventsSeed, Dispersal, NicheDistance,
  InterventionPatchType, InterventionPatchSeed, InterventionTimeType,
  InterventionTimeSeed, InterventionDispersal, InterventionNicheDistance
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
    "Beta Jaccard", "Gamma Richness", "Alpha Richness"
  )
) %>% dplyr::mutate(
  Measurement2 = dplyr::case_when(
    Measurement == "Beta Jaccard" ~ "Spatial Diss.",
    Measurement == "Gamma Richness" & Aggregation == "Gamma" ~ "Regional Rich.",
    Measurement == "Gamma Richness" & Aggregation == "Mean" ~ "Local Rich.", # Panel
    Measurement == "Alpha Richness"  ~ "Local Rich.", # Otherwise
    TRUE ~ Measurement
  )
)

### Presence: #################################################################
presences <- do.call(rbind, lapply(
  presences,
  function(p, pPTU = pointsPerTimeUnit) {
    if ("FullID" %in% names(p$Ellipsis)) {
      id <- strsplit(
        strsplit(p$Ellipsis$FullID, "_", fixed = TRUE)[[1]], # Split seeds off.
        "-", fixed = TRUE)
    } else if ("ID" %in% names(p$Ellipsis)) {
      id <- strsplit(
        strsplit(p$Ellipsis$ID, "_", fixed = TRUE)[[1]], # Split seeds off.
        "-", fixed = TRUE)
    } else {
      id <- strsplit(
        strsplit(
          strsplit(p$File, ".", fixed = TRUE)[[1]][1], # Remove .RData.
          "_", fixed = TRUE)[[1]][3:4], # Remove TSTS_Type and split seeds off.
        "-", fixed = TRUE # Separate out the id values.
      )
    }

    if (length(id) < 3) {
      # I.e., no intervention.
      id[[3]] <- rep(NA, 4)
      id[[4]] <- rep(NA, 2)
    }


  if ("TimeFloor" %in% names(p$SpeciesPresences)) {
    p$SpeciesPresences <- p$SpeciesPresences %>% dplyr::rename(
      Time = TimeFloor
    )
  }

  # retval <- p$SpeciesPresences %>% dplyr::group_by(
  #   Species, round(Time), Environment, dplyr::contains("Affinity"), Size
  # ) %>% dplyr::summarise(
  #   Abundance = mean(Abundance),
  #   Biomass = Abundance *  Size
  # ) %>% dplyr::rename(
  #   Time = `round(Time)`
  # )

  # Balanced Thinning

    minTimeDiff <- p$SpeciesPresences %>% dplyr::arrange(
      Time
    ) %>% dplyr::mutate( # Thin according to weighted time grouping.
      TimeGroup = floor(Time * pPTU) / pPTU
    ) %>% dplyr::group_by(
      TimeGroup, Species, Environment, Size, Type, Affinity, EnvAffinity
    ) %>% dplyr::mutate(
      Weights = c(diff(Time), NA)
    ) %>% dplyr::pull(Weights) %>% min(na.rm = TRUE)

  retval <- p$SpeciesPresences %>% dplyr::arrange(
    Time
  ) %>% dplyr::mutate( # Thin according to weighted time grouping.
    TimeGroup = floor(Time * pPTU) / pPTU
  ) %>% dplyr::group_by(
    TimeGroup, Species, Environment, Size, Type, Affinity, EnvAffinity
  ) %>% dplyr::group_modify(
    .f = function(.x, .y) {
      ## Add beginning and end of time group:
      rbind(
        if(!unname(.y$TimeGroup) %in% .x$Time)
          data.frame(Time = unname(.y$TimeGroup),
                     Abundance = 0), # in presence, only present if nonzero.
        #                          but we're now feeling in for averaging.
        .x,
        if(!any(.x$Time > unname(.y$TimeGroup) + 0.99/pPTU))
          data.frame(Time = c(max(.x$Time) + minTimeDiff,
                              unname(.y$TimeGroup) + 0.99/pPTU),
                     Abundance = 0) # If Far => 0, if Near => keep prev value, but
        #                       # if near, then there should be a value nearby.
      ) # Solution isn't ideal: e.g. 1, 2, 0, ...., 0, 1 gives 2 high weight.
    }   # This does cover 1, 2, 0, ..., 0, 0 though.
  ) %>% dplyr::ungroup(
  ) %>% dplyr::group_by(
    Species, Environment, EnvAffinity
  ) %>% dplyr::mutate(
    Weights = c(diff(Time), NA)
  ) %>% dplyr::ungroup(
  ) %>% dplyr::filter(!is.na(Weights))

  retval %>% dplyr::group_by(
    TimeGroup, Species, Environment, Size, Type, Affinity, EnvAffinity
  ) %>% dplyr::filter(
    sum(!is.na(Weights)) > 0
  ) %>% dplyr::summarise(
    Abundance = Hmisc::wtd.quantile(Abundance, Weights, normwt = TRUE, probs = 0.5),
    Abundance = ifelse(Abundance < 1e-4, 0, Abundance),
    #Value = sum(Weights * Value) / sum(Weights), # Mean
    Time = unique(TimeGroup)[1],
    .groups = "drop"
  ) %>% dplyr::select(
    -TimeGroup
  ) %>% dplyr::mutate(
    PoolPatchAffinity = id[[1]][1],
    PoolPatchAffinitySeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    Dispersal = id[[1]][4],
    NicheDistance = id[[1]][5],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
  )
}))

diversity_cuts <- diversities %>% dplyr::select(
  Affinity,
  PoolPatchAffinity, PoolPatchAffinitySeed,
  Interactions, InteractionsSeed,
  Events, EventsSeed,
  Dispersal, NicheDistance,
  InterventionPatchType, InterventionPatchSeed, InterventionTimeType,
  InterventionTimeSeed, InterventionDispersal, InterventionNicheDistance
) %>% dplyr::distinct() %>% dplyr::mutate(
  Affinity = unlist(Affinity)
) %>% dplyr::filter(
  !is.na(Affinity)
  ) %>% tidyr::separate(
    sep = ",", col = Affinity, into = c("Affinity.Lower", "Affinity.Upper"),
    fill = "right"
  ) %>% dplyr::mutate(
    Affinity.Upper = ifelse(is.na(Affinity.Upper),
                            Affinity.Lower,
                            Affinity.Upper),
    Affinity.Lower = gsub(pattern = "([(]|[]])", replacement = "", Affinity.Lower),
    Affinity.Upper = gsub(pattern = "([(]|[]])", replacement = "", Affinity.Upper),
    Affinity.Lower = as.numeric(Affinity.Lower),
    Affinity.Upper = as.numeric(Affinity.Upper)
  )

presences <- presences %>% dplyr::left_join(
  diversity_cuts, by = c("PoolPatchAffinity", "PoolPatchAffinitySeed",
                         "Interactions", "InteractionsSeed",
                         "Events", "EventsSeed",
                         "Dispersal", "NicheDistance",
                         "InterventionPatchType", "InterventionPatchSeed",
                         "InterventionTimeType", "InterventionTimeSeed",
                         "InterventionDispersal", "InterventionNicheDistance")
) %>% dplyr::filter(
  Affinity.Lower <= Affinity,
  Affinity <= Affinity.Upper
)

# Plotting: ###################################################################

### Diversity Plots: ##########################################################
# This works, but just barely...
# TODO: Get separate scales (-> convert back to facet_wrap, maybe + patchwork).
# TODO: Decide how to separate the nichedistance, which is sometimes important.
# TODO: Color scales.
# TODO: Add back in the intervals and median lines.
PLOT_B <- ggplot2::ggplot(
  diversitiesRounded %>% dplyr::filter(
    Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    is.na(Aggregation) | Aggregation == "Gamma"
  )  %>% dplyr::mutate(
    Affinity = unlist(Affinity)
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = interaction(Dispersal, NicheDistance)
  )
) + ggplot2::geom_line(
  # alpha = 0.4,
  mapping = ggplot2::aes(
    group = paste(Dispersal, NicheDistance,
                  PoolPatchAffinity, PoolPatchAffinitySeed, Affinity,
                  Environment, Environment2,
                  Measurement2),
    alpha = ifelse(Measurement2 == "Regional Rich.", 1, 0.4)
  )
# ) + ggplot2::geom_line(
#   data = diversitiesRounded %>% dplyr::filter(
#     Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
#     Aggregation == "Mean"
#   ),
#   size = 1.5
# ) + ggplot2::geom_ribbon(
#   data = dplyr::bind_rows(
#     DiversityRibbons,
#     DiversityRibbons_Gamma
#   ) %>% dplyr::mutate(
#     Measurement2 = dplyr::case_when(
#       Measurement == "Jaccard" ~ "Spatial Diss.",
#       Measurement == "Richness" ~ "Local Rich.",
#       TRUE ~ Measurement
#     )
#   ),
#   mapping = ggplot2::aes(
#     ymin = Low,
#     ymax = High,
#     x = Time,
#     fill = pasteCustom(Dispersal, Space)
#   ),
#   alpha = 0.4,
#   inherit.aes = FALSE
) + ggplot2::theme_bw(
) + ggplot2::labs(
  y = "Value", # Number of Species",
  x = paste0("Time (Characteristic Scale)"),
  color = "Dispersal\nNicheDistance"#,
  # tag = "(b)"
  # x = ""
# ) + ggplot2::theme(
#   plot.tag.position = c(0.02, 0.98),
#   plot.tag = ggplot2::element_text(face = "bold"),
#   strip.text.x = ggplot2::element_text(size = 8)
# ) + ggplot2::scale_color_manual(
#   name = legend_bl_name,
#   values = c("darkorange", "plum1", "cyan")
# ) + ggplot2::scale_fill_manual(
#   name = legend_bl_name,
#   values = c("darkorange4", "plum4", "cyan4")
) + ggplot2::facet_grid(
  factor(
    Measurement2, ordered = T,
    levels = c("Local Rich.", "Regional Rich.", "Spatial Diss.")
  ) ~ PoolPatchAffinity + Affinity ,# ncol = 3,
  scales = "free_y"
) + ggplot2::scale_alpha(guide = "none") + ggplot2::coord_cartesian(
  ylim = c(0, NA)
)

### Presence Plots: ###########################################################

# As with dispersal plots, this is the bare minimum to be functional.
PLOT_T <- ggplot2::ggplot(
  presences %>% dplyr::filter(
    Abundance > 0
  ) %>% dplyr::group_by(
    Species, Size, Type, Affinity.Lower, Affinity.Upper,
    Time,
    PoolPatchAffinity, PoolPatchAffinitySeed,
    Interactions, InteractionsSeed,
    Events, EventsSeed,
    Dispersal, NicheDistance
  ) %>% dplyr::summarise(
    Count = dplyr::n(),
    .groups = "drop"
  ),
  ggplot2::aes(x = Time, y = Size, color = Count)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1,
  limits = c(1, 10)
) + ggplot2::facet_grid(
  Dispersal ~ paste(PoolPatchAffinitySeed,
              NicheDistance,#) ~ paste(
                PoolPatchAffinity, Affinity.Lower, Affinity.Upper)
# ) + ggplot2::geom_hline(
#   yintercept = 1/3 * (ncol(results[[1]]$Abundance) - 1) / results[[1]]$NumEnvironments, color = "red"
) + ggplot2::labs(
  y = "Species Size",
  x = "Time (Characteristic Scale)"#,
  # tag = "a)"
) + ggplot2::theme_bw(
) + ggplot2::theme(
  # axis.text.x = ggplot2::element_blank(),
  # plot.tag.position = c(0.02, 0.98)
) + ggplot2::scale_y_continuous(
  limits = c(0, 1)
)
