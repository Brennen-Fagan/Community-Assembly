# Introduction: ###############################################################
# Follows TimeSpaceAndTimeSeries-6c-Diversities.R applied to 7b and 7d.
#   1) Plots corresponding to single patch systems.
#   2) Plots corresponding to double patch systems.

datfolders <- dir(pattern = "TSTS_Simulations_")

savedir <- file.path(datfolders, "Images")
lapply(savedir, function(savedir) if (!dir.exists(savedir)) {
  dir.create(savedir)
})
# Problems with X11
options(bitmapType = "cairo")

# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer) # Shading: stackoverflow.com/a/24436825
library(patchwork)

library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Dictionaries.R")
source('TimeSpaceAndTimeSeries-0-Functions.R')

# Load Data: ##################################################################
diversities <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Diversity"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = dirname(x), "File" = basename(x)))
  })


# Desired density:
pointsPerTimeUnit <- 1/10

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
  ) %>% dplyr::filter(!is.na(Weights), !is.na(Value))
  
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
) %>% tidyr::separate(
  Measurement, into = c("Measurement", "Species Layer"), sep = "_",
  fill = "right"
)

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
  ),
  Intervention = unlist(mapply(FUN = function(ppa, ipt) {
    ppDO <- poolpatchDictionaryOrigin[ppa, ]
    
    if (explicit <- grepl(pattern = "rep", ppDO$PatchAffinities)) {
      initState <- 
        paste0("(", paste( # NOT PRETTY FOR 10, MAY WANT TO JUST REPORT FUNC CALL
          vals <- retrieveFunction(ppDO$PatchAffinities)(ppDO$NumberEnvironments), 
          collapse = ", "), ")")
    } else {
      initState <- 
        paste0(ppDO$PatchAffinities, "(", ppDO$NumberEnvironments, ")")
    }
    
    if(is.na(ipt)) {return(initState)}
    
    ipDO <- interventionPatchDictionaryOrigin[ipt, ]
    
    if (is.na(ipDO$InterventionLocation) || 
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
  }, PoolPatchAffinity, InterventionPatchType)),
  Alignment = dplyr::case_when(
    is.na(Affinity) ~ "All Species",
    Affinity == 0 ~ "Type 0 Species",
    Affinity == 1 ~ "Type 1 Species"
  )
)

# Plotting: ###################################################################

plotDiversityOverview <- function(dRounded) {
  ggplot2::ggplot(
    dRounded %>% dplyr::filter(
      Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
      is.na(Aggregation) | Aggregation == "Gamma"
    )  %>% dplyr::mutate(
      Affinity = unlist(Affinity)
    ),
    ggplot2::aes(
      x = Time,
      y = Value,
      color = interaction(Intervention),
      alpha = Alpha
    )
  ) + ggplot2::geom_line(
    # alpha = 0.4,
    mapping = ggplot2::aes(
      group = paste(Environment, Environment2, Affinity, Measurement, Aggregation,
                    PoolPatchAffinity, PoolPatchAffinitySeed,
                    Interactions, InteractionsSeed,
                    Events, EventsSeed, Dispersal, NicheDistance,
                    InterventionPatchType, InterventionPatchSeed,
                    InterventionTimeType, InterventionTimeSeed,
                    InterventionDispersal, InterventionNicheDistance),
      # alpha = ifelse(Measurement2 == "Regional Rich.", 1, 0.4)#,
      # size = ifelse(Measurement2 == "Regional Rich.", 1.2, 1)
    ),
    # alpha = 0.6
  ) + ggplot2::geom_line(
    data = dRounded %>% dplyr::filter(
      Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
      Aggregation %in% c("Mean", "Gamma")
    ),
    size = 1
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
    y = "Value", # Number of Species",
    x = paste0("Time (Characteristic Scale)"),
    color = "Intervention",
    fill = "Intervention"#,
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
  # ) + ggplot2::facet_grid(
  #   factor(
  #     Measurement2, ordered = T,
  #     levels = c("Local Rich.", "Regional Rich.", "Spatial Diss.")
  #   ) ~ #PoolPatchAffinity +
  #     #Affinity ,# ncol = 3,
  #     Alignment,
  #   scales = "free_y"
  # ) + ggplot2::scale_alpha(
  #   guide = "none"
  ) + ggplot2::scale_size(
    guide = "none"
  ) + ggplot2::coord_cartesian(
    ylim = c(0, NA)
  )}
