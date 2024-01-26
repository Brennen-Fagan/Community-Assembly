files <- c(
file.path("Data_2023-07-06",
          "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData"),
  file.path("Data_2024-01-17",
            "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData"),
  file.path("Data_2024-01-18",
            "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData")
)




cores <- 3


# Libraries: ##################################################################
library(dplyr)
library(RMTRCode2)
library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

# Load Data: ##################################################################

Diversity <- foreach::foreach(
  x = iterators::iter(
    files
  ), .packages = c("dplyr", "RMTRCode2")
) %op% {
  # x_properties <- strsplit(basename(x), split = "_")
  # stopifnot(length(x_properties) == 1,
  #           x_properties[[1]][1] == "TSTS",
  #           x_properties[[1]][2] == "Simulation")
  # filename <- file.path(
  #   dirname(x),
  #   paste0(x_properties[[1]][1],
  #          "_Diversity_",
  #          x_properties[[1]][3])
  # )
  # 
  # if(file.exists(filename)) {
  #   load(filename)
  # } else {
    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects
    
    numberOfSpecies <- (ncol(loaded$Abundance) - 1)/loaded$NumEnvironments
    Diversity <- list(
      Diversities = RMTRCode2::Calculate_Diversity(
        loaded,
        nspecies = c(Basal = 34, Consumer = 66) * numberOfSpecies / 100
        # My standard approach for nspecies.
      ),
      Ellipsis = x
    )
    
    # save(Diversity, file = filename)
  # }
  
  return(Diversity)
}


# Cleanup: ####################################################################
if (exists("clust")) {
  parallel::stopCluster(clust)
}

diversities2 <- do.call(rbind, lapply(Diversity, function(d) {
  id <- d$Ellipsis
  d$Diversities %>% dplyr::mutate(
    Directory = dirname(id),
    File = basename(id)
  ) %>% tidyr::separate(
    # Acknowledge that we have two patch types and that can be important.
    Environment, into = c("Environment1", "Environment2"),
    sep = " ", remove = FALSE, fill = "right"
  )
}))

# Computationally does not seem feasible to run on the entire thing!!
diversitiesRounded <- diversities2 %>% dplyr::filter(
  round(Time) <= 10000
) %>%  dplyr::group_by(
  round(Time), # Time
  Environment,Environment1, Environment2,
  Directory, File,
  Measurement # Measurement
) %>% dplyr::summarise(
  Value = median(Value)
) %>% dplyr::rename(
  Time = `round(Time)`
)

diversityRibbons <- diversitiesRounded %>% dplyr::filter(
  !(Environment %in% c("Mean", "Gamma")), # Not Gamma
  Measurement == "Richness" | Measurement == "Jaccard" # Not PoolType Specific
) %>%  dplyr::group_by(
  Time, # Time
  Directory, File,
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
  Directory, File,
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

# Plotting:
PLOT_B <- ggplot2::ggplot(
  rbind(diversitiesRounded %>% dplyr::filter(
    Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    Environment == "Mean" | Environment == "Gamma"
  ), diversitiesRounded %>% dplyr::filter(
    Measurement2 %in% c("Local Rich.")
  ) %>% dplyr::group_by(
    Time,
    Directory, File,
    Measurement2
  ) %>% dplyr::summarise(
    Value = mean(Value)
  )),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = interaction(Directory, File)
  )
) + ggplot2::geom_line(
  # alpha = 0.4,
  mapping = ggplot2::aes(
    group = interaction(Directory, File)#,
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
    fill = interaction(Directory, File)
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
# ) + ggplot2::scale_color_manual(
#   name = "Run",
  # values = colors,
  # breaks = paste0(colorNameKeys, ".1"),
  # labels = colorNameKeys,
  # aesthetics = c("colour", "fill")
) + ggplot2::facet_grid(
  factor(
    Measurement2, ordered = T,
    levels = c("Local Rich.", "Regional Rich.", "Spatial Diss.")
  ) ~ Directory, scales = "free_y"
) + ggplot2::scale_alpha(guide = "none") + ggplot2::coord_cartesian(
  ylim = c(0, NA)
)
