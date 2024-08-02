# Plots used for BES Macro Presentation.
options(bitmapType = "cairo")

# Currently only set up for one at a time...
datfolders <- c(
  "TSTS_Simulations_18-1_6-6_2024-05-23"
)

targetFiles <- c(
  "18-1-4-15-2_6-6-6.RData", "18-1-4-15-2_6-6-6_12-1-p-3_1-1.RData"
)

targettimes <- c(0,
                 200) # 0 is intervention, ... is timespan, 50% is symmetry.

Div_rounding <- .25

timeType <- " Time Series " # "Time\nFor Time"
spaceType <- "Space\nFor Time"

histFillColor <- "#70ad47" # green
lineNullColor <- "#5b9bd5" # blue
linePertColor <-  "#ed7d31" # orange

# Libraries: ##################################################################
library(dplyr)
library(tidyr)
library(slider) # rolling windows. tried runner, but seemed buggy on a test.

library(ggplot2)
library(RColorBrewer) # Shading: stackoverflow.com/a/24436825
library(ggpubr) # replace patchwork with ggarrange.

library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R")

# Load Data: ##################################################################

samples <- unlist(lapply(targetFiles, function(tf)
  dir(datfolders, full.names = TRUE,
      pattern = paste0("Sampling_", tf))
))
samples <- lapply(
  samples, function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

diversities <- unlist(lapply(targetFiles, function(tf)
  dir(datfolders, full.names = TRUE,
      pattern = paste0("Diversity_", tf))
))
diversities <- lapply(
  diversities, function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = x))
  })

# Plotting: ###################################################################
# Parameters
targetpatches <- c(2, 1) # first patch is focal, second patch is control

### Plot of maximally distant true diversities, not separated by niche. #######
# We use the precomputed diversities. We expect experimental patch above,
# control patch below. Each panel should have two lines, one for the no
# intervention case and one for the pool swap case. We make three plots in
# total, one for each distance.

diversitiesPlottable <- do.call(
  rbind,
  diversities %>% lapply(function(d) {
    fileproperties <- strsplit(basename(d$Dir), split = "_", fixed = TRUE)[[1]]
    
    filepropertiesInit <- strsplit(fileproperties[3], split = "-",
                                   fixed = TRUE)[[1]]
    
    filepropertiesExt <- filepropertiesInit[3]
    filepropertiesExt <- ifelse(
      filepropertiesExt == "2", 1,
      ifelse(filepropertiesExt == "4", 0.9,
             ifelse(filepropertiesExt == "6", 0,
                    NA))
    )
    
    filepropertiesSpace <- filepropertiesInit[4]
    filepropertiesSpace <- ifelse(
      filepropertiesSpace == "10", 0,
      ifelse(filepropertiesSpace == "15", 5,
             ifelse(filepropertiesSpace == "NA", "Inf", NA))
    )
    
    if(length(fileproperties) == 6) {
      filepropertiesIntervention <- strsplit(fileproperties[5], split = "-",
                                             fixed = TRUE)[[1]]
      filepropertiesInterventionType <- filepropertiesIntervention[1]
      filepropertiesIntervention <- filepropertiesIntervention[4]
      filepropertiesIntervention <- ifelse(
        filepropertiesIntervention == "2", "2 ^ (1 - 2 euclid(m, n))",
        ifelse(filepropertiesIntervention == "3", "10 ^ (1 - 2 euclid(m, n))",
               "NA")
      )
    } else {
      filepropertiesIntervention <- "NA"
      filepropertiesInterventionType <- "NA"
    }
    
    # Need to unify Diversities Format
    d$Diversities$alpha <- d$Diversities$alpha %>% dplyr::select(
      -Species
    ) %>% tidyr::pivot_longer(
      cols = Richness:Richness_Consumer,
      names_to = "Aggregation", values_to = "Measurement"
    ) %>% dplyr::mutate(
      Environment = as.character(Environment)
    )
    
    d$Diversities$beta <- dplyr::bind_rows(
      d$Diversities$beta
    ) %>% dplyr::mutate(
      Environment = paste(Env1, Env2),
      Aggregation = "Jaccard",
      Jaccard = as.double(Jaccard)
    ) %>% dplyr::rename(
      Measurement = Jaccard
    ) %>% dplyr::select(
      -Env1, -Env2
    )
    
    d$Diversities$gamma <- d$Diversities$gamma %>% dplyr::mutate(
      Aggregation = ifelse(Aggregation == "Gamma", "Richness", Aggregation)
    ) %>% dplyr::rename(
      Agg = Aggregation,
      `_Gamma` = Richness,
      `_Gamma_Basal` = Basals,
      `_Gamma_Consumer` = Consumers
    ) %>% tidyr::pivot_longer(
      cols = `_Gamma`:`_Gamma_Consumer`,
      names_to = "Agg2", values_to = "Measurement"
    ) %>% dplyr::mutate(
      Environment = NA,
      Aggregation = paste0(Agg, Agg2)
    ) %>% dplyr::select(
      -Agg, -Agg2
    )
    
    d$Diversities <- dplyr::bind_rows(
      d$Diversities
    ) %>% dplyr::rename(
      Value = Measurement,
      Measurement = Aggregation
    )
    
    if (d$Ellipsis$Timescale == "Simulation") {
      d$Diversities$Time <- d$Diversities$Time / d$Ellipsis$ReactionTime
    }
    
    d$Diversities %>% dplyr::mutate(
      Time = round(Time / Div_rounding) * Div_rounding
    ) %>% dplyr::group_by(
      Time, Environment, Measurement
    ) %>% dplyr::summarise(
      Value = median(Value),
      Frequency = dplyr::n(),
      .groups = "drop"
    ) %>% dplyr::mutate(
      Space = filepropertiesSpace,
      ExtirpationProportion = filepropertiesExt,
      Intervention = ifelse(
        length(fileproperties) == 6,
        ifelse(filepropertiesInterventionType=="40",
               "(0.5, 0.5) -> (0, 1)",
               ifelse(filepropertiesInterventionType=="12",
                      "(0.5, 0.5) -> (0.5, 1)",
                      "UNLISTED")),
        "No Intervention"),
      InterventionIntensity = filepropertiesIntervention,
      ID = d$Dir
    )
  })
)

# checkGelmanRubin <- function(x) {
#   x %>% tidyr::pivot_wider(
#     id_cols = c(Time, Space, ExtirpationProportion, 
#                 Intervention, InterventionIntensity, ID),
#     names_from = c(Environment, Measurement),
#     values_from = Value
#   ) %>% dplyr::group_by(
#     Space, ExtirpationProportion, 
#     Intervention, InterventionIntensity, ID
#   ) %>% dplyr::group_map(
#     .f = function(x, y) {
#       coda::mcmc(data = x %>% dplyr::select(-Time))
#       # Thin doesn't show up in gelman.diag documentation
#     }
#     # ) %>% coda::mcmc.list(
#     # ) %>% coda::gelman.diag(
#   )
# }

# Need to duplicate data from the non-intervention case to the intervention.
postInterventionStart <- diversitiesPlottable %>% dplyr::group_by(
  Space, ExtirpationProportion, 
  Intervention, InterventionIntensity, ID
) %>% dplyr::filter(
  Intervention != "No Intervention"
) %>% dplyr::summarise(
  Time = min(Time), .groups = "drop"
)

preInterventionData <- postInterventionStart %>% dplyr::rowwise(
) %>% dplyr::group_modify(
  .f = function(x, y) {
    diversitiesPlottable %>% dplyr::filter(
      Time < x$Time,
      Space == x$Space,
      ExtirpationProportion == x$ExtirpationProportion
    ) %>% dplyr::mutate(
      Intervention = x$Intervention,
      InterventionIntensity = x$InterventionIntensity,
      ID = x$ID
    )
  }
)

windowSizes <- 2^c(1,3,5,7,9)
timeSeriesWindows <- dplyr::bind_rows(lapply(
  windowSizes, function(windowSize, dat) {
    # dat %>% runner::run_by(
    #   k = windowSize, na_pad = TRUE
    # )
    dat %>% dplyr::mutate(
      WindowSize = windowSize * Div_rounding,
      WindowAverage = slider::slide_index2(
        .x = Value, .y = Time, .f = function(value, time) {
          widths <- diff(c(time, time[1] + windowSize * Div_rounding))
          heights <- value
          return(sum(widths * heights) / (windowSize * Div_rounding))
        }, 
        .i = Time, .after = windowSize * Div_rounding, .complete = TRUE
      ),
      WindowChange = slider::slide_index(
        .x = Value, .f = function(y) {y[length(y)] - y[1]},
        .i = Time, .after = windowSize * Div_rounding, .complete = TRUE
      ),
      WindowSlope = slider::slide_index2(
        .x = Value, .y = Time, .f = 
          function(value, time) {coefficients(lm(value ~ time))[2]},
        .i = Time, .after = windowSize * Div_rounding, .complete = TRUE
      ),
      WindowTime = Time + WindowSize # Location to plot as at.
    ) %>% dplyr::mutate(
      dplyr::across(
        WindowAverage:WindowSlope, 
        .fns = function(x) {
          unlist(lapply(x, function(y) ifelse(is.null(y), NA, y)))
        }
      )
    ) %>% dplyr::filter(
      !is.na(WindowAverage)
    ) %>% dplyr::mutate(
      dplyr::across(WindowAverage:WindowSlope, 
                    .fns = function(x) ifelse(abs(x) < 1e-10, 0, x))
    )
  }, dat = dplyr::bind_rows(
    preInterventionData,
    diversitiesPlottable
  ) %>% dplyr::arrange(
    Time
  ) %>% dplyr::group_by(
    Environment, Space, ExtirpationProportion, 
    Intervention, InterventionIntensity, ID
  ) %>% dplyr::filter(
    Measurement == "Richness",
    Environment %in% targetpatches
  )
))

spaceForTimeWindows <- dplyr::bind_rows(lapply(
  windowSizes, function(windowSize, dat) {
    # In reality, we are trying to investigate the differences, 
    # so windowed average of differences makes the most sense.
    # But in practice, we'd receive a window of data which could then be either
    #   paired, then differenced, then averaged/sloped or
    #   averaged/sloped, then paired, then differenced.
    
    # Pair -> Difference -> Window -> Average/Slope
    v1 <- dat %>% dplyr::group_by(
      .add = TRUE, Time
    ) %>% dplyr::summarise(
      Value = sum(ifelse(Environment == targetpatches[2], -Value, Value)),
      .groups = "drop_last"
    ) %>% dplyr::mutate(
      WindowSize = windowSize * Div_rounding,
      DifferencesWindowAverages = slider::slide_index2(
        .x = Value, .y = Time, .f = function(value, time) {
          widths <- diff(c(time, time + WindowSize))
          heights <- value
          return(sum(widths * heights) / WindowSize)
        }, 
        .i = Time, .after = windowSize * Div_rounding, .complete = TRUE
      ),
      DifferencesWindowSlopes = slider::slide_index2(
        .x = Value, .y = Time, .f = 
          function(value, time) {coefficients(lm(value ~ time))[2]},
        .i = Time, .after = windowSize * Div_rounding, .complete = TRUE
      ),
      WindowTime = Time + WindowSize # Location to plot as at.
    )
    
    # Window -> Pair -> Difference -> Average/Slope
    # More "realistic", should match v1
    v2 <- dat %>% tidyr::pivot_wider(
      names_from = Environment,
      names_prefix = "Env",
      names_sep = "_",
      values_from = Value
    )
    colnames(v2)[colnames(v2) == paste0("Env_", targetpatches[1])] <- "Env_1"
    colnames(v2)[colnames(v2) == paste0("Env_", targetpatches[2])] <- "Env_2"
    
    v2 <- v2 %>% dplyr::mutate(
      WindowSize = windowSize * Div_rounding,
      WindowDifferencesAverages = slider::slide_index2(
        .x = Env_1, .y = Env_2, .f = function(x, y) {
          # Focal - Control
          heights <- Env_1 - Env_2
          widths <- diff(c(time, time + WindowSize))
          return(sum(widths * heights) / WindowSize)
        }
        .i = Time, .after = windowSize * Div_rounding, .complete = TRUE
      ),
      WindowDifferencesSlopes = slider::pslide_index(
        .l = list(Env_1, Env_2, Time), .f = function(x, y, Time) {
          # Focal - Control
          Value <- Env_1 - Env_2
          return(coefficients(lm(Value ~ Time))[2])
        }
        .i = Time, .after = windowSize * Div_rounding, .complete = TRUE
      ),
      WindowTime = Time + WindowSize # Location to plot as at.
    )
    
    # Window -> Average/Slope -> Pair -> Difference 
    v3 <- dat %>% dplyr::mutate(
      Value = ifelse(Environment == targetpatches[2], -Value, Value),
      WindowSize = windowSize * Div_rounding,
      WindowAveragesDifferences = runner::runner(
        x = ., f = function(y) {
          averages <- with(y, by(Value, Environment, mean))
          sum(averages)
        }, idx = Time, k = windowSize, na_pad = TRUE
      ),
      WindowSlopesDifferences = runner::runner(
        x = ., f = function(y) {
          slopes <- with(y, by(y, Environment, 
                               function(z) {
                                 coefficients(lm(Value ~ Time, data = z))[2]
                               } ))
          sum(slopes)
        },
        idx = Time, k = windowSize, na_pad = TRUE
      ),
      WindowTime = Time + WindowSize # Location to plot as at.
    )
    
    dat %>% dplyr::mutate(
      WindowAverage = runner::runner(
        x = Value, f = mean, idx = Time, k = windowSize, na_pad = TRUE
      ),
      WindowChange = runner::runner(
        x = Value, f = function(y) {y[length(y)] - y[1]},
        idx = Time, k = windowSize, na_pad = TRUE
      ),
      WindowSlope = runner::runner(
        x = ., f = function(y) coefficients(lm(Value ~ Time, data = y))[2],
        idx = Time, k = windowSize, na_pad = TRUE
      ),
      WindowSize = windowSize * Div_rounding
    ) %>% dplyr::mutate(
      dplyr::across(WindowAverage:WindowSlope, 
                    .fns = function(x) ifelse(abs(x) < 0.01, 0, x))
    ) %>% dplyr::filter(
      !is.na(WindowSlope)
    )
  }, dat = dplyr::bind_rows(
    preInterventionData,
    diversitiesPlottable
  ) %>% dplyr::arrange(
    Time
  ) %>% dplyr::group_by(
    Space, ExtirpationProportion, 
    Intervention, InterventionIntensity, ID
  ) %>% dplyr::filter(
    Measurement == "Richness",
    Environment %in% targetpatches
  )
))


ggplot2::ggplot(
  timeSeriesWindows %>% tidyr::pivot_longer(
    WindowAverage:WindowSlope, 
    names_to = "Window Function", 
    values_to = "Window Value" 
  ) %>% dplyr::group_by(
    WindowSize, `Window Function`, Intervention, Environment
  ) %>% dplyr::mutate(
    Alpha = 
      (abs(`Window Value`) > quantile(
        p = 0.99, x = abs(`Window Value`))) * 0.2 + 0.01
  ),
  ggplot2::aes(x = WindowTime,
               y = `Window Value`,
               color = factor(WindowSize, 
                              levels = windowSizes* Div_rounding, 
                              ordered = TRUE),
               group = WindowSize)
) + ggplot2::geom_vline(
  xintercept = postInterventionStart$Time,
  linetype = "dotted"
  # ) + ggplot2::geom_point(
  # ggplot2::aes(alpha = Alpha)
  #   shape = '.'
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  `Window Function` ~ Intervention + Environment  ,
  scales = "free_y"
) + ggplot2::labs(
  color = "Window Size"
  # ) + ggplot2::scale_y_continuous(
  #   transform = scales::pseudo_log_trans(sigma = 0.1)
) + ggplot2::guides(
  color = ggplot2::guide_legend(override.aes = list(alpha = 1))
)

ggplot2::ggplot(
  timeSeriesWindows %>% tidyr::pivot_longer(
    WindowAverage:WindowSlope, 
    names_to = "Window Function", 
    values_to = "Window Value" 
  ) %>% dplyr::group_by(
    WindowSize, `Window Function`, Intervention, Environment
  ) %>% dplyr::mutate(
    Alpha = 
      (abs(`Window Value`) > quantile(
        p = 0.99, x = abs(`Window Value`))) * 0.2 + 0.01,
    `Window Value` = `Window Value` / max(abs(`Window Value`))
  ),
  ggplot2::aes(x = WindowTime,
               y = `Window Value`,
               color = factor(WindowSize, 
                              levels = windowSizes* Div_rounding, 
                              ordered = TRUE),
               group = WindowSize)
) + ggplot2::geom_vline(
  xintercept = postInterventionStart$Time,
  linetype = "dotted"
) + ggplot2::geom_point(
  ggplot2::aes(alpha = Alpha)
  #   shape = '.'
  # ) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  `Window Function` ~ Intervention + Environment  ,
  scales = "free_y"
) + ggplot2::labs(
  color = "Window Size"
  # ) + ggplot2::scale_y_continuous(
  #   transform = scales::pseudo_log_trans(sigma = 0.1)
) + ggplot2::guides(
  color = ggplot2::guide_legend(override.aes = list(alpha = 1))
)

# Examine more closely the Change and Slope for sensitivity ~ WindowSize.
ggplot2::ggplot(
  timeSeriesWindows %>% tidyr::pivot_longer(
    WindowChange:WindowSlope, 
    names_to = "Window Function", 
    values_to = "Window Value" 
  ),
  ggplot2::aes(x = WindowTime,
               y = `Window Value` * WindowSize,
               color = factor(WindowSize, 
                              levels = windowSizes* Div_rounding, 
                              ordered = TRUE),
               group = WindowSize)
) + ggplot2::geom_vline(
  xintercept = postInterventionStart$Time,
  linetype = "dotted"
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  `Window Function` ~ Intervention + Environment  ,
  scales = "free_y"
) + ggplot2::labs(
  color = "Window Size"
)
# Slope magnitude maximised when even on either side it looks like.
# This indicates it's more of a late warning signal.
