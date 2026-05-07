# Follows MNA-Image-VikingOut-Figure4.R

# Tuftian Small Multiples Plot: ################################################
concludingplotdata_wide <- concludingplotdata %>% dplyr::select(
  -FacetPanel
) %>% tidyr::pivot_wider(
  id_cols = Set:Dispersal,
  names_from = Area, values_from = Value,
  names_expand = TRUE, values_fill = NA
)

# https://stackoverflow.com/a/53086215 for the color-group error.
addNoise <- function(x, sd = 0.1){
  xna <- is.na(x)
  xx <- x[!xna]
  anyBelow0 <- any(xx < 0)
  anyAbove1 <- any(xx > 1)
  noisy <- rnorm(length(x), mean = x, sd = sd)
  if (!anyBelow0) noisy[noisy < 0] <- 0
  if (!anyAbove1) noisy[noisy > 1] <- 1
  noisy
}

SmallMultiples <- function(
  data = concludingplotdata_wide,
  columns = 9:13, columnOffset = 8, # Lazy, could be better...
  columnLimits =
    concludingplotdata_wide[, columns] %>% lapply(range, na.rm = TRUE),
  color = "Neutral"
) {
  #columns <- 9:13
  #columnOffset <- 8
  measurements <- colnames(data)[columns]
  measurements_q <- paste0("`", measurements, "`")
  objs <- rep(list(vector("list", length = length(columns))), length(columns))
  colors <- sort(unique(data[[color]]))
  if (is.numeric(colors)) {
    scale <- ggplot2::scale_color_continuous(
      breaks = colors
    )
  } else {
    scale <- ggplot2::scale_color_discrete(
      breaks = colors
    )
  }

  for (i in columns) {
    for (j in columns) {
      if (i == j) {
        allDisps <- length(unique(data$Dispersal)) == 11
        # Dispersal Plot
        objs[[i - columnOffset]][[j - columnOffset]] <- ggplot2::ggplot(
          data,
          ggplot2::aes_string(x = "Dispersal",
                              y = measurements_q[i - columnOffset])
        ) + ggplot2::geom_boxplot(
        ) + ggplot2::theme_bw(
        ) + ggplot2::scale_y_continuous(
          limits = columnLimits[[j - columnOffset]]
        )

        if (allDisps) {
          objs[[i - columnOffset]][[j - columnOffset]] <-
            objs[[i - columnOffset]][[j - columnOffset]]  +
            ggplot2::scale_x_discrete(
              breaks = labelsforplot$Dispersal[c(1,2,4,6,8,10,11)],
              labels = c(as.character(labelsforplot$Label[c(1,2,4,6,8)]),
                         "0.59   " ,"   0.9999")
            )
        }
      } else if(i < j) {
        # Scatter Plot
        objs[[i - columnOffset]][[j - columnOffset]] <- ggplot2::ggplot(
          data %>% dplyr::mutate(
            xx = addNoise(!!rlang::sym(measurements[i - columnOffset])),
            yy = addNoise(!!rlang::sym(measurements[j - columnOffset]))
          ),
          ggplot2::aes_string(x = "xx",
                              y = "yy",
                              color = color, group = color)
        ) + ggplot2::geom_point(
          size = 1,
          alpha = 0.25#, ggplot2::aes(color = Neutral)
        ) + ggplot2::theme_bw(
        ) + ggplot2::labs(
          x = measurements[i - columnOffset],
          y = measurements[j - columnOffset]
        ) + ggplot2::scale_x_continuous(
          limits = columnLimits[[i - columnOffset]]
        ) + ggplot2::scale_y_continuous(
          limits = columnLimits[[j - columnOffset]]
        ) + scale
      } else {
        # Contours
        objs[[i - columnOffset]][[j - columnOffset]] <- ggplot2::ggplot(
          data %>% dplyr::mutate(
            xx = addNoise(!!rlang::sym(measurements[i - columnOffset])),
            yy = addNoise(!!rlang::sym(measurements[j - columnOffset]))
          ),
          ggplot2::aes_string(x = "xx",
                              y = "yy",
                              color = color, group = color)
        ) + ggplot2::geom_density_2d(
          na.rm = TRUE, bins = 3
        ) + ggplot2::theme_bw(
        ) + ggplot2::labs(
          x = measurements[i - columnOffset],
          y = measurements[j - columnOffset]
        ) + ggplot2::scale_x_continuous(
          limits = columnLimits[[i - columnOffset]]
        ) + ggplot2::scale_y_continuous(
          limits = columnLimits[[j - columnOffset]]
        ) + scale
      }
    }
  }

  objs
}

SmallMultiples_All <- SmallMultiples()
SmallMultiples_Neutral <- concludingplotdata_wide %>% dplyr::group_split(
  Neutral, .keep = TRUE
) %>% lapply(SmallMultiples, color = "Dispersal")
# e.g. wrap_plots(unlist(SmallMultiples_Neutral[[2]], recursive = FALSE), guides = "collect")

neutralorder <- concludingplotdata_wide %>% dplyr::group_split(
  Neutral, .keep = TRUE
) %>% lapply(function(x) {x %>% dplyr::pull(Neutral) %>% unique })

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "MetricComparison-SmallMultiples-All.png"
  ),
  plot = patchwork::wrap_plots(unlist(SmallMultiples_All, recursive = FALSE),
                               guides = "collect"),
  width = plot_width*2, height = plot_height*2, units = plot_units,
  dpi = plot_dpi
)

lapply(
  seq_along(SmallMultiples_Neutral),
  function(i) {
    ggplot2::ggsave(
      filename = file.path(
        outputLocation,
        paste0("MetricComparison-SmallMultiples-", neutralorder[[i]], ".png")
      ),
      plot = patchwork::wrap_plots(
        unlist(SmallMultiples_Neutral[[i]], recursive = FALSE),
        guides = "collect"),
      width = plot_width*2, height = plot_height*2, units = plot_units,
      dpi = plot_dpi
    )
  }
)

## End of Simulation Richnesses: ###############################################
# We can't use concludingplotdata since that already has the medians built in.

concluding_EOS <- dplyr::bind_rows(
  DiversitiesAlphaGamma %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = "")
    # ) %>% tidyr::unite( # we do not need, since we'll group by statistic.
    #   col = "group", Neutral, Space, remove = FALSE
  ) %>% dplyr::group_by(
    Set, Number, History, Pool, Noise, Neutral, Space, Dispersal
  ) %>% dplyr::filter(
    Time < time_burnout
  ) %>% dplyr::filter(
    Time == max(Time) # "End of Simulation", last pre-burn-out time
  ) %>% dplyr::summarise(
    # Summarising over Time and Environment so 1 value per simulation.
    `Local Richness` = median(`Richness, Alpha`, na.rm = TRUE),
    `Regional Richness` = median(`Richness, Gamma`, na.rm = TRUE),
    .groups = "drop"
  ) %>% tidyr::pivot_longer(
    names_to = "Area",
    values_to = "Value",
    cols = c("Local Richness", "Regional Richness")
  ),
  Invadabilities %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = ""),
    Area = "Median Invasibility", #"Average Invadability",
    Value = LocalMdn # Choose LocalMdn for Medians, LocalAvg for Means.
  ) %>% dplyr::select(
    -Occupancy, -LocalAvg, -LocalMdn
  ),
  TimeJaccards %>% dplyr::filter(
    Measurement == "JaccardTemporal"
  ) %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = ""),
    Area = "Temporal Jaccard",
    Value = Mdn # Choose Mdn for Medians, Avg for Means.
  ) %>% dplyr::select(
    -Measurement, -Avg, -Mdn
  ),
  DiversitiesBeta %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = "")
  ) %>% dplyr::group_by(
    Set, Number, History, Pool, Noise, Neutral, Space, Dispersal
  ) %>% dplyr::filter(
    Time < time_burnout
  ) %>% dplyr::filter(
    Time == max(Time) # "End of Simulation", last pre-burn-out time
  ) %>% dplyr::summarise(
    # Summarising over Time and Environment so 1 value per simulation.
    JaccardSpatial = median(`Jaccard`, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::mutate(
    Area = "Spatial Jaccard",
    Value = JaccardSpatial
  ) %>% dplyr::select(
    -JaccardSpatial
  )
) %>% dplyr::distinct(
) %>% dplyr::mutate(
  FacetPanel = dplyr::case_when(
    Area == "Local Richness" ~ FacetPanels[1],
    Area == "Regional Richness" ~ FacetPanels[1],
    Area == "Temporal Jaccard" ~ FacetPanels[2],
    Area == "Spatial Jaccard" ~ FacetPanels[2],
    Area == "Median Invasibility" ~ FacetPanels[3],
    Area == "Average Invasibility" ~ FacetPanels[3]
  )
)

concluding_EOS_wide <- concluding_EOS %>% dplyr::select(
  -FacetPanel
) %>% tidyr::pivot_wider(
  id_cols = Set:Dispersal,
  names_from = Area, values_from = Value,
  names_expand = TRUE, values_fill = NA
)

SmallMultiples_All <- SmallMultiples(concluding_EOS_wide)
SmallMultiples_Neutral <- concluding_EOS_wide %>% dplyr::group_split(
  Neutral, .keep = TRUE
) %>% lapply(SmallMultiples, color = "Dispersal")
# e.g. wrap_plots(unlist(SmallMultiples_Neutral[[2]], recursive = FALSE), guides = "collect")

neutralorder <- concluding_EOS_wide %>% dplyr::group_split(
  Neutral, .keep = TRUE
) %>% lapply(function(x) {x %>% dplyr::pull(Neutral) %>% unique })

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "MetricComparison-SmallMultiples-All-EOS.png"
  ),
  plot = patchwork::wrap_plots(unlist(SmallMultiples_All, recursive = FALSE),
                               guides = "collect"),
  width = plot_width*2, height = plot_height*2, units = plot_units,
  dpi = plot_dpi
)

lapply(
  seq_along(SmallMultiples_Neutral),
  function(i) {
    ggplot2::ggsave(
      filename = file.path(
        outputLocation,
        paste0("MetricComparison-SmallMultiples-", neutralorder[[i]], "-EOS.png")
      ),
      plot = patchwork::wrap_plots(
        unlist(SmallMultiples_Neutral[[i]], recursive = FALSE),
        guides = "collect"),
      width = plot_width*2, height = plot_height*2, units = plot_units,
      dpi = plot_dpi
    )
  }
)
