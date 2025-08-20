GammaOverAlpha <- function(
  i, d, levels,
  spacelabs, neutrallabs,
  burnin = time_burnin,
  burnout = time_burnout,
  # For each parameter set, there are 10 Pools, and in each...
  # Pool (for which there are 10 Histories),
  # History (for which there are 10 Environments), Environment
  NumberSelect = NULL ,#1,
  HistorySelect = NULL, #1,
  EnvirSelect = NULL #1
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]
  
  d <- d %>% dplyr::filter(
    Pool == pbs,
    Noise == inm
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  )
  
  if (!is.null(NumberSelect) &&
      !is.null(HistorySelect) &&
      !is.null(EnvirSelect)) {
    number <- d %>% dplyr::select(
      Space, Neutral, Number, History, Environment
    ) %>% dplyr::distinct() %>% dplyr::arrange(
      Number
    ) %>% dplyr::group_by(
      Space, Neutral
    ) %>% dplyr::filter(
      Number %in% (min(Number) - 1 + NumberSelect)
    ) %>% dplyr::group_by(
      Space, Neutral, Number
    ) %>% dplyr::filter(
      History %in% (min(History) - 1 + HistorySelect)
    ) %>% dplyr::group_by(
      Space, Neutral, Number, History
    ) %>% dplyr::filter(
      Environment %in% (min(Environment) - 1 + EnvirSelect)
    )
    
    dTraj <- d %>% dplyr::semi_join(number, by = c(
      "Space", "Neutral", "Number", "History", "Environment"
    ))
    
    dTrajStop <- dTraj %>% dplyr::group_by(
      Number, History, Environment
    ) %>% dplyr::slice_max(Time)
    
  } else if (!is.null(NumberSelect) ||
             !is.null(HistorySelect) ||
             !is.null(EnvirSelect)) {
    warning("Not all of Selects are null or non-null. Skipping trajectory.")
  }
  
  dga <- d %>% dplyr::filter(
    Time > burnin, Time < burnout
  )
  
  obj <- ggplot2::ggplot(
    dga,
    ggplot2::aes(x = `Richness, Alpha`,
                 y = `Richness, Gamma`)
  ) + ggplot2::geom_bin2d(
    breaks = seq(-0.5, max(dga$`Richness, Alpha`,
                           dga$`Richness, Gamma`) + 1, 1)
  ) + ggplot2::geom_point(
    data = dga %>% group_by(
      Immigration, Extinction, Space, Dispersal
    ) %>% summarise(
      `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
      `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
      .groups = "drop"
    ),
    shape = 4, size = 3, color = "red", alpha = 0.8
  ) + ggplot2::scale_fill_viridis_c(
    direction = -1, trans = "log10"
  ) + ggplot2::facet_grid(
    Immigration + Extinction ~ Space + Dispersal,
    labeller = ggplot2::label_parsed
  ) + ggplot2::theme_bw(
  ) + ggplot2::geom_abline(
    slope = 1, intercept = 0
  ) + ggplot2::labs(
    title = if (exists("dTraj")) {
      "Sample trajectories in diversity counts"
    } else {
      "Gamma by Alpha"
    },
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm)
  )
  
  if (exists("dTraj")) {
    obj <- obj + ggplot2::geom_path(
      data = dTraj,
      mapping = ggplot2::aes(
        group = interaction(Environment, Number, History),
        color = Time
      ), alpha = 0.3
    ) + ggplot2::geom_point(
      data = dTrajStop,
      shape = 1, size = 3, color = "yellow", alpha = 0.8
    ) + ggplot2::scale_color_viridis_c(
      option = "B"
    )
  }
  
  return(obj)
}