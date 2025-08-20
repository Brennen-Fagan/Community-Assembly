StatisticOverTime <- function(
  i, d, da, levels,
  spacelabs, neutrallabs,
  ytarget, ylabel, title,
  burnin = time_burnin,
  burnout = time_burnout,
  units = time_units
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]
  
  # Attach labels.
  d <- d %>% dplyr::filter(
    pbs == Pool,
    inm == Noise
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  ) %>% tidyr::unite(
    col = "group", !Time & !dplyr::all_of(ytarget), remove = FALSE
  )
  
  da <- da %>% dplyr::filter(
    pbs == Pool,
    inm == Noise
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  ) %>% tidyr::unite(
    col = "group", !Time & !Low & !High & !Central, remove = FALSE
  )
  
  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Time,
                 group = group)
  ) + ggplot2::geom_line(
    mapping = ggplot2::aes_string(
      y = ytarget
    ),
    alpha = 0.05
  ) + ggplot2::geom_line(
    data = da,
    mapping = ggplot2::aes(
      y = Central
    ),
    color = "red", size = 1.5
  ) + ggplot2::geom_ribbon(
    data = da,
    mapping = ggplot2::aes(
      ymin = Low,
      ymax = High,
      x = Time
    ),
    alpha = 0.4,
    fill = "blue",
    inherit.aes = FALSE
  ) + ggplot2::geom_vline(
    xintercept = c(burnin, burnout)
  ) + ggplot2::theme_bw(
  ) + ggplot2::facet_grid(
    Immigration + Extinction ~ Space + Dispersal,
    labeller = ggplot2::label_parsed
  ) + ggplot2::labs(
    title = title,
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm),
    ylab = ylabel
  ) + ggplot2::xlab(
    bquote("Time, " * 10^.(log10(time_units)) * " Units")
  )
  
  return(obj)
}