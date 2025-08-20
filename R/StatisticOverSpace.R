StatisticOverSpace <- function(
  i, d, levels,
  spacelabs,
  # ytarget, SET TO Richness
  ylabel, title,
  burnin = time_burnin,
  burnout = time_burnout
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]
  
  # Attach labels.
  d <- d %>% dplyr::filter(
    pbs == Pool,
    inm == Noise,
    Time > burnin, Time < burnout
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = "")
  ) %>% tidyr::unite(
    col = "group", Neutral, Space, remove = FALSE
  )
  
  xlabs <- d %>% dplyr::select(
    Space, Dispersal
  ) %>% dplyr::distinct() %>% dplyr::arrange(Space) %>% dplyr::mutate(
    Space = as.character(Space),
    Label = paste(Space, Dispersal, sep = "\n")
  )
  
  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = as.character(Space), y = Richness,
                 fill = Neutral,
                 group = group)
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.5), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.5),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
    title = title,
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm),
    y = ylabel
  ) + ggplot2::scale_x_discrete(
    name = "Inter-community distance and unit proportion dispersing",
    breaks = xlabs$Space,
    labels = xlabs$Label
  ) + ggplot2::scale_fill_discrete(
    name = "(Imm., Ext.)\nModifiers"
  )
  
  return(obj)
}