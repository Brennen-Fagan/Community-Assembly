# Follows MNA-Image-VikingOut-Figure4.R

forHeatMaps <- dplyr::bind_rows(
  Attributes %>% dplyr::mutate(Source = "Div"),
  AttrInv %>% dplyr::mutate(Source = "Inv"),
  AttrTimeJacs %>% dplyr::mutate(Source = "Jac")
) %>% dplyr::group_by(
  Pool, Noise, Neutral, Space, Set, CaseNumber, History, Part
) %>% dplyr::summarise(
  sources = toString(unique(Source)),
  count = dplyr::n()
) %>% dplyr::filter(
  sources == "Div, Inv, Jac", count == 3
) %>% dplyr::mutate(
  SetPaper = dplyr::case_when(
    Pool %in% "0, 0, 0, 0" &
      Noise %in% "1" &
      Neutral %in% c("1, 0", "0.1, 0.1", "10, 0.1", "0.1, 10", "10, 10") ~ 1,
    Pool %in% "0, 0, 0, 0" &
      Noise %in% "0" &
      Neutral %in% c("1, 0", "1, 1") ~ 2,
    Pool %in% c("0, 0, 0, 0", "-1, 0, 0, 0", "0, 0, 0, 1", "-1, 0, 0, 1") &
      Noise %in% "1" &
      Neutral %in% c("1, 1", "10, 1", "0.1, 1", "1, 10", "1, 0.1") ~ 3,
    TRUE ~ 4
  )
) %>% tidyr::separate(
  Neutral, into = c("Immigration", "Extirpation"), sep = ", "
) %>% tidyr::unite(
  "PoolNoise", Pool, Noise, sep = " & "
)

forHeatMaps <- forHeatMaps %>% dplyr::left_join(labelsforplot, by = "Space")

((
  (
    ggplot2::ggplot(
      forHeatMaps,
      ggplot2::aes(x = Dispersal, y = Immigration)
    ) + ggplot2::theme_bw() + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
  ) +
    (ggplot2::ggplot(
      forHeatMaps,
      ggplot2::aes(x = Dispersal, y = Extirpation)
    ) + ggplot2::theme_bw() + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
    )
) / (
  (
    ggplot2::ggplot(
      forHeatMaps,
      ggplot2::aes(x = Dispersal, y = `PoolNoise`)
    ) + ggplot2::ylab(
      "Pool & Noise"
    ) + ggplot2::theme_bw() + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
  ) +
    ggplot2::ggplot(
      forHeatMaps,
      ggplot2::aes(x = Immigration, y = Extirpation)
    ) + ggplot2::theme_bw()
) / (
  ggplot2::ggplot(
    forHeatMaps,
    ggplot2::aes(x = Immigration, y = `PoolNoise`)
  ) + ggplot2::ylab(
    "Pool & Noise"
  ) + ggplot2::theme_bw() + ggplot2::ggplot(
    forHeatMaps,
    ggplot2::aes(x = Extirpation, y = `PoolNoise`)
  ) + ggplot2::ylab(
    "Pool & Noise"
  ) + ggplot2::theme_bw()
)) & ggplot2::geom_bin2d(
) & ggplot2::scale_fill_viridis_c(
  trans = "log10", direction = -1
)

ggplot2::ggsave(
  filename = file.path(outputLocation, "SimsCompleted-Absolute.png"),
  width = plot_width, height = plot_height, units = plot_units,
  dpi = plot_dpi
)

forHeatMaps_Desired <- dplyr::bind_rows(
  expand_grid(
    PoolNoise = c("0, 0, 0, 0 & 1"),
    Neutral = c("1, 0", "0.1, 0.1", "10, 0.1", "0.1, 10", "10, 10"),
    Space = 10^c(Inf, 9:3, 0)
  ),
  expand_grid(
    PoolNoise = c("0, 0, 0, 0 & 0"),
    Neutral = c("1, 0", "1, 1"),
    Space = 10^c(Inf, 9:3, 0)
  ),
  expand_grid(
    PoolNoise = c("0, 0, 0, 0 & 1", "-1, 0, 0, 0 & 1", "0, 0, 0, 1 & 1", "-1, 0, 0, 1 & 1"),
    Neutral = c("1, 1", "1, 0.1", "1, 10", "0.1, 1", "10, 1"),
    Space = 10^c(Inf, 9:3, 0)
  ),
  expand_grid(
    PoolNoise = c("0, 0, 0, 0 & 1", "-1, 0, 0, 0 & 1", "0, 0, 0, 1 & 1", "-1, 0, 0, 1 & 1", "0, 0, 0, 0 & 0"),
    Neutral = c("1, 1", "1, 0.1", "1, 10", "0.1, 1", "10, 1",
                "1, 0", "0.1, 0.1", "10, 0.1", "0.1, 10", "10, 10"),
    Space = 10^c(2, 1)
  )
) %>% dplyr::rename(
  `PoolNoise` = PoolNoise
) %>% tidyr::separate(
  Neutral, sep = ", ", into = c("Immigration", "Extirpation")
) %>% dplyr::mutate(
  Space = as.character(Space)
) %>% dplyr::left_join(labelsforplot, by = "Space")

forHeatMaps_Relative <- forHeatMaps %>% dplyr::group_by(
  `PoolNoise`, Immigration, Extirpation, Space, Dispersal
) %>% dplyr::summarise(
  Count = dplyr::n()
) %>% dplyr::full_join(
  forHeatMaps_Desired %>% dplyr::mutate(Desired = 100),
  by = c("PoolNoise", "Immigration", "Extirpation", "Space", "Dispersal")
) %>% dplyr::mutate(
  Count = ifelse(is.na(Count), 0, Count)
)

heatMaps_Relative <- vector("list", 6)
heatMaps_Relative_choices <- c(
  "Dispersal", "Immigration", "Extirpation", "PoolNoise"
)
heatMaps_Relative_counter <- 1

for(x in 1:3) {
  for (y in (x + 1):4) {
    heatMaps_Relative[[heatMaps_Relative_counter]] <- ggplot2::ggplot(
      forHeatMaps_Relative %>% dplyr::group_by(
        dplyr::across(dplyr::all_of(c(
          heatMaps_Relative_choices[x],
          heatMaps_Relative_choices[y]
        )))
      ) %>% dplyr::summarise(
        Percentage = sum(Count) / sum(Desired) # Scales does * 100
      ),
      mapping = ggplot2::aes_string(
        x = heatMaps_Relative_choices[x],
        y = heatMaps_Relative_choices[y],
        fill = "Percentage"
      )
    ) + ggplot2::theme_bw(
    )

    if(heatMaps_Relative_choices[x] == "Dispersal")
      heatMaps_Relative[[heatMaps_Relative_counter]] <-
        heatMaps_Relative[[heatMaps_Relative_counter]] + ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)
        )
    if(heatMaps_Relative_choices[y] == "PoolNoise")
      heatMaps_Relative[[heatMaps_Relative_counter]] <-
        heatMaps_Relative[[heatMaps_Relative_counter]] + ggplot2::ylab(
          "Pool & Noise"
        )
    heatMaps_Relative_counter <- heatMaps_Relative_counter + 1
  }
}

(
  ((
    heatMaps_Relative[[1]] +
      heatMaps_Relative[[2]]
  ) / (
    heatMaps_Relative[[3]] +
      heatMaps_Relative[[4]]
  ) / (
    heatMaps_Relative[[5]] +
      heatMaps_Relative[[6]]
  )) & ggplot2::geom_raster(
  ) & ggplot2::scale_fill_viridis_c(
    direction = -1,
    limits = c(0, 1),
    labels = scales::percent_format()
  )
) + patchwork::plot_layout(guides = "collect")

ggplot2::ggsave(
  filename = file.path(outputLocation, "SimsCompleted-Relative.png"),
  width = plot_width, height = plot_height, units = plot_units,
  dpi = plot_dpi
)
