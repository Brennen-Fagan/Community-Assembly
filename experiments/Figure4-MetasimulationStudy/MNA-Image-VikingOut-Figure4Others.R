# Follows MNA-Image-VikingOut-Figure4.R
# Deprecated plots not used in the main text.

## Grouped by Extirpation: #####################################################
objs <- lapply(FacetPanels, function(Fac) {
  d <- concludingplotdata %>% dplyr::filter(
    FacetPanel == Fac
  ) %>% tidyr::separate(
    Neutral, sep = ", ", into = c("Imm.", "Ext.")
  ) %>% dplyr::mutate(
    Imm. = substr(Imm., 2, nchar(Imm.)),
    Ext. = substr(Ext., 1, nchar(Ext.) - 1)
  )
  fills <- unique(d$Area)

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Ext., y = Value,
                 fill = Area,
                 group = interaction(Area, Ext.))
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.8),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
  ) + ggplot2::scale_x_discrete(
    name = "Extirpation Rate Multipliers"
  ) + ggplot2::geom_line(
    data = d %>% dplyr::group_by(
      Area, FacetPanel, Ext.
    ) %>% dplyr::summarise(
      Value = mean(Value, na.rm = TRUE), .groups = "drop"
    ),
    mapping = ggplot2::aes(
      x = Ext., y = Value,
      group = interaction(Area),
      color = Area
      #linetype = Area
    ), inherit.aes = FALSE
  ) + ggplot2::scale_fill_manual(
    name = "", aesthetics = c("colour", "fill"),
    values = c(
      "Local Richness" = colorscheme[1],
      "Regional Richness" = colorscheme[2],
      "Temporal Jaccard" = colorscheme[3],
      "Spatial Jaccard" = colorscheme[4],
      "Median Invasibility" = colorscheme[5],
      "Average Invasibility" = colorscheme[5]
    ),
    breaks = fills
  ) + ggplot2::theme(
    # legend.position = #c(0.87, 0.93), Wide
    #   c(0.125, 0.225), # Square (half width, twice height) or single column
    # legend.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14)
  ) + ggplot2::ylab(
    if (Fac == FacetPanels[1]) {
      "Species\nRichness"
    } else if (Fac == FacetPanels[2]) {
      "Jaccard\nTurnover"
    } else if (Fac == FacetPanels[3]) {
      "End of Simulation\nInvasibility"
    } else {
      "Value"
    }
  )# + ggplot2::facet_grid(Imm. ~ Ext.)

  if (Fac != FacetPanels[3]) {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  } else {
  }

  obj
})

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "RichnessComparison-Invadability-ByExtirpation.png"
  ),
  plot = (
    objs[[1]]/objs[[2]]/objs[[3]] +
      patchwork::plot_layout(ncol = 1)#, heights = c(3, 3, 3, 0.5))#, 1, 1))
  ),
  width = plot_width, height = plot_height, units = plot_units,
  dpi = plot_dpi
)

## Grouped by Ext. + Dispersal: ################################################
objs <- lapply(FacetPanels, function(Fac) {
  d <- concludingplotdata %>% dplyr::filter(
    FacetPanel == Fac
  ) %>% tidyr::separate(
    Neutral, sep = ", ", into = c("Imm.", "Ext.")
  ) %>% dplyr::mutate(
    Imm. = substr(Imm., 2, nchar(Imm.)),
    Ext. = substr(Ext., 1, nchar(Ext.) - 1)
  )
  fills <- unique(d$Area)

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Ext., y = Value,
                 fill = Area,
                 group = interaction(Area, Ext.))
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.8),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
  ) + ggplot2::scale_x_discrete(
    name = "Extirpation Rate Multipliers"
  ) + ggplot2::geom_line(
    data = d %>% dplyr::group_by(
      Area, FacetPanel, Ext., Dispersal
    ) %>% dplyr::summarise(
      Value = mean(Value, na.rm = TRUE), .groups = "drop"
    ),
    mapping = ggplot2::aes(
      x = Ext., y = Value,
      group = interaction(Area, Dispersal),
      color = Area
      #linetype = Area
    ), inherit.aes = FALSE
  ) + ggplot2::scale_fill_manual(
    name = "", aesthetics = c("colour", "fill"),
    values = c(
      "Local Richness" = colorscheme[1],
      "Regional Richness" = colorscheme[2],
      "Temporal Jaccard" = colorscheme[3],
      "Spatial Jaccard" = colorscheme[4],
      "Median Invasibility" = colorscheme[5],
      "Average Invasibility" = colorscheme[5]
    ),
    breaks = fills
  ) + ggplot2::theme(
    legend.position = #c(0.87, 0.93), Wide
      c(0.85, 0.90), # Square (half width, twice height) or single column
    legend.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14)
  ) + ggplot2::ylab(
    if (Fac == FacetPanels[1]) {
      "Species\nRichness"
    } else if (Fac == FacetPanels[2]) {
      "Jaccard\nTurnover"
    } else if (Fac == FacetPanels[3]) {
      "End of Simulation\nInvasibility"
    } else {
      "Value"
    }
  ) + ggplot2::facet_grid(. ~ Dispersal)

  if (Fac != FacetPanels[3]) {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  } else {
  }

  obj
})

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "RichnessComparison-Invadability-ExtFacetDispersal.png"
  ),
  plot = (
    objs[[1]]/objs[[2]]/objs[[3]] +
      patchwork::plot_layout(ncol = 1)#, heights = c(3, 3, 3, 0.5))#, 1, 1))
  ),
  width = plot_width, height = plot_height, units = plot_units,
  dpi = plot_dpi
)

## Grouped by Imm. + Dispersal: ################################################
objs <- lapply(FacetPanels, function(Fac) {
  d <- concludingplotdata %>% dplyr::filter(
    FacetPanel == Fac
  ) %>% tidyr::separate(
    Neutral, sep = ", ", into = c("Imm.", "Ext.")
  ) %>% dplyr::mutate(
    Imm. = substr(Imm., 2, nchar(Imm.)),
    Ext. = substr(Ext., 1, nchar(Ext.) - 1)
  )
  fills <- unique(d$Area)

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Imm., y = Value,
                 fill = Area,
                 group = interaction(Area, Imm.))
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.8),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
  ) + ggplot2::scale_x_discrete(
    name = "Immigration Rate Multipliers"
  ) + ggplot2::geom_line(
    data = d %>% dplyr::group_by(
      Area, FacetPanel, Imm., Dispersal
    ) %>% dplyr::summarise(
      Value = mean(Value, na.rm = TRUE), .groups = "drop"
    ),
    mapping = ggplot2::aes(
      x = Imm., y = Value,
      group = interaction(Area, Dispersal),
      color = Area
      #linetype = Area
    ), inherit.aes = FALSE
  ) + ggplot2::scale_fill_manual(
    name = "", aesthetics = c("colour", "fill"),
    values = c(
      "Local Richness" = colorscheme[1],
      "Regional Richness" = colorscheme[2],
      "Temporal Jaccard" = colorscheme[3],
      "Spatial Jaccard" = colorscheme[4],
      "Median Invasibility" = colorscheme[5],
      "Average Invasibility" = colorscheme[5]
    ),
    breaks = fills
  ) + ggplot2::theme(
    legend.position = #c(0.87, 0.93), Wide
      c(0.85, 0.90), # Square (half width, twice height) or single column
    legend.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14)
  ) + ggplot2::ylab(
    if (Fac == FacetPanels[1]) {
      "Species\nRichness"
    } else if (Fac == FacetPanels[2]) {
      "Jaccard\nTurnover"
    } else if (Fac == FacetPanels[3]) {
      "End of Simulation\nInvasibility"
    } else {
      "Value"
    }
  ) + ggplot2::facet_grid(. ~ Dispersal)

  if (Fac != FacetPanels[3]) {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  } else {
  }

  obj
})

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "RichnessComparison-Invadability-ImmFacetDispersal.png"
  ),
  plot = (
    objs[[1]]/objs[[2]]/objs[[3]] +
      patchwork::plot_layout(ncol = 1)#, heights = c(3, 3, 3, 0.5))#, 1, 1))
  ),
  width = plot_width, height = plot_height, units = plot_units,
  dpi = plot_dpi
)

## Grouped by Pool & Noise + Dispersal: ########################################

objs <- lapply(FacetPanels, function(Fac) {
  d <- concludingplotdata %>% dplyr::filter(
    FacetPanel == Fac
  ) %>% tidyr::unite(
    Pool, Noise, col = "Pool & Noise", sep = " & "
  )

  fills <- unique(d$Area)

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = `Pool & Noise`, y = Value,
                 fill = Area,
                 group = interaction(Area, `Pool & Noise`))
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.8),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
  ) + ggplot2::scale_x_discrete(
    name = "Pool & Noise"
  ) + ggplot2::geom_line(
    data = d %>% dplyr::group_by(
      Area, FacetPanel, `Pool & Noise`, Dispersal
    ) %>% dplyr::summarise(
      Value = mean(Value, na.rm = TRUE), .groups = "drop"
    ),
    mapping = ggplot2::aes(
      x = `Pool & Noise`, y = Value,
      group = interaction(Area, Dispersal),
      color = Area
      #linetype = Area
    ), inherit.aes = FALSE
  ) + ggplot2::scale_fill_manual(
    name = "", aesthetics = c("colour", "fill"),
    values = c(
      "Local Richness" = colorscheme[1],
      "Regional Richness" = colorscheme[2],
      "Temporal Jaccard" = colorscheme[3],
      "Spatial Jaccard" = colorscheme[4],
      "Median Invasibility" = colorscheme[5],
      "Average Invasibility" = colorscheme[5]
    ),
    breaks = fills
  ) + ggplot2::theme(
    legend.position = #c(0.87, 0.93), Wide
      c(0.85, 0.90), # Square (half width, twice height) or single column
    legend.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14)
  ) + ggplot2::ylab(
    if (Fac == FacetPanels[1]) {
      "Species\nRichness"
    } else if (Fac == FacetPanels[2]) {
      "Jaccard\nTurnover"
    } else if (Fac == FacetPanels[3]) {
      "End of Simulation\nInvasibility"
    } else {
      "Value"
    }
  ) + ggplot2::facet_grid(. ~ Dispersal)

  if (Fac != FacetPanels[3]) {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  } else {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1)
    )
  }

  obj
})

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "RichnessComparison-Invadability-PoolNoiseFacetDispersal.png"
  ),
  plot = (
    objs[[1]]/objs[[2]]/objs[[3]] +
      patchwork::plot_layout(ncol = 1)#, heights = c(3, 3, 3, 0.5))#, 1, 1))
  ),
  width = plot_width, height = plot_height, units = plot_units,
  dpi = plot_dpi
)
