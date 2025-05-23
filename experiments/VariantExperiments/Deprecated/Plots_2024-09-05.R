options(bitmapType = "cairo")
imageHeight <- 2400; imageWidth <- 3600; imageUnits <- "px"
saveImages <- FALSE
dispersal <- "15" # "NA", "10", "15"

# 0, 0.5, 1's with Intensity 2 ################################################
dplyr::bind_rows(
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c(
      "(0, 0)", "(0.5, 0.5)", "(1, 1)", "(0.5, 0.5)->(0.5, 1)", 
      "(0.5, 0.5)->(0.5, 0)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)"
    ),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 2
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = dplyr::case_when(
      Intervention %in% c("(0.5, 0.5)", "(0.5, 0.5)->(0.5, 1)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0, 0)", "(1, 1)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)") ~ 2,
      TRUE ~ 3
    ),
    FacetCol = dplyr::case_when(
      Intervention %in% c("(0, 0)", "(0, 0)->(0, 0.5)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0.5, 0.5)->(0.5, 1)", "(1, 1)",  "(1, 1)->(1, 0.5)") ~ 2,
      Intervention %in% c("(0.5, 0.5)") ~ 1
    ),
    Alpha = ifelse(
      PoolPatchAffinitySeed == min(PoolPatchAffinitySeed) & EventsSeed == min(EventsSeed),
      0.25, 1
    )
  ),
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c("(0.5, 0.5)"),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 2
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = 1,
    FacetCol = 2,
    Alpha = ifelse(
      PoolPatchAffinitySeed == min(PoolPatchAffinitySeed) & EventsSeed == min(EventsSeed),
      0.25, 1
    )
  )
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0, 0)" = "black",
    "(0.5, 0.5)" = "black",
    "(1, 1)" = "black",
    "(0.5, 0.5)->(0.5, 1)" = "pink", 
    "(0.5, 0.5)->(0.5, 0)" = "purple",
    "(1, 1)->(1, 0.5)" = "red",
    "(0, 0)->(0, 0.5)" = "blue"
  )
)
if (saveImages)
ggplot2::ggsave(
  filename = "Plot_Alpha_InterventionsAround05_Intens2.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

# 0, 0.5, 1's with Intensity 3 ################################################
dplyr::bind_rows(
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c(
      "(0, 0)", "(0.5, 0.5)", "(1, 1)", "(0.5, 0.5)->(0.5, 1)", 
      "(0.5, 0.5)->(0.5, 0)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)"
    ),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 3
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = dplyr::case_when(
      Intervention %in% c("(0.5, 0.5)", "(0.5, 0.5)->(0.5, 1)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0, 0)", "(1, 1)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)") ~ 2,
      TRUE ~ 3
    ),
    FacetCol = dplyr::case_when(
      Intervention %in% c("(0, 0)", "(0, 0)->(0, 0.5)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0.5, 0.5)->(0.5, 1)", "(1, 1)",  "(1, 1)->(1, 0.5)") ~ 2,
      Intervention %in% c("(0.5, 0.5)") ~ 1
    ),
    Alpha = ifelse(
      PoolPatchAffinitySeed == min(PoolPatchAffinitySeed) & EventsSeed == min(EventsSeed),
      1, 0
    )
  ),
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c("(0.5, 0.5)"),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 3
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = 1,
    FacetCol = 2,
    Alpha = ifelse(
      PoolPatchAffinitySeed == min(PoolPatchAffinitySeed) & EventsSeed == min(EventsSeed),
      1, 0
    )
  )
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0, 0)" = "black",
    "(0.5, 0.5)" = "black",
    "(1, 1)" = "black",
    "(0.5, 0.5)->(0.5, 1)" = "pink", 
    "(0.5, 0.5)->(0.5, 0)" = "purple",
    "(1, 1)->(1, 0.5)" = "red",
    "(0, 0)->(0, 0.5)" = "blue"
  )
)

if (saveImages)
ggplot2::ggsave(
  filename = "Plot_Alpha_InterventionsAround05_Intens3.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

# 0.25, 0.75's with Intensity 2 ###############################################

diversitiesRounded %>% dplyr::filter(
  is.na(`Species Layer`), 
  Intervention %in% c(
    "(0.25, 0.25)", "(0.75, 0.75)",
    "(0.25, 0.25)->(0.25, 0.75)", "(0.75, 0.75)->(0.75, 0.25)"
  ),
  Alignment == "All Species", 
  Dispersal == dispersal, 
  Measurement == "Alpha Richness", 
  Environment == 2,
  NicheDistance == 2
) %>% dplyr::ungroup() %>% dplyr::group_by(
  PoolPatchAffinity,
  InterventionPatchType # NA or not
) %>% dplyr::mutate(
  FacetRow = 1,
  FacetCol = dplyr::case_when(
    Intervention %in% c("(0.25, 0.25)", "(0.25, 0.25)->(0.25, 0.75)") ~ 1,
    Intervention %in% c("(0.75, 0.75)", "(0.75, 0.75)->(0.75, 0.25)") ~ 2,
    TRUE ~ 3
  ),
  Alpha = ifelse(
    PoolPatchAffinitySeed == min(PoolPatchAffinitySeed) & EventsSeed == min(EventsSeed),
    1, 0
  )
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0.25, 0.25)" = "black",
    "(0.75, 0.75)" = "black",
    "(0.25, 0.25)->(0.25, 0.75)" = "green", 
    "(0.75, 0.75)->(0.75, 0.25)" = "orange"
  )
)

if (saveImages)
ggplot2::ggsave(
  filename = "Plot_Alpha_SymmetriesAround05_Intens2.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

# 0.25, 0.75's with Intensity 3 ###############################################

diversitiesRounded %>% dplyr::filter(
  is.na(`Species Layer`), 
  Intervention %in% c(
    "(0.25, 0.25)", "(0.75, 0.75)",
    "(0.25, 0.25)->(0.25, 0.75)", "(0.75, 0.75)->(0.75, 0.25)"
  ),
  Alignment == "All Species", 
  Dispersal == dispersal, 
  Measurement == "Alpha Richness", 
  Environment == 2,
  NicheDistance == 3
) %>% dplyr::ungroup() %>% dplyr::group_by(
  PoolPatchAffinity,
  InterventionPatchType # NA or not
) %>% dplyr::mutate(
  FacetRow = 1,
  FacetCol = dplyr::case_when(
    Intervention %in% c("(0.25, 0.25)", "(0.25, 0.25)->(0.25, 0.75)") ~ 1,
    Intervention %in% c("(0.75, 0.75)", "(0.75, 0.75)->(0.75, 0.25)") ~ 2,
    TRUE ~ 3
  ),
  Alpha = ifelse(
    PoolPatchAffinitySeed == min(PoolPatchAffinitySeed) & EventsSeed == min(EventsSeed),
    1, 0
  )
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0.25, 0.25)" = "black",
    "(0.75, 0.75)" = "black",
    "(0.25, 0.25)->(0.25, 0.75)" = "green", 
    "(0.75, 0.75)->(0.75, 0.25)" = "orange"
  )
)

if (saveImages)
ggplot2::ggsave(
  filename = "Plot_Alpha_SymmetriesAround05_Intens3.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

# 0, 0.5, 1's with Intensity 2, Means #########################################
dplyr::bind_rows(
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c(
      "(0, 0)", "(0.5, 0.5)", "(1, 1)", "(0.5, 0.5)->(0.5, 1)", 
      "(0.5, 0.5)->(0.5, 0)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)"
    ),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 2
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity, Time,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = dplyr::case_when(
      Intervention %in% c("(0.5, 0.5)", "(0.5, 0.5)->(0.5, 1)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0, 0)", "(1, 1)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)") ~ 2,
      TRUE ~ 3
    ),
    FacetCol = dplyr::case_when(
      Intervention %in% c("(0, 0)", "(0, 0)->(0, 0.5)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0.5, 0.5)->(0.5, 1)", "(1, 1)",  "(1, 1)->(1, 0.5)") ~ 2,
      Intervention %in% c("(0.5, 0.5)") ~ 1
    ),
    Value = mean(Value),
    Alpha = 1
  ),
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c("(0.5, 0.5)"),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 2
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity, Time,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = 1,
    FacetCol = 2,
    Value = mean(Value),
    Alpha = 1
  )
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0, 0)" = "black",
    "(0.5, 0.5)" = "black",
    "(1, 1)" = "black",
    "(0.5, 0.5)->(0.5, 1)" = "pink", 
    "(0.5, 0.5)->(0.5, 0)" = "purple",
    "(1, 1)->(1, 0.5)" = "red",
    "(0, 0)->(0, 0.5)" = "blue"
  )
)

if (saveImages)
ggplot2::ggsave(
  filename = "Plot_AlphaMean_InterventionsAround05_Intens2.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

# 0, 0.5, 1's with Intensity 3, Means #########################################
dplyr::bind_rows(
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c(
      "(0, 0)", "(0.5, 0.5)", "(1, 1)", "(0.5, 0.5)->(0.5, 1)", 
      "(0.5, 0.5)->(0.5, 0)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)"
    ),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 3
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity, Time,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = dplyr::case_when(
      Intervention %in% c("(0.5, 0.5)", "(0.5, 0.5)->(0.5, 1)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0, 0)", "(1, 1)", "(0, 0)->(0, 0.5)", "(1, 1)->(1, 0.5)") ~ 2,
      TRUE ~ 3
    ),
    FacetCol = dplyr::case_when(
      Intervention %in% c("(0, 0)", "(0, 0)->(0, 0.5)", "(0.5, 0.5)->(0.5, 0)") ~ 1,
      Intervention %in% c("(0.5, 0.5)->(0.5, 1)", "(1, 1)",  "(1, 1)->(1, 0.5)") ~ 2,
      Intervention %in% c("(0.5, 0.5)") ~ 1
    ),
    Value = mean(Value),
    Alpha = 1
  ),
  diversitiesRounded %>% dplyr::filter(
    is.na(`Species Layer`), 
    Intervention %in% c("(0.5, 0.5)"),
    Alignment == "All Species", 
    Dispersal == dispersal, 
    Measurement == "Alpha Richness", 
    Environment == 2,
    NicheDistance == 3
  ) %>% dplyr::ungroup() %>% dplyr::group_by(
    PoolPatchAffinity, Time,
    InterventionPatchType # NA or not
  ) %>% dplyr::mutate(
    FacetRow = 1,
    FacetCol = 2,
    Value = mean(Value),
    Alpha = 1
  )
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0, 0)" = "black",
    "(0.5, 0.5)" = "black",
    "(1, 1)" = "black",
    "(0.5, 0.5)->(0.5, 1)" = "pink", 
    "(0.5, 0.5)->(0.5, 0)" = "purple",
    "(1, 1)->(1, 0.5)" = "red",
    "(0, 0)->(0, 0.5)" = "blue"
  )
)

if (saveImages)
ggplot2::ggsave(
  filename = "Plot_AlphaMean_InterventionsAround05_Intens3.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

# 0.25, 0.75's with Intensity 2, Means #########################################

diversitiesRounded %>% dplyr::filter(
  is.na(`Species Layer`), 
  Intervention %in% c(
    "(0.25, 0.25)", "(0.75, 0.75)",
    "(0.25, 0.25)->(0.25, 0.75)", "(0.75, 0.75)->(0.75, 0.25)"
  ),
  Alignment == "All Species", 
  Dispersal == dispersal, 
  Measurement == "Alpha Richness", 
  Environment == 2,
  NicheDistance == 2
) %>% dplyr::ungroup() %>% dplyr::group_by(
  PoolPatchAffinity, Time,
  InterventionPatchType # NA or not
) %>% dplyr::mutate(
  FacetRow = 1,
  FacetCol = dplyr::case_when(
    Intervention %in% c("(0.25, 0.25)", "(0.25, 0.25)->(0.25, 0.75)") ~ 1,
    Intervention %in% c("(0.75, 0.75)", "(0.75, 0.75)->(0.75, 0.25)") ~ 2,
    TRUE ~ 3
  ),
  Value = mean(Value),
  Alpha = 1
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0.25, 0.25)" = "black",
    "(0.75, 0.75)" = "black",
    "(0.25, 0.25)->(0.25, 0.75)" = "green", 
    "(0.75, 0.75)->(0.75, 0.25)" = "orange"
  )
)

if (saveImages)
ggplot2::ggsave(
  filename = "Plot_AlphaMean_SymmetriesAround05_Intens2.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

# 0.25, 0.75's with Intensity 3, Means #########################################

diversitiesRounded %>% dplyr::filter(
  is.na(`Species Layer`), 
  Intervention %in% c(
    "(0.25, 0.25)", "(0.75, 0.75)",
    "(0.25, 0.25)->(0.25, 0.75)", "(0.75, 0.75)->(0.75, 0.25)"
  ),
  Alignment == "All Species", 
  Dispersal == dispersal, 
  Measurement == "Alpha Richness", 
  Environment == 2,
  NicheDistance == 3
) %>% dplyr::ungroup() %>% dplyr::group_by(
  PoolPatchAffinity, Time,
  InterventionPatchType # NA or not
) %>% dplyr::mutate(
  FacetRow = 1,
  FacetCol = dplyr::case_when(
    Intervention %in% c("(0.25, 0.25)", "(0.25, 0.25)->(0.25, 0.75)") ~ 1,
    Intervention %in% c("(0.75, 0.75)", "(0.75, 0.75)->(0.75, 0.25)") ~ 2,
    TRUE ~ 3
  ),
  Value = mean(Value),
  Alpha = 1
) %>% plotDiversityOverview() + ggplot2::facet_grid(
  ggplot2::vars(FacetRow), ggplot2::vars(FacetCol)
) + ggplot2::scale_color_manual(
  values = c(
    "(0.25, 0.25)" = "black",
    "(0.75, 0.75)" = "black",
    "(0.25, 0.25)->(0.25, 0.75)" = "green", 
    "(0.75, 0.75)->(0.75, 0.25)" = "orange"
  )
)

if (saveImages)
ggplot2::ggsave(
  filename = "Plot_AlphaMean_SymmetriesAround05_Intens3.png",
  height = imageHeight, width = imageWidth, units = imageUnits
)

################################################################################
