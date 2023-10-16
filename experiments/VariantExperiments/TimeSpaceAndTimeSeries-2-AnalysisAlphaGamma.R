### Alpha:#####################################################################
# Note the order of operations is very relevant here.
# We imagine that the researchers consider all samples across times to be
# samples from the same population (even if it isn't in reality), so long as
# the control and experiment are respected.
# We then believe they would look at the difference between the control and the
# experiment.
# Edit: This reflects In\'es's idea that it is the difference between number of
# species between control and experiment.
bootstrapSamplesDeltaAlpha <- bootstrapSamples %>% dplyr::group_by(
  # Average across patches first.
  Type, Control, Bootstrap
) %>% dplyr::summarise(
  AverageSamplingAlpha = mean(SamplingAlpha),
  AverageSamplingAlphaNative = mean(SamplingAlphaNative),
  AverageSamplingAlphaInvasive = mean(SamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::group_by(
  # Then perform difference by converting control to negatives and adding.
  Type, Bootstrap
) %>% dplyr::mutate(
  dplyr::across(
    .cols = AverageSamplingAlpha : AverageSamplingAlphaInvasive,
    .fns = ~ ifelse(Control == "Control", -.x, .x)
  )
) %>% dplyr::summarise(
  DeltaAverageSamplingAlpha = sum(AverageSamplingAlpha),
  DeltaAverageSamplingAlphaNative = sum(AverageSamplingAlphaNative),
  DeltaAverageSamplingAlphaInvasive = sum(AverageSamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::rename(
  `Overall` = DeltaAverageSamplingAlpha,
  `Detected in Control` = DeltaAverageSamplingAlphaNative,
  `Not Detected in Control` = DeltaAverageSamplingAlphaInvasive
) %>% tidyr::pivot_longer(
  cols = c(`Overall`, `Detected in Control`, `Not Detected in Control`),
  names_to = "Species Subset",
  values_to = "Difference of Average Number of Species in Patch"
)

bootstrapSamplesPairedAlpha <- bootstrapSamples %>% dplyr::mutate(
  TimeSinceStart = round(TimeSinceStart, digits = 1) # Rare 1 != 1 issue.
) %>% dplyr::filter( # Excessive, but explicit:
  (Type == "Time series" & Control == "Control" & 
     TimeSinceStart == min(TimeSinceStart)) |
    (Type == "Time series" & Control != "Control" & 
       TimeSinceStart == max(TimeSinceStart)) |
    (Type != "Time series" & TimeSinceStart == max(TimeSinceStart))
) %>% dplyr::group_by(
  # Average across patches first.
  Type, Control, Bootstrap, TimeSinceStart
) %>% dplyr::mutate(
  DistanceFromCenter = dplyr::case_when(
    min(Patch) == 1 & max(Patch) == 10 ~ {
      ifelse(Patch < 5, Patch + 10, Patch) - median(
        ifelse(Patch < 5, Patch + 10, Patch)
      )
    },
    TRUE ~ Patch - median(Patch)
  ),
  DistanceFromCenterExpRev = ifelse(Control == "Experiment" & 
                                      Type == "Space for time", 
                                   -DistanceFromCenter,
                                   DistanceFromCenter)
) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  DistanceFromCenterExpRev, Bootstrap, TimeSinceStart, Type
) %>% dplyr::group_modify(
  .f = ~ computeSpeciesInControl(.x)
) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  DistanceFromCenterExpRev, Bootstrap, TimeSinceStart, Type, Control
) %>% dplyr::summarise(
  AverageSamplingAlpha = mean(SamplingAlpha),
  AverageSamplingAlphaNative = mean(SamplingAlphaNative),
  AverageSamplingAlphaInvasive = mean(SamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::group_by(
  # Then perform difference by converting control to negatives and adding.
  DistanceFromCenterExpRev, Bootstrap, Type
) %>% dplyr::mutate(
  dplyr::across(
    .cols = AverageSamplingAlpha : AverageSamplingAlphaInvasive,
    .fns = ~ ifelse(Control == "Control", -.x, .x)
  )
) %>% dplyr::summarise(
  DeltaAverageSamplingAlpha = sum(AverageSamplingAlpha),
  DeltaAverageSamplingAlphaNative = sum(AverageSamplingAlphaNative),
  DeltaAverageSamplingAlphaInvasive = sum(AverageSamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::rename(
  `Overall` = DeltaAverageSamplingAlpha,
  `Detected in Control` = DeltaAverageSamplingAlphaNative,
  `Not Detected in Control` = DeltaAverageSamplingAlphaInvasive
) %>% tidyr::pivot_longer(
  cols = c(`Overall`, `Detected in Control`, `Not Detected in Control`),
  names_to = "Species Subset",
  values_to = "Difference of Average Number of Species in Patch"
)

### Gamma: ####################################################################
# (Note: no guarantee of agreement due to differing sampling.)
bootstrapSamplesDeltaGamma <- bootstrapSamples %>% dplyr::group_by(
  # "Average" (Collect) across patches first.
  Type, Control, Bootstrap
) %>% dplyr::group_modify(
  .f = function(x, y) {
    IDsNative <- unique(unlist(strsplit(
      x$SamplingIDsNative, ", ", fixed = TRUE
    )))
    IDsInvasive <- unique(unlist(strsplit(
      x$SamplingIDsInvasive, ", ", fixed = TRUE
    )))
    
    data.frame(
      SamplingGammaNative = length(IDsNative),
      IDsNative = toString(IDsNative),
      SamplingGammaInvasive = length(IDsInvasive),
      IDsInvasive = toString(IDsInvasive)
    ) %>% dplyr::mutate(
      SamplingGamma = SamplingGammaNative + SamplingGammaInvasive
    )
  }
) %>% dplyr::group_by(
  # Then perform difference by converting control to negatives and adding.
  Type, Bootstrap
) %>% dplyr::mutate(
  dplyr::across(
    .cols = dplyr::starts_with("SamplingGamma"),
    .fns = ~ ifelse(Control == "Control", -.x, .x)
  )
) %>% dplyr::summarise(
  DeltaSamplingGamma = sum(SamplingGamma),
  DeltaSamplingGammaNative = sum(SamplingGammaNative),
  DeltaSamplingGammaInvasive = sum(SamplingGammaInvasive),
  .groups = "drop"
)  %>% dplyr::rename(
  `Overall` = DeltaSamplingGamma,
  `Detected in Control` = DeltaSamplingGammaNative,
  `Not Detected in Control` = DeltaSamplingGammaInvasive
) %>% tidyr::pivot_longer(
  cols = c(`Overall`, `Detected in Control`, `Not Detected in Control`),
  names_to = "Species Subset",
  values_to = "Difference of Number of Species across Patches"
)

# Plotting: ###################################################################
### Plot 1: Change in Local Richness: #########################################
plot_1_DeltaAlpha <- ggplot2::ggplot(
  bootstrapSamplesDeltaAlpha,
  ggplot2::aes(
    x = Type,
    y = `Difference of Average Number of Species in Patch`
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = 0.1,
  notch = TRUE
) + ggplot2::geom_line(
  data = bootstrapSamplesDeltaAlpha %>% dplyr::left_join(
    bootstrapSamplesDeltaAlpha %>% dplyr::group_by(
    Bootstrap
  ) %>% dplyr::filter(
    `Species Subset` == "Overall"
  ) %>% dplyr::arrange(Type) %>% dplyr::summarise(
    Slope = diff(`Difference of Average Number of Species in Patch`)
  ),
  by = "Bootstrap"
  ),
  ggplot2::aes(group = Bootstrap, color = Slope),
  alpha = 0.2
) + ggplot2::facet_wrap(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ .
) + ggplot2::labs(
  title = "Delta Alpha",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", file_result)
) + ggplot2::scale_color_viridis_c(option = "C")

# Note: Can try to establish if there is correlation here.
# with(
#   bootstrapSamplesDeltaAlpha %>% dplyr::filter(
#     `Species Subset` == "Overall"
#   ) %>% tidyr::pivot_wider(
#     names_from = "Type",
#     values_from = "Difference of Average Number of Species in Patch",
#     id_cols = "Bootstrap"
#   ), cor.test(`Space for time`, `Time series`)
# )

### Plot 2: ###################################################################
plot_2_DeltaGamma <- ggplot2::ggplot(
  bootstrapSamplesDeltaGamma,
  ggplot2::aes(
    x = Type,
    y = `Difference of Number of Species across Patches`
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = .1,
  notch = TRUE
) + ggplot2::geom_line(
  ggplot2::aes(group = Bootstrap),
  alpha = 0.1
) + ggplot2::facet_wrap(
  `Species Subset` ~ .
) + ggplot2::labs(
  title = paste0("Delta Gamma"),
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", file_result)
)

### Plot 3: ###################################################################
plot_1_PairedAlpha <- ggplot2::ggplot(
  bootstrapSamplesPairedAlpha,
  ggplot2::aes(
    x = Type,
    y = `Difference of Average Number of Species in Patch`
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = 0.1,
  notch = TRUE
) + ggplot2::geom_line(
  data = bootstrapSamplesPairedAlpha %>% dplyr::left_join(
    bootstrapSamplesPairedAlpha %>% dplyr::group_by(
      Bootstrap, DistanceFromCenterExpRev
    ) %>% dplyr::filter(
      `Species Subset` == "Overall"
    ) %>% dplyr::arrange(Type) %>% dplyr::summarise(
      Slope = diff(`Difference of Average Number of Species in Patch`)
    ),
    by = c("Bootstrap", "DistanceFromCenterExpRev")
  ),
  ggplot2::aes(group = Bootstrap, color = Slope),
  alpha = 0.2
) + ggplot2::facet_grid(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ DistanceFromCenterExpRev
) + ggplot2::labs(
  title = "Paired Alpha",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", file_result)
) + ggplot2::scale_color_viridis_c(option = "C")
