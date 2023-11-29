# Fix naming differences and bootstraps vs interventions differences.
if (exists("result")) {
  preferredTimeGap <- which.max(timediffs >= result$ReactionTime*10)
  result_ <- result
  files_ <- file_result
  timeColumn <- "Time"
} else if (exists("resultBase")) {
  preferredTimeGap <- which.max(timediffs >= 10) # resultBase's time rescaled.
  result_ <- resultBase
  files_ <- paste(fileBaseResult, fileInterventionResult, collapse = ", ")
  timeColumn <- "TimeBase"
}

ylimRichnessChange <- c(-10, 10)

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
  TimeSinceStart = round(# Rare 1 != 1 issue.
    TimeSinceStart,
    digits = ceiling(-log10(min(diff(result_$Abundance[, 1])))))
) %>% dplyr::group_by(
  Type, Control, Bootstrap, Patch
) %>% dplyr::arrange(
  TimeSinceStart
) %>% dplyr::mutate(
  TimeGapNumber = seq_along(TimeSinceStart),
  TimeGapNumber = ifelse(
    Type == "Time series" & Control == "Control",
    rev(TimeGapNumber), TimeGapNumber
  ) # i.e.  5, 4, 3, 2, 1, ___, 1, 2, 3, 4, 5
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
  DistanceFromCenterExpRev, Bootstrap, TimeGapNumber, Type
) %>% dplyr::group_modify(
  .f = ~ computeSpeciesInControl(.x, Time = timeColumn)
) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  DistanceFromCenterExpRev, Bootstrap, TimeGapNumber, Type, Control
) %>% dplyr::summarise( # Effectively a rename, but useful struct. for later.
  AverageSamplingAlpha = mean(SamplingAlpha),
  AverageSamplingAlphaNative = mean(SamplingAlphaNative),
  AverageSamplingAlphaInvasive = mean(SamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::group_by(
  # Then perform difference by converting control to negatives and adding.
  DistanceFromCenterExpRev, Bootstrap, TimeGapNumber, Type
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

bootstrapSamplesTimedAlpha <- bootstrapSamples %>% dplyr::mutate(
  TimeSinceStart = round(# Rare 1 != 1 issue.
    TimeSinceStart,
    digits = ceiling(-log10(min(diff(result_$Abundance[, 1])))))
) %>% dplyr::group_by(
  Type, Control, Bootstrap, Patch
) %>% dplyr::arrange(
  TimeSinceStart
) %>% dplyr::mutate(
  TimeGapNumber = seq_along(TimeSinceStart),
  TimeGapNumber = ifelse(
    Type == "Time series" & Control == "Control",
    rev(TimeGapNumber), TimeGapNumber
  ) # i.e.  5, 4, 3, 2, 1, ___, 1, 2, 3, 4, 5
) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  Bootstrap, TimeGapNumber, Type
) %>% dplyr::group_modify(
  .f = ~ computeSpeciesInControl(.x, Time = timeColumn)
) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  Bootstrap, TimeGapNumber, Type, Control
) %>% dplyr::summarise( # Average over Space.
  AverageSamplingAlpha = mean(SamplingAlpha),
  AverageSamplingAlphaNative = mean(SamplingAlphaNative),
  AverageSamplingAlphaInvasive = mean(SamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::group_by(
  # Then perform difference by converting control to negatives and adding.
  Bootstrap, TimeGapNumber, Type
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

### Beta: #####################################################################
# Break into pairs across time and space that test the effect of temporal and
# spatial distance.
bootstrapSamplesPairedBeta_source <- bootstrapSamples %>% addDistanceColumns(
  mindig = ceiling(-log10(min(diff(result_$Abundance[, 1])))),
  Time = timeColumn
) %>% dplyr::mutate(
  Patch = ifelse(
    Control == "Experiment",
    paste0(Patch, "E"),
    paste0(Patch, "C")
  )
)

##### Perform Beta for all, control, !control: ################################
# Note we want to do this twice:
#    once for fixed patches and
#    once for fixed times.

# Note the ", " combination is reserved and required in ConvertPreparedToBeta.
postfixes <- c(", All", ", Det. In Control", ", Not Det. In Control")
betacolumns <- c("SamplingIDs", "SamplingIDsNative", "SamplingIDsInvasive")

bootstrapSamplesPairedBeta_Jaccard <- ConvertPreparedToBeta(
  bootstrapSamplesPairedBeta_source %>% dplyr::filter(
    TimeGapNumber == preferredTimeGap
  ) %>% dplyr::group_by(
    DistanceFromCenterExpRev, Bootstrap, Type
  ),
  columns = betacolumns,
  method = "Jaccard", presenceabsence = TRUE,
  postfixes = postfixes,
  indicator = "DistanceFromCenterExpRev"
)

bootstrapSamplesPairedBeta_BrayCurtis <- ConvertPreparedToBeta(
  bootstrapSamplesPairedBeta_source %>% dplyr::filter(
    TimeGapNumber == preferredTimeGap
  ) %>% dplyr::group_by(
    DistanceFromCenterExpRev, Bootstrap, Type
  ),
  columns = betacolumns,
  method = "Bray", presenceabsence = FALSE,
  postfixes = postfixes,
  indicator = "DistanceFromCenterExpRev"
)

bootstrapSamplesTimedBeta_Jaccard <- ConvertPreparedToBeta(
  bootstrapSamplesPairedBeta_source %>% dplyr::filter(
    DistanceFromCenterExpRev == 0
  ) %>% dplyr::group_by(
    TimeGapNumber, Bootstrap, Type
  ),
  columns = betacolumns,
  method = "Jaccard", presenceabsence = TRUE,
  postfixes = postfixes,
  indicator = "TimeGapNumber"
)

bootstrapSamplesTimedBeta_BrayCurtis <- ConvertPreparedToBeta(
  bootstrapSamplesPairedBeta_source %>% dplyr::filter(
    DistanceFromCenterExpRev == 0
  ) %>% dplyr::group_by(
    TimeGapNumber, Bootstrap, Type
  ),
  columns = betacolumns,
  method = "Bray", presenceabsence = FALSE,
  postfixes = postfixes,
  indicator = "TimeGapNumber"
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
  caption = paste0("file: ", files_)
) + ggplot2::scale_color_viridis_c(
  option = "C"
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

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

target <- bootstrapSamplesPairedAlpha %>% dplyr::filter(
  DistanceFromCenterExpRev == 0,
  TimeGapNumber == preferredTimeGap
)

plot_1_PairedAlphaStart <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = Type,
    y = `Difference of Average Number of Species in Patch`
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = 0.1,
  notch = TRUE
) + ggplot2::geom_line(
  data = target %>% dplyr::left_join(
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
) + ggplot2::facet_wrap(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ .#DistanceFromCenterExpRev + TimeGapNumber
) + ggplot2::labs(
  title = "Paired Alpha, Maximal Distance, Gap = 10 * lambda",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", files_)
) + ggplot2::scale_color_viridis_c(
  option = "C"
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

target <- bootstrapSamplesPairedAlpha %>% dplyr::filter(
  DistanceFromCenterExpRev == 0,
  TimeGapNumber %in% c(1, 25, 50, 75, 100)
)

plot_1_PairedAlphaTimeGaps <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = TimeGapNumber,
    y = `Difference of Average Number of Species in Patch`,
    fill = Type, group = interaction(Type, TimeGapNumber)
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(10)
) + ggplot2::geom_boxplot(
  position = ggplot2::position_dodge(10),
  width = 2,
  notch = TRUE
) + ggplot2::facet_wrap(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ .#DistanceFromCenterExpRev + TimeGapNumber
) + ggplot2::labs(
  title = "Paired Alpha, Maximal Distance, Varying Gaps",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", files_)
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

target <- bootstrapSamplesPairedAlpha %>% dplyr::filter(
  TimeGapNumber %in% c(preferredTimeGap)
)

plot_1_PairedAlphaPairGaps <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = DistanceFromCenterExpRev,
    y = `Difference of Average Number of Species in Patch`,
    fill = Type, group = interaction(Type, DistanceFromCenterExpRev)
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(1)
) + ggplot2::geom_boxplot(
  position = ggplot2::position_dodge(1),
  width = 0.1,
  notch = TRUE
) + ggplot2::facet_wrap(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ .#DistanceFromCenterExpRev + TimeGapNumber
) + ggplot2::labs(
  title = "Paired Alpha, Varying Distance, Gap = 10 * lambda",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", files_)
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

target <- bootstrapSamplesTimedAlpha %>% dplyr::filter(
  TimeGapNumber %in% c(1, 25, 50, 75, 100)
)

plot_1_TimedAlphaTimeGaps <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = TimeGapNumber,
    y = `Difference of Average Number of Species in Patch`,
    fill = Type, group = interaction(Type, TimeGapNumber)
  )
) + ggplot2::geom_violin(
  position = ggplot2::position_dodge(10)
) + ggplot2::geom_boxplot(
  position = ggplot2::position_dodge(10),
  width = 2,
  notch = TRUE
) + ggplot2::facet_wrap(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ .#DistanceFromCenterExpRev + TimeGapNumber
) + ggplot2::labs(
  title = paste0(
    "Paired Alpha, Average over Patches, Varying Gaps",
    if (logarithmicTimeScale) ", Log Scale X")
  ,
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", files_)
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

target <- bootstrapSamplesPairedAlpha %>% dplyr::filter(
  DistanceFromCenterExpRev == 0
) %>% dplyr::group_by(
  Type, TimeGapNumber, `Species Subset`
) %>% dplyr::rename(
  y = `Difference of Average Number of Species in Patch`
) %>% dplyr::summarise(
  yminp = min(y),
  ymin1 = quantile(y, probs = 0.025),
  ymin2 = quantile(y, probs = 0.25),
  ymedian = quantile(y, probs = 0.5),
  ymean = mean(y),
  ymax1 = quantile(y, probs = 1 - 0.25),
  ymax2 = quantile(y, probs = 1 - 0.025),
  ymaxp = max(y),
  .groups = "drop"
)

plot_1_PairedAlphaTimeGaps_Ribbon <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = TimeGapNumber,
    fill = Type, group = Type, color = Type
  )
) + ggplot2::geom_ribbon( # Inner
  ggplot2::aes(ymin = ymin2,
               ymax = ymax1), alpha = 0.3
) + ggplot2::geom_ribbon( # Outer
  ggplot2::aes(ymin = ymin1,
               ymax = ymax2), alpha = 0.2
) + ggplot2::geom_line(
  ggplot2::aes(y = ymedian), linetype = "dotted"
) + ggplot2::geom_line(
  ggplot2::aes(y = ymean), linetype = "dashed"
) + ggplot2::geom_point( # Minimum detected
  ggplot2::aes(y = yminp), shape = 22, alpha = 0.4
) + ggplot2::geom_point( # Maximum detected
  ggplot2::aes(y = ymaxp), shape = 22, alpha = 0.4
) + ggplot2::facet_grid(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ Type
) + ggplot2::labs(
  title = paste0(
    "Paired Alpha, Maximal Distance, Varying Gaps",
    if (logarithmicTimeScale) ", Log Scale X")
  ,
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", files_),
  ylab = "Difference of Average Number of Species in Patch"
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

target <- bootstrapSamplesPairedAlpha %>% dplyr::filter(
  # DistanceFromCenterExpRev == 0
  TimeGapNumber == preferredTimeGap
) %>% dplyr::group_by(
  Type, DistanceFromCenterExpRev, `Species Subset`
) %>% dplyr::rename(
  y = `Difference of Average Number of Species in Patch`
) %>% dplyr::summarise(
  yminp = min(y),
  ymin1 = quantile(y, probs = 0.025),
  ymin2 = quantile(y, probs = 0.25),
  ymedian = quantile(y, probs = 0.5),
  ymean = mean(y),
  ymax1 = quantile(y, probs = 1 - 0.25),
  ymax2 = quantile(y, probs = 1 - 0.025),
  ymaxp = max(y),
  .groups = "drop"
)

plot_1_PairedAlphaPairGaps_Ribbon <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = DistanceFromCenterExpRev,
    fill = Type, group = Type, color = Type
  )
) + ggplot2::geom_ribbon( # Inner
  ggplot2::aes(ymin = ymin2,
               ymax = ymax1), alpha = 0.3
) + ggplot2::geom_ribbon( # Outer
  ggplot2::aes(ymin = ymin1,
               ymax = ymax2), alpha = 0.2
) + ggplot2::geom_line(
  ggplot2::aes(y = ymedian), linetype = "dotted"
) + ggplot2::geom_line(
  ggplot2::aes(y = ymean), linetype = "dashed"
) + ggplot2::geom_point( # Minimum detected
  ggplot2::aes(y = yminp), shape = 22, alpha = 0.4
) + ggplot2::geom_point( # Maximum detected
  ggplot2::aes(y = ymaxp), shape = 22, alpha = 0.4
) + ggplot2::facet_grid(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ Type
) + ggplot2::labs(
  title = "Paired Alpha, Varying Distance, Gap = 10 * lambda",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", files_),
  ylab = "Difference of Average Number of Species in Patch"
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

target <- bootstrapSamplesTimedAlpha %>% dplyr::filter(
) %>% dplyr::group_by(
  Type, TimeGapNumber, `Species Subset`
) %>% dplyr::rename(
  y = `Difference of Average Number of Species in Patch`
) %>% dplyr::summarise(
  yminp = min(y),
  ymin1 = quantile(y, probs = 0.025),
  ymin2 = quantile(y, probs = 0.25),
  ymedian = quantile(y, probs = 0.5),
  ymean = mean(y),
  ymax1 = quantile(y, probs = 1 - 0.25),
  ymax2 = quantile(y, probs = 1 - 0.025),
  ymaxp = max(y),
  .groups = "drop"
)

plot_1_TimedAlphaTimeGaps_Ribbon <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = TimeGapNumber,
    fill = Type, group = Type, color = Type
  )
) + ggplot2::geom_ribbon( # Inner
  ggplot2::aes(ymin = ymin2,
               ymax = ymax1), alpha = 0.3
) + ggplot2::geom_ribbon( # Outer
  ggplot2::aes(ymin = ymin1,
               ymax = ymax2), alpha = 0.2
) + ggplot2::geom_line(
  ggplot2::aes(y = ymedian), linetype = "dotted"
) + ggplot2::geom_line(
  ggplot2::aes(y = ymean), linetype = "dashed"
) + ggplot2::geom_point( # Minimum detected
  ggplot2::aes(y = yminp), shape = 22, alpha = 0.4
) + ggplot2::geom_point( # Maximum detected
  ggplot2::aes(y = ymaxp), shape = 22, alpha = 0.4
) + ggplot2::facet_grid(
  factor(`Species Subset`, ordered = TRUE, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"
  )) ~ Type
) + ggplot2::labs(
  title = "Paired Alpha, Average over Patches, Varying Gaps",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", files_),
  ylab = "Difference of Average Number of Species in Patch"
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

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
  caption = paste0("file: ", files_)
) + ggplot2::coord_cartesian(
  ylim = ylimRichnessChange
)

### Plot 3: ###################################################################
target <- bootstrapSamplesPairedBeta_Jaccard %>% dplyr::filter(
  !Meaningless
) %>% dplyr::group_by(
  Type, DistanceFromCenterExpRev, `Subset`
) %>% dplyr::rename(
  y = `Jaccard`
) %>% dplyr::summarise(
  yminp = min(y),
  ymin1 = quantile(y, probs = 0.025),
  ymin2 = quantile(y, probs = 0.25),
  ymedian = quantile(y, probs = 0.5),
  ymean = mean(y),
  ymax1 = quantile(y, probs = 1 - 0.25),
  ymax2 = quantile(y, probs = 1 - 0.025),
  ymaxp = max(y),
  .groups = "drop"
)

plot_1_PairedBetaDistGaps_Ribbon <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = DistanceFromCenterExpRev,
    fill = Type, group = Type, color = Type
  )
) + ggplot2::geom_ribbon( # Inner
  ggplot2::aes(ymin = ymin2,
               ymax = ymax1), alpha = 0.3
) + ggplot2::geom_ribbon( # Outer
  ggplot2::aes(ymin = ymin1,
               ymax = ymax2), alpha = 0.2
) + ggplot2::geom_line(
  ggplot2::aes(y = ymedian), linetype = "dotted"
) + ggplot2::geom_line(
  ggplot2::aes(y = ymean), linetype = "dashed"
) + ggplot2::geom_point( # Minimum detected
  ggplot2::aes(y = yminp), shape = 22, alpha = 0.4
) + ggplot2::geom_point( # Maximum detected
  ggplot2::aes(y = ymaxp), shape = 22, alpha = 0.4
) + ggplot2::facet_grid(
  factor(`Subset`, ordered = TRUE, levels = c(
    "All", "Det. In Control", "Not Det. In Control"
  )) ~ Type
) + ggplot2::labs(
  title = "Paired Beta, Varying Distance, Gap = 10 * lambda",
  subtitle = "Presence/Absence Beta Dissimilarity",
  caption = paste0("file: ", files_),
  ylab = "Jaccard Dissimilarity"
) + ggplot2::coord_cartesian(
  ylim = c(0, 1)
)

target <- bootstrapSamplesTimedBeta_Jaccard %>% dplyr::filter(
  !Meaningless
) %>% dplyr::group_by(
  Type, TimeGapNumber, `Subset`
) %>% dplyr::rename(
  y = `Jaccard`
) %>% dplyr::summarise(
  yminp = min(y),
  ymin1 = quantile(y, probs = 0.025),
  ymin2 = quantile(y, probs = 0.25),
  ymedian = quantile(y, probs = 0.5),
  ymean = mean(y),
  ymax1 = quantile(y, probs = 1 - 0.25),
  ymax2 = quantile(y, probs = 1 - 0.025),
  ymaxp = max(y),
  .groups = "drop"
)

plot_1_TimedBetaTimeGaps_Ribbon <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = TimeGapNumber,
    fill = Type, group = Type, color = Type
  )
) + ggplot2::geom_ribbon( # Inner
  ggplot2::aes(ymin = ymin2,
               ymax = ymax1), alpha = 0.3
) + ggplot2::geom_ribbon( # Outer
  ggplot2::aes(ymin = ymin1,
               ymax = ymax2), alpha = 0.2
) + ggplot2::geom_line(
  ggplot2::aes(y = ymedian), linetype = "dotted"
) + ggplot2::geom_line(
  ggplot2::aes(y = ymean), linetype = "dashed"
) + ggplot2::geom_point( # Minimum detected
  ggplot2::aes(y = yminp), shape = 22, alpha = 0.4
) + ggplot2::geom_point( # Maximum detected
  ggplot2::aes(y = ymaxp), shape = 22, alpha = 0.4
) + ggplot2::facet_grid(
  factor(`Subset`, ordered = TRUE, levels = c(
    "All", "Det. In Control", "Not Det. In Control"
  )) ~ Type
) + ggplot2::labs(
  title = paste0(
    "Paired Beta, Maximal Distance, Varying Time Gap",
    if (logarithmicTimeScale) ", Log Scale X")
  ,
  subtitle = "Presence/Absence Beta Dissimilarity",
  caption = paste0("file: ", files_),
  ylab = "Jaccard Dissimilarity"
) + ggplot2::coord_cartesian(
  ylim = c(0, 1)
)

target <- bootstrapSamplesPairedBeta_BrayCurtis %>% dplyr::filter(
  !Meaningless
) %>% dplyr::group_by(
  Type, DistanceFromCenterExpRev, `Subset`
) %>% dplyr::rename(
  y = `Bray`
) %>% dplyr::summarise(
  yminp = min(y),
  ymin1 = quantile(y, probs = 0.025),
  ymin2 = quantile(y, probs = 0.25),
  ymedian = quantile(y, probs = 0.5),
  ymean = mean(y),
  ymax1 = quantile(y, probs = 1 - 0.25),
  ymax2 = quantile(y, probs = 1 - 0.025),
  ymaxp = max(y),
  .groups = "drop"
)

plot_1_PairedBetaDistGaps_Ribbon_BrayCurtis <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = DistanceFromCenterExpRev,
    fill = Type, group = Type, color = Type
  )
) + ggplot2::geom_ribbon( # Inner
  ggplot2::aes(ymin = ymin2,
               ymax = ymax1), alpha = 0.3
) + ggplot2::geom_ribbon( # Outer
  ggplot2::aes(ymin = ymin1,
               ymax = ymax2), alpha = 0.2
) + ggplot2::geom_line(
  ggplot2::aes(y = ymedian), linetype = "dotted"
) + ggplot2::geom_line(
  ggplot2::aes(y = ymean), linetype = "dashed"
) + ggplot2::geom_point( # Minimum detected
  ggplot2::aes(y = yminp), shape = 22, alpha = 0.4
) + ggplot2::geom_point( # Maximum detected
  ggplot2::aes(y = ymaxp), shape = 22, alpha = 0.4
) + ggplot2::facet_grid(
  factor(`Subset`, ordered = TRUE, levels = c(
    "All", "Det. In Control", "Not Det. In Control"
  )) ~ Type
) + ggplot2::labs(
  title = "Paired Beta, Varying Distance, Gap = 10 * lambda",
  subtitle = "Beta Dissimilarity",
  caption = paste0("file: ", files_),
  ylab = "Bray-Curtis Dissimilarity"
) + ggplot2::coord_cartesian(
  ylim = c(0, 1)
)

target <- bootstrapSamplesTimedBeta_BrayCurtis %>% dplyr::filter(
  !Meaningless
) %>% dplyr::group_by(
  Type, TimeGapNumber, `Subset`
) %>% dplyr::rename(
  y = `Bray`
) %>% dplyr::summarise(
  yminp = min(y),
  ymin1 = quantile(y, probs = 0.025),
  ymin2 = quantile(y, probs = 0.25),
  ymedian = quantile(y, probs = 0.5),
  ymean = mean(y),
  ymax1 = quantile(y, probs = 1 - 0.25),
  ymax2 = quantile(y, probs = 1 - 0.025),
  ymaxp = max(y),
  .groups = "drop"
)

plot_1_TimedBetaTimeGaps_Ribbon_BrayCurtis <- ggplot2::ggplot(
  target,
  ggplot2::aes(
    x = TimeGapNumber,
    fill = Type, group = Type, color = Type
  )
) + ggplot2::geom_ribbon( # Inner
  ggplot2::aes(ymin = ymin2,
               ymax = ymax1), alpha = 0.3
) + ggplot2::geom_ribbon( # Outer
  ggplot2::aes(ymin = ymin1,
               ymax = ymax2), alpha = 0.2
) + ggplot2::geom_line(
  ggplot2::aes(y = ymedian), linetype = "dotted"
) + ggplot2::geom_line(
  ggplot2::aes(y = ymean), linetype = "dashed"
) + ggplot2::geom_point( # Minimum detected
  ggplot2::aes(y = yminp), shape = 22, alpha = 0.4
) + ggplot2::geom_point( # Maximum detected
  ggplot2::aes(y = ymaxp), shape = 22, alpha = 0.4
) + ggplot2::facet_grid(
  factor(`Subset`, ordered = TRUE, levels = c(
    "All", "Det. In Control", "Not Det. In Control"
  )) ~ Type
) + ggplot2::labs(
  title = paste0(
    "Paired Beta, Maximal Distance, Varying Time Gap",
    if (logarithmicTimeScale) ", Log Scale X")
  ,
  subtitle = "Beta Dissimilarity",
  caption = paste0("file: ", files_),
  ylab = "Bray-Curtis Dissimilarity"
) + ggplot2::coord_cartesian(
  ylim = c(0, 1)
)

