
newplot2_dataB <- diversitiesAll %>% newplot2_filtration(
) %>% tidytable::left_join(
  endTimes %>% dplyr::select(-Times)
) %>% tidytable::filter(
  Time > Start, Time < Stop,
  Metric == "Alpha Abundance",
  !is.na(Subset)
) %>% tidytable::separate_wider_delim(
  cols = Subset, names = c("Guild", "AffinityBin"), delim = "_"
) %>% tidytable::group_by(
  Environment1, Environment2, Metric, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal, Guild, Time
) %>% tidytable::summarise(# Across Affinity Bins
  Value = sum(Value, na.rm = TRUE), .groups = "drop_last"
) %>% tidytable::summarise(# Across Time
  Mean = mean(Value),
  q025 = quantile(Value, probs = 0.25),
  q075 = quantile(Value, probs = 0.75)
) %>% tidytable::pivot_wider(
  names_from = Guild, values_from = c(Mean, q025, q075)
)

# Change in amount of abundance between basals.
# Note we are pairing time-averages of simulations, then dividing, then
# averaging, but this ignores the internal (averaged-out) time dynamics, which
# judging by the inner 50% of values over time, appears to be quite variable.
newplot2_dataB %>% tidytable::select(
  -Affinity, -AffinitySeed, -InterventionInitial, -InterventionFinal
) %>% tidytable::ungroup()  %>% tidytable::pivot_wider(
  values_from = c(q025_Basal, Mean_Basal, q075_Basal,
                  q025_Consumer, Mean_Consumer, q075_Consumer),
  names_from = Intervention
) %>% tidytable::distinct(
) %>% tidytable::mutate(
  RatioB_1_05 = `Mean_Basal_(1)` / `Mean_Basal_(0.5)`,
  RatioB_0_1 = `Mean_Basal_(0)` / `Mean_Basal_(1)`,
  RatioC_05_1 = `Mean_Consumer_(0.5)` / `Mean_Consumer_(1)`,
  RatioC_0_1 = `Mean_Consumer_(0)` / `Mean_Consumer_(1)`
) %>% tidytable::summarise(
  MinB_1_05 = min(RatioB_1_05),
  RatioB_1_05 = mean(RatioB_1_05),
  MaxB_1_05 = max(RatioB_1_05),
  MinB_0_1 = min(RatioB_0_1),
  RatioB_0_1 = mean(RatioB_0_1),
  MaxB_0_1 = max(RatioB_0_1),
  MinC_05_1 = min(RatioC_05_1),
  RatioC_05_1 = mean(RatioC_05_1),
  MaxC_05_1 = max(RatioC_05_1),
  MinC_0_1 = min(RatioC_0_1),
  RatioC_0_1 = mean(RatioC_0_1),
  MaxC_0_1 = max(RatioC_0_1)
)


newplot2_b <- ggplot2::ggplot(
  newplot2_dataB,
  ggplot2::aes(
    x = Mean_Basal,
    y = Mean_Consumer,
    color = Intervention)
) + ggplot2::geom_point(
  show.legend = FALSE
  # ) + ggplot2::geom_errorbar(
  #   inherit.aes = FALSE,
  #   mapping = ggplot2::aes(
  #     x = Mean_Basal, ymin = q025_Consumer, ymax = q075_Consumer,
  #     color = Intervention
  #   )
  # ) + ggplot2::geom_errorbarh(
  #   inherit.aes = FALSE,
  #   mapping = ggplot2::aes(
  #     y = Mean_Consumer, xmin = q025_Basal, xmax = q075_Basal,
  #     color = Intervention
  #   )
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::coord_fixed(
) + ggplot2::labs(
  x = "Basal\nAbundance",
  y = "Consumer Abundance"
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
)

# Trying to support a breakdown of richness and abundance to show what is
# happening to, and because of, the species that have the wrong land-use
# preference for the given habitat.
newplot3_dataB <- diversitiesAll %>% newplot3_filtration(
) %>% tidytable::left_join(
  endTimes %>% dplyr::select(-Times)
) %>% tidytable::filter(
  Time > Start, Time < Stop,
  Metric == "Alpha Biomass",
  !is.na(Subset)
) %>% tidytable::separate_wider_delim(
  cols = Subset, names = c("Guild", "AffinityBin"), delim = "_"
) %>% tidytable::group_by(
  Environment1, Environment2, Metric, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal, Guild, AffinityBin
) %>% tidytable::summarise(# Across Time
  Mean = mean(Value)
) %>% tidytable::pivot_wider(
  names_from = "Guild", values_from = "Mean"
)

newplot3_b <- ggplot2::ggplot(
  newplot3_dataB ,
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = Intervention,
    shape = AffinityBin
  )
) + ggplot2::geom_point(
  show.legend = FALSE
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::coord_fixed(
) + ggplot2::labs(
  x = "Basal Abundance",
  y = "Consumer Abundance"
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme(
  plot.tag.position = c(0.05, 1)
)


newplot4_dataB <- diversitiesAll %>% newplot4_filtration(
) %>% tidytable::left_join(
  endTimes %>% dplyr::select(-Times)
) %>% tidytable::filter(
  Time > Start, Time < Stop,
  Metric == "Alpha Biomass",
  !is.na(Subset)
) %>% tidytable::separate_wider_delim(
  cols = Subset, names = c("Guild", "AffinityBin"), delim = "_"
) %>% tidytable::group_by(
  Environment1, Environment2, Metric, PoolPatch, PoolPatchSeed,
  Interactions, InteractionsSeed, Events, EventsSeed,
  InitialConditions, InitialConditionsSeed, Dispersal, NicheDistance,
  Affinity, AffinitySeed, InterventionPatchType, InterventionPatchSeed,
  InterventionTimeType, InterventionTimeSeed, InterventionDispersal,
  InterventionNicheDistance, Intervention, SpeciesAffinity,
  InterventionInitial, InterventionFinal, Guild, AffinityBin
) %>% tidytable::summarise(# Across Time
  Mean = mean(Value)
) %>% tidytable::pivot_wider(
  names_from = "Guild", values_from = "Mean"
) %>% tidytable::separate(
  col = "AffinityBin", into = c("Left", "Right"), sep = ","
) %>% tidytable::mutate(
  Left = round(as.numeric(gsub(pattern = "^.", replacement = "", x = Left))*5)/5,
  Right = round(as.numeric(gsub(pattern = ".$", replacement = "", x = Right))*5)/5,
  AffinityBin = paste0("(", Left, ", ", Right, "]")
)

newplot4_b <- ggplot2::ggplot(
  newplot4_dataB ,
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = Intervention,
    shape = AffinityBin
  )
) + ggplot2::geom_point(
  show.legend = TRUE
) + ggplot2::scale_x_log10(
) + ggplot2::scale_y_log10(
) + ggplot2::coord_fixed(
) + ggplot2::labs(
  x = "Basal Abundance",
  y = "Consumer Abundance"
) + ggplot2::theme_minimal(
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Land-use"
) + ggplot2::labs(
  tag = "b)"
) + ggplot2::theme(
  plot.tag.position = c(0.02, 1)
) + ggplot2::guides(
  color = "none", fill = "none"
) + ggplot2::labs(
  shape = "Land-use\nPreference"
)