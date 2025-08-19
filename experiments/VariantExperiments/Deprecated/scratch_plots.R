ggplot2::ggplot(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("rep_0"),
    Intervention %in% c("(0.75)",  "(0)->(0.75)",  "(0.25)->(0.75)",
                        "(0.5)->(0.75)", "(0.75)->(0.75)", "(1)->(0.75)"),
    NicheDistance == 7, #PoolPatchSeed == "341",
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  ggplot2::aes(
    x = Time, y = Value, color = PoolPatchSeed
  )
) + ggplot2::geom_line(
) + ggplot2::facet_grid(
  Metric~Intervention#,
  # nrow = 5
) + ggplot2::theme(
  legend.position = "none"
  # ) + ggplot2::geom_boxplot(
  #   ggplot2::aes(x = 34000),
  #   width = 1000
) + ggplot2::labs(
  # title = "Pool: 0; Time > 30,000; Effect Magnitude: 10"
)
#

diversities[[1]]$Presence %>% ggplot2::ggplot(
  ggplot2::aes(
    x = Time,
    y = log10(Size),
    color = log10(Abundance)
  )
) + ggplot2::geom_hline(
  color = "red",
  yintercept = -1
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1
) # but rbind the presences, filter to the appropriate subset, then facet_wrap.




load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_341-341_2025-01-26/TSTS_ColExt_142486-4929-28-1-NA-3-17_341-341-384-387-400_115-1-p-p_1-1.RData")
ggplot2::ggplot(ColExt$Events %>% dplyr::mutate(GainLoss = (Type.x %in% c("Arrival", "Dispersal"))), ggplot2::aes(x = interaction(GainLoss, AffinityBins), color = AffinityBins, fill = interaction(Type.x, Type.y))) + ggplot2::geom_bar()


# Final bar type plot
ggplot2::ggplot(
  ColExt$Events %>% dplyr::filter(
    Type.x != "Present",
    # Type.y == "Consumer"
  ) %>% dplyr::mutate(
    Sign = ifelse(Type.x == "Arrival", 1, -1)
  ) %>% dplyr::group_by(
    # Affinity,
    AffinityBins, Type.x, Type.y
  ) %>% dplyr::summarise(
    height = sum(Sign),
    .groups = "drop"
  ), ggplot2::aes(
    x = AffinityBins,
    fill = Type.x,
    group = interaction(Type.y, Type.x),
    color = Type.y,
    y = height
  )
) + ggplot2::geom_bar(
  stat = "identity", position = "dodge"
) + ggplot2::geom_hline(
  yintercept = 0
) + ggplot2::geom_bar(
  data = ColExt$Events %>% dplyr::filter(
    Type.x != "Present",
    # Type.y == "Consumer"
  ) %>% dplyr::mutate(
    Sign = ifelse(Type.x == "Arrival", 1, -1)
  ) %>% dplyr::group_by(
    # Affinity,
    AffinityBins, Type.y
  ) %>% dplyr::summarise(
    height = sum(Sign),
    .groups = "drop"
  ),
  alpha = 0.7,
  fill = "black",
  inherit.aes = FALSE,
  mapping = ggplot2::aes(
    x = AffinityBins,
    group = Type.y,
    y = height
  ),
  stat = "identity"
)


# Ines suggested just lines through time
ggplot2::ggplot(
  ColExt$Events %>% dplyr::filter(
    Type.x != "Present"
  ) %>% dplyr::group_by(
    Times, Type.x, Type.y
  ) %>% dplyr::summarise(
    Count = dplyr::n(),
    .groups = "drop"
  ) %>% dplyr::group_by(
    Type.x, Type.y
  ) %>% dplyr::arrange(
    Times
  ) %>% dplyr::mutate(
    CountTotal = cumsum(Count),
    RateTotal = CountTotal/Times,
    CountRecent = slider::slide_index_sum(Count, Times, before = 1000),
    RateRecent = CountRecent/1000
  ) %>% dplyr::select(
    -Count
  ) %>% tidyr::pivot_longer(
    CountTotal:RateRecent, names_to = "Measurement", values_to = "Value"
  ),
  ggplot2::aes(
    x = Times, y = Value,
    color = Type.x,
    linetype = Type.y
  )
) + ggplot2::geom_line(
) + ggplot2::facet_wrap(
  .~Measurement, scales = "free_y"
)

# I think the most useful is perhaps CountTotal (showing the "race")
# and then the local average rate RateRecent (showing the "speed").
# Either way, I should also account for the intervention, meaning I need to
# staple the appropriate ColExt's together.
load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_348-348_2025-01-21/TSTS_ColExt_142486-4929-28-1-NA-3-13_348-348-391-394-510.RData")
ColExtOld <- ColExt
load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_348-348_2025-01-21/TSTS_ColExt_142486-4929-28-1-NA-3-13_348-348-391-394-510_115-1-p-p_1-1.RData")
ColExtNew <- ColExt
interventionTime <- min(ColExtNew$Events$Times) #Approximation!
ColExtNew$Ellipsis$GrandparentRun == ColExtOld$Ellipsis$ParentRun
ColExtNew$Events <- rbind(ColExtOld$Events %>% dplyr::filter(Times < min(ColExtNew$Events$Times)), ColExtNew$Events)


tempold <- ggplot2::ggplot(
  ColExtOld$Events %>% dplyr::filter(
    Type.x != "Present"
  ) %>% dplyr::group_by(
    Times, Type.x, Type.y
  ) %>% dplyr::summarise(
    Count = dplyr::n(),
    .groups = "drop"
  ) %>% dplyr::group_by(
    Type.x, Type.y
  ) %>% dplyr::arrange(
    Times
  ) %>% dplyr::mutate(
    CountTotal = cumsum(Count),
    RateTotal = CountTotal/Times,
    CountRecent = slider::slide_index_sum(Count, Times, before = 1000),
    RateRecent = CountRecent/1000
  ) %>% dplyr::select(
    -Count
  ) %>% tidyr::pivot_longer(
    CountTotal:RateRecent, names_to = "Measurement", values_to = "Value"
  ),
  ggplot2::aes(
    x = Times, y = Value,
    color = Type.x,
    linetype = Type.y
  )
) + ggplot2::geom_line(
) + ggplot2::facet_wrap(
  .~Measurement, scales = "free_y"
)

tempnew <- ggplot2::ggplot(
  ColExtNew$Events %>% dplyr::filter(
    Type.x != "Present"
  ) %>% dplyr::group_by(
    Times, Type.x, Type.y
  ) %>% dplyr::summarise(
    Count = dplyr::n(),
    .groups = "drop"
  ) %>% dplyr::group_by(
    Type.x, Type.y
  ) %>% dplyr::arrange(
    Times
  ) %>% dplyr::mutate(
    CountTotal = cumsum(Count),
    RateTotal = CountTotal/Times,
    CountRecent = slider::slide_index_sum(Count, Times, before = 1000),
    RateRecent = CountRecent/1000
  ) %>% dplyr::select(
    -Count
  ) %>% tidyr::pivot_longer(
    CountTotal:RateRecent, names_to = "Measurement", values_to = "Value"
  ),
  ggplot2::aes(
    x = Times, y = Value,
    color = Type.x,
    linetype = Type.y
  )
) + ggplot2::geom_line(
) + ggplot2::facet_wrap(
  .~Measurement, scales = "free_y"
) + ggplot2::geom_vline(xintercept = interventionTime)

ggpubr::ggarrange(tempold,tempnew, nrow = 1, common.legend = TRUE, legend = "right")
