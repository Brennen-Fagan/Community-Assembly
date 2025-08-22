##### In/Out Statistics: ######################################################
# Problem here: need to translate to the newer formatting, and then make the
# plots nicer for human consumption, since the bar charts weren't doing a great
# job iirc.



# plotValueChart <- function(
#   data, facets = as.formula(Intervention ~ SpeciesPreferences)
# ) {
#   ggplot2::ggplot(
#     data,
#     ggplot2::aes(x = interaction(InType,
#                                  OutType,
#                                  sep = "\n"),#"",
#                  y = ChartValue,
#                  fill = interaction(InType,
#                                     OutType,
#                                     sep = "/"))
#   ) + ggplot2::geom_col(
#     show.legend = FALSE
#   ) + ggplot2::facet_grid(
#     facets
#   ) + ggplot2::theme_minimal(
#   ) + ggplot2::labs(
#     x = "Enter/Exit Type", y = "Count"#, fill = "In/Out"
#   )
# }

# .divs, # diversitiesAll
# .ps, # Pers
# .ces, # ColExt
# .ets, # endTimes
# .gs, # Singleton Graphs
# ... # commonfilters
# plotValueChart(
#   rbind(
#     .ps %>% tidytable::filter(
#       ...,
#       Persistence > 0,
#       InType != externalNames["Dispersal"],
#       In < Stop, Out > Start
#     ) %>% tidytable::group_by(
#       SpeciesPreferences, InType, OutType, Intervention
#     ) %>% tidytable::summarise(
#       ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
#     ) %>% dplyr::mutate( # Tidytable renders as character again!
#       Intervention =
#         factor(Intervention,
#                levels = rev(c("(0)", "(0)->(0.5)", "(0)->(1)",
#                               "(0.5)",
#                               "(1)")),
#                ordered = TRUE)
#     ),
#     .ces %>% tidytable::filter(
#       ...,
#       !Success | EventType == "Present",
#       Times > Start, Times < Stop
#     ) %>% tidytable::mutate(
#       InType = externalNames[
#         ifelse(EventType == "Arrival", "Failed Arrival", "Present")
#         ],
#       OutType = externalNames["NA"]
#     ) %>% tidytable::group_by(
#       SpeciesPreferences, InType, OutType, Intervention
#     ) %>% tidytable::summarise(
#       ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
#     ) %>% dplyr::mutate( # Tidytable renders as character again!
#       Intervention =
#         factor(Intervention,
#                levels = rev(c("(0)", "(0)->(0.5)", "(0)->(1)",
#                               "(0.5)",
#                               "(1)")),
#                ordered = TRUE)
#     )
#   ),
#   facets = as.formula(Intervention ~ .)
# ) + ggplot2::theme(
#   legend.position = "bottom"
# ) + ggplot2::guides(
#   fill = ggplot2::guide_legend(nrow = 2)
# ) + ggplot2::ylab(
#   "Average Count in a Simulation"
# )

rbind(
  Pers |> tidytable::filter(
    ...,
    Persistence > 0,
    InType != externalNames["Dispersal"],
    In < Stop, Out > Start
  ) |> tidytable::group_by(
    SpeciesPreferences, InType, OutType, Intervention
  ) |> tidytable::summarise(
    ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
  ),
  ColExt |> tidytable::left_join(
    # Start and Stop aren't already present in this version
    endTimes
  ) |> tidytable::mutate(
    SpeciesPreferences =
      affinityDictionaryOrigin$SpeciesAffinities[as.numeric(Affinity)],
    Intervention = unlist(mapply(
      FUN = interventionNamingScheme,
      Affinity, PoolPatch, InterventionPatchType
    ))
  ) |> changeAffinityLevels(
  ) |> changeInterventionLevels(
  ) |> tidytable::filter(
    ...,
    !Success | EventType == "Present",
    Times > Start, Times < Stop
  ) |> tidytable::mutate(
    InType = externalNames[
      ifelse(EventType == "Arrival", "Failed Arrival", "Present")
      ],
    OutType = externalNames["NA"]
  ) |> tidytable::group_by(
    SpeciesPreferences, InType, OutType, Intervention
  ) |> tidytable::summarise(
    ChartValue = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
  )
)
