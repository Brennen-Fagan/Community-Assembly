# Setup: ######################################################################


source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

### 7 Supplement: #############################################################
##### 5s Idea 2: ##############################################################
# Another interesting option: we plot vertical kde's for each intervention at
# "-1", 10, 100, 200, and 400 time units from the intervention across all sims
# with that land-use/intervention and species affinity setup. Then, *inside the
# kde's* we plot all of points with little-no sizes along the axis (or maybe
# along the edge of the kde?). More imporantly, we plot the *edges* with a fixed
# alpha and size so that we can see the spread of interactions across the kdes.

####### Mockup 1: #############################################################
newplot5_as_Specification <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  # SpeciesPreferences == "100% 0",
  # SpeciesPreferences == "50% 0, 50% 1",
  SpeciesPreferences == "Uniform(0, 1)",
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)"),
  is.na(Subset)
) |> tidytable::group_by(
  PoolPatch:InterventionFinal
) |> tidytable::filter(
  Time %% 1 != 0 | # I.e., time recorded just before the intervention!
    dplyr::lag(Time) %% 1 != 0 | # First time after the intervention (v9 !!)
    dplyr::lag(Time, n = 11) %% 1 != 0 | # 10 * 10 + 1 time steps
    dplyr::lag(Time, n = 21) %% 1 != 0 |
    dplyr::lag(Time, n = 41) %% 1 != 0
) |> tidytable::mutate(
  Time2 = tidytable::case_when(
    Time %% 1 != 0 ~ -1,
    dplyr::lag(Time) %% 1 != 0 ~ 10,
    dplyr::lag(Time, n = 2) %% 1 != 0 ~ 100,
    dplyr::lag(Time, n = 3) %% 1 != 0 ~ 200,
    dplyr::lag(Time, n = 4) %% 1 != 0 ~ 400
  )
) |> tidytable::ungroup(
)

newplot5_as_Networks <- generateNetworks(newplot5_as_Specification)

newplot5_kdes <- lapply(
  newplot5_as_Networks$Envs, function(e) {
    e$trophics$EdgeVertexLists[[1]][[1]]$Vertices |> tidytable::select(
      node, Type, Size, N
    ) |> cbind(
      e$Row |> tidytable::select(Time, Time2, PoolPatch:InterventionFinal)
    ) |> tidytable::mutate(
      AffinityVals = e$result$Ellipsis$Affinity$SpeciesAffinities[
        as.numeric(substring(node, 2))
        ]
    )
  }
) |> tidytable::bind_rows(
  # ) |> tidytable::group_by(
  #   PoolPatch:InterventionFinal
  # ) |> tidytable::arrange(
  #   Time
  # ) |> tidytable::mutate( # fix for having not done it ahead of time...
  #   Time2 = tidytable::case_when(
  #     Time == unique(Time)[1] ~ -1,
  #     Time == unique(Time)[2] ~ 10,
  #     Time == unique(Time)[3] ~ 100,
  #     Time == unique(Time)[4] ~ 200,
  #     Time == unique(Time)[5] ~ 400
  #   )
  # ) |> tidytable::ungroup(
)

# newplot5_graph <- lapply(newplot5_as_Networks$Envs, function(e)
#   e$graphs[[1]] %N>% select(
#     node, Type, Size
#   ) %N>% mutate(
#     e$Row |> select(Time, PoolPatch:InterventionFinal)
#   ) %E>% mutate(
#     e$Row |> select(Time, PoolPatch:InterventionFinal)
#   )
# ) |> tidygraph::bind_graphs(
# ) %N>% tidygraph::group_by( # Only supports manual specification
#   PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
#   EventsSeed, InitialConditions, InitialConditionsSeed, Dispersal,
#   NicheDistance, Affinity, AffinitySeed, InterventionPatchType,
#   InterventionPatchSeed, InterventionTimeType, InterventionTimeSeed,
#   InterventionDispersal, InterventionNicheDistance, Intervention,
#   SpeciesPreferences, InterventionInitial, InterventionFinal
# ) %N>% tidygraph::arrange(
#   Time
# ) %N>% tidygraph::mutate( # fix for having not done it ahead of time...
#   Time2 = dplyr::case_when(
#     Time == unique(Time)[1] ~ -1,
#     Time == unique(Time)[2] ~ 10,
#     Time == unique(Time)[3] ~ 100,
#     Time == unique(Time)[4] ~ 200,
#     Time == unique(Time)[5] ~ 400
#   )
# ) %N>% tidygraph::ungroup(
# ) %E>% tidygraph::group_by( # Only supports manual specification
#   PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
#   EventsSeed, InitialConditions, InitialConditionsSeed, Dispersal,
#   NicheDistance, Affinity, AffinitySeed, InterventionPatchType,
#   InterventionPatchSeed, InterventionTimeType, InterventionTimeSeed,
#   InterventionDispersal, InterventionNicheDistance, Intervention,
#   SpeciesPreferences, InterventionInitial, InterventionFinal
# ) %E>% tidygraph::arrange(
#   Time
# ) %E>% tidygraph::mutate( # fix for having not done it ahead of time...
#   Time2 = dplyr::case_when(
#     Time == unique(Time)[1] ~ -1,
#     Time == unique(Time)[2] ~ 10,
#     Time == unique(Time)[3] ~ 100,
#     Time == unique(Time)[4] ~ 200,
#     Time == unique(Time)[5] ~ 400
#   )
# ) %E>% tidygraph::ungroup(
# )

# ggraph::ggraph(
#   ggraph::create_layout(newplot5_graph, "manual",
#                         y = newplot5_graph %N>% tidygraph::pull(Size), x = -0.5)
ggplot2::ggplot(
) + ggplot2::geom_density(
  data = newplot5_kdes,
  mapping = ggplot2::aes(
    y = Size, fill = Type, color = Intervention
  ),
  trim = TRUE
  # ) + ggplot2::geom_rug(
  #   ggplot2::aes(color = Type)
  # ) + ggraph::geom_edge_arc(
  #   mapping = aes(
  #     color = Type,
  #     #color = node1.Type, # but then exploit+ between consumers is orange.
  #     linetype = Type,
  #     alpha = log10(effectNormalised),
  #     end_cap = circle(2, 'pt')
  #   ),
  #   arrow = arrow(length = unit(2, 'mm')),
  #   alpha = 0.2,
  #   show.legend = FALSE
  # ) + ggraph::geom_node_point(
  #   mapping = aes(
  #     color = Type
  #   ),
  #   show.legend = FALSE
  #   # ) + ggplot2::geom_hline(
  #   #   yintercept = -1, linetype = "dashed", color = "black"
) + ggplot2::facet_grid(
  Intervention + SpeciesPreferences ~ Time2
) + ggplot2::scale_y_log10(
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
) + ggplot2::theme_minimal(
)

ggplot2::ggsave(
  ggplot2::ggplot(
  ) + ggplot2::geom_density(
    data = newplot5_kdes,
    mapping = ggplot2::aes(
      y = Size, fill = Type, color = Intervention
    ),
    trim = TRUE
  ) + ggplot2::facet_grid(
    Intervention + SpeciesPreferences ~ Time2
  ) + ggplot2::scale_y_log10(
  ) + ggplot2::scale_color_manual(
    values = colorPalette,
    name = "Habitat Land-use"
  ) + ggplot2::scale_fill_manual(
    values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
  ) + ggplot2::theme_minimal(
  ),
  # filename = "Figure5s1_Prototype1.png", # 100% 0
  # filename = "Figure5s2_Prototype1.png", # 50% 0, 50% 1
  filename = "Figure5s3_Prototype1.png", # Uniform(0, 1)
  units = "cm", width = 6.5*3, height = 6.5*2
)

ggplot2::ggsave(
  ggplot2::ggplot(
  ) + ggplot2::geom_density_2d(
    # ) + ggplot2::geom_bin_2d(
    data = newplot5_kdes,
    mapping = ggplot2::aes(
      x = AffinityVals, y = Size,
      # fill = Type,
      color = Intervention,
      group = interaction(Type, Intervention)
    ),
    # bins = 10
    alpha = 0.4,
    contour_var = "count",
    adjust = 0.7
    # trim = TRUE
  ) + ggplot2::geom_hline(
    yintercept = 0.1, color = "red", show.legend = FALSE
  ) + ggplot2::facet_grid(
    Intervention + SpeciesPreferences ~ Time2
  ) + ggplot2::scale_y_log10(
  ) + ggplot2::scale_color_manual(
    values = colorPalette,
    name = "Habitat Land-use"
  ) + ggplot2::scale_fill_manual(
    values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
    # ) + ggplot2::scale_fill_viridis_c(
  ) + ggplot2::theme_minimal(
  ) + ggplot2::xlab(
    "Land-use Type"
  ),
  filename = "Figure5s4_Prototype1.png", # Uniform(0, 1)
  units = "cm", width = 10*3, height = 10*2
)
