# newplot2_a_time <- 25000
newplot7_Spec <- diversitiesRichness |> tidytable::select(c(
  # Which network:
  "Time", "Environment1",
  # Which File (Base):
  "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
  "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
  "Dispersal", "NicheDistance", "Affinity", "AffinitySeed",
  # Which File (Intervention):
  "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
  "InterventionTimeSeed", "InterventionDispersal", "InterventionNicheDistance",
  # Ease of Use
  "SpeciesAffinity", "Intervention"
)) |> tidytable::filter(
  (
    SpeciesAffinity == "100% 0" &
      NicheDistance == defaultNicheDistance &
      Intervention %in% c("(0)", "(0.5)", "(1)") &
      PoolPatchSeed %in% as.character(343:386) &
      Time == newplot2_a_time
    # ) | (
    
  )
) |> tidytable::distinct(
)

newplot7_Nets <- generateNetworks(newplot7_Spec)

# Scratch work for what we are trying to do upon getting to v10.
exampleNetworks$Index |> split(
  1:nrow(exampleNetworks$Index)
) |> tidytable::map_dfr(
  .f = function(spec) {
    env <- exampleNetworks$Envs[[spec$ID]]
    edges <- env$trophics$EdgeVertexLists[[1]][[1]]$Edges |> tidytable::select(
      from, to, effectActual
    ) |> tidytable::mutate(
      Weight = abs(effectActual) # Want the overall effect, but not clear if
      #                            normalisation is a concern here.
    ) |> tidytable::left_join(
      env$trophics$EdgeVertexLists[[1]][[1]]$Vertices |> tidytable::select(
        node, Size
      ),
      by = c("from" = "node")
    ) |> tidytable::rename(
      fromSize = Size
    ) |> tidytable::left_join(
      env$trophics$EdgeVertexLists[[1]][[1]]$Vertices |> tidytable::select(
        node, Size
      ),
      by = c("to" = "node")
    ) |> tidytable::rename(
      toSize = Size
    ) |> tidytable::bind_cols(
      spec
    )
  }
) |> tidytable::rowwise(
) |> tidytable::mutate(
  InSize = min(fromSize, toSize),
  OutSize = max(fromSize, toSize)
) |> tidytable::pivot_longer(
  cols = c(InSize, OutSize),
  names_to = "inout",
  values_to = "Size"
) |> tidytable::mutate(
  Weight = ifelse(inout == "InSize", +Weight, -Weight)
) |> tidytable::group_by(
  PoolPatchSeed, SpeciesAffinity, NicheDistance, Intervention,
  Size
) |> tidytable::summarise(
  Total = sum(Weight),
  .groups = "drop"
) |> tidytable::arrange(
  Size
) |> tidytable::group_by(
  PoolPatchSeed, SpeciesAffinity, NicheDistance, Intervention
) |> tidytable::mutate(
  Total = cumsum(Total)
) |> ggplot2::ggplot(
  ggplot2::aes(x = Size, y = Total, color = Intervention)
) + ggplot2::geom_line(
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_x_log10(
) + ggplot2::geom_vline(
  xintercept = 0.1
)