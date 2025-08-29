# Setup: ######################################################################
# Whereas Supplement 9 asks what is going on with the species, we now want to
# ask what we can do to generalise the interactions between the species
# Effectively, we are taking a flux-like approach: how much total "effect" is
# being had in different parts of the foodweb if we reduce it to a line along
# the sizes? So every new edge is a step up, every loss of an edge is a step
# down. Because of size-structure, we shouldn't lose any edges until there are
# consumers, so the maximum location is trivial. Rather, we are focusing on the
# general shape.

# Too many graphs to do all at once; need to do multiple runs.
targetPrefIndex <- 2
targetPref <-
  c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)")[targetPrefIndex]

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source(file.path("R", "flattenDiversity.R")) # Req'd by below
source(file.path("R", "generateNetworks.R")) # To create inset graphs.

supplement10 <- list()
supplement10$graph$time <- 25000

### 10 Supplement: ############################################################
supplement10$graph$specification <- diversitiesRichness |> tidytable::select(c(
  # Which network:
  "Time", "Environment1",
  # Which File (Base):
  "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
  "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
  "Dispersal", "NicheDistance", "SpeciesAffinity", "SpeciesAffinitySeed",
  "PatchAffinity", "PatchAffinitySeed",
  # Which File (Intervention):
  "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
  "InterventionTimeSeed", "InterventionDispersal", "InterventionNicheDistance",
  # Ease of Use
  "SpeciesPreferences", "Intervention"
)) |> tidytable::filter(
  SpeciesPreferences == targetPref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.5)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Time == supplement10$graph$time
) |> tidytable::distinct(
)

supplement10$graph$networks <-
  generateNetworks(supplement10$graph$specification,
                   Date = "2025-07-30")

# Scratch work for what we are trying to do upon getting to v10.
supplement10$plot <- supplement10$graph$networks$Index |> split(
  1:nrow(supplement10$graph$networks$Index)
) |> tidytable::map_dfr(
  .f = function(spec) {
    env <- supplement10$graph$networks$Envs[[spec$ID]]
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
) |> tidytable::filter(
  from != to
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
  PoolPatchSeed,
  SpeciesPreferences, NicheDistance, Intervention,
  Size
) |> tidytable::summarise(
  Total = sum(Weight),
  .groups = "drop"
) |> tidytable::arrange(
  Size
) |> tidytable::group_by(
  # PoolPatchSeed,
  SpeciesPreferences, NicheDistance, Intervention
) |> tidytable::mutate(
  Total = cumsum(Total)/length(unique(PoolPatchSeed)) # Average 'Flux'
) |> ggplot2::ggplot(
  ggplot2::aes(x = Size, y = Total, color = Intervention,
               group = Intervention
               # group = interaction(Intervention, PoolPatchSeed)
  )
) + ggplot2::geom_step(
  show.legend = FALSE
) + ggplot2::facet_wrap(
  .~SpeciesPreferences
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_x_log10(
) + ggplot2::geom_vline(
  xintercept = 0.1
) + ggplot2::ylab(
  "Averaged Effect Flux"
) + ggplot2::theme_minimal(
)

ggplot2::ggsave(
  plot = supplement10$plot,
  filename = paste0("Figure_supplement10_v1_", targetPrefIndex, ".png"),
  units = "cm", width = 6.5*3, height = 6.5*2
)

