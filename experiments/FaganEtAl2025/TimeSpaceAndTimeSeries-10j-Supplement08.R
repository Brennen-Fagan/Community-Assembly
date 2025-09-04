# Setup: ######################################################################
# Plot the generalisation of Figure 5 by comparing and contrasting size KDEs.
# For the Uniform case, 2-D is needed to properly convey the shift, so we also
# plot the contours.

# Too many graphs to do all at once; need to do multiple runs.
targetPrefIndex <- 1
targetPref <-
  c("100% 0", "50% 0, 50% 1", "Uniform(0, 1)")[targetPrefIndex]

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source(file.path("R", "flattenDiversity.R")) # Req'd by below
source(file.path("R", "generateNetworks.R")) # To create inset graphs.
load("TSTS_Interventions_10a1.RData")

supplement8 <- list()

### 8 Supplement: #############################################################
####### Mockup 1: #############################################################
supplement8$graph$specification <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  SpeciesPreferences == targetPref,
  Intervention %in% c("(0)->(0.5)", "(0.5)->(0)"),
  is.na(Subset)
) |> tidytable::left_join(
  InterventionTimes |> tidytable::select(
    TimeIntervention, PoolPatch:PatchAffinitySeed
  ),
  by = c("PoolPatch", "PoolPatchSeed", "Interactions",
         "InteractionsSeed", "Events",
         "EventsSeed", "InitialConditions",
         "InitialConditionsSeed", "Dispersal",
         "NicheDistance", "SpeciesAffinity",
         "SpeciesAffinitySeed", "PatchAffinity",
         "PatchAffinitySeed")
) |> tidytable::mutate(
  TimeSinceIntervention = Time - TimeIntervention
) |> tidytable::group_by(
  PoolPatch:InterventionFinal
) |> tidytable::filter(
  TimeSinceIntervention < 0 | # time recorded before the intervention
    abs(TimeSinceIntervention - 20) < 1e-6 | # time 20 CTUs after.
    round((TimeSinceIntervention - 200)/10) == 0 | # time ~200 CTUs after
    round((TimeSinceIntervention - 500)/10) == 0 | # time ~500 CTUs after
    round((TimeSinceIntervention - 900)/10) == 0  # time ~900 CTUs after
) |> tidytable::arrange(
  TimeSinceIntervention
) |> tidytable::mutate(
  TimeLabel = tidytable::case_when(
    TimeSinceIntervention < 0 ~ "Pre-Intervention",
    dplyr::lag(TimeSinceIntervention) < 0 ~ "20 CTUs",
    dplyr::lag(TimeSinceIntervention, n = 2) < 0 ~ "~200 CTUs",
    dplyr::lag(TimeSinceIntervention, n = 3) < 0 ~ "~500 CTUs",
    dplyr::lag(TimeSinceIntervention, n = 4) < 0 ~ "~900 CTUs"
  ),
  TimeLabel = factor(TimeLabel, levels = c(
    "Pre-Intervention", "20 CTUs", "~200 CTUs", "~500 CTUs", "~900 CTUs"
  ), ordered = TRUE)
) |> tidytable::ungroup(
)

supplement8$graph$networks <- generateNetworks(supplement8$graph$specification)

supplement8$kdes <- lapply(
  supplement8$graph$networks$Envs, function(e) {
    e$trophics$EdgeVertexLists[[1]][[1]]$Vertices |> tidytable::select(
      node, Type, Size, N
    ) |> cbind(
      e$Row |> tidytable::select(Time, TimeLabel, PoolPatch:InterventionFinal)
    ) |> tidytable::mutate(
      AffinityVals = e$result$Ellipsis$Affinity$SpeciesAffinities[
        as.numeric(substring(node, 2))
        ]
    )
  }
) |> tidytable::bind_rows(
)

supplement8$PlotDensity <- ggplot2::ggplot(
) + ggplot2::geom_density(
  data = supplement8$kdes,
  mapping = ggplot2::aes(
    y = Size, fill = Type, color = Intervention
  ),
  trim = TRUE
) + ggplot2::facet_grid(
  Intervention + SpeciesPreferences ~ TimeLabel
) + ggplot2::scale_y_log10(
) + ggplot2::scale_color_manual(
  values = colorPalette,
  name = "Habitat Land-use"
) + ggplot2::scale_fill_manual(
  values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
) + ggplot2::theme_minimal(
)

ggplot2::ggsave(
  supplement8$PlotDensity,
  filename = paste0("Figure_supplement8_v1_", targetPrefIndex, "_Dens.pdf"),
  units = "cm", width = 6.5*3, height = 6.5*2
)

if (targetPref != "100% 0") {
  supplement8$PlotContour <-
    ggplot2::ggplot(
    ) + ggplot2::geom_density_2d(
      # ) + ggplot2::geom_bin_2d(
      data = supplement8$kdes,
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
      Intervention + SpeciesPreferences ~ TimeLabel
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
    )

  ggplot2::ggsave(
    supplement8$PlotContour,
    filename = paste0("Figure_supplement8_v1_", targetPrefIndex, "_Contour.pdf"),
    units = "cm", width = 10*3, height = 10*2
  )
}
