# Setup: ######################################################################
# Show how abundance/biomass change. This might benefit from paths for each
# poolpatchseed following the interventionFinal from (0) to (1). Maybe a 1:1
# line.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

# This is better as an environment, but that's more opaque.
supplement15 <- list()
supplement15$preferences <- 2 # 1 Abundance 2 Biomass

# Main Plots: #################################################################
### Plot 16:###################################################################
supplement15$data <- diversitiesAbund |> tidytable::filter(
  Metric == switch(supplement15$preferences,
                   "Alpha Abundance", # TOTAL Abundance
                   "Alpha Biomass" # TOTAL Biomass
  ),
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Start <= Time, Time <= Stop,
  !is.na(Subset)
) |> tidytable::separate_wider_delim(
  Subset, delim = "_", names = c("Guild", "AffinityBins")
) |> unifyAffinityBins(
)

supplement15$yMax <- supplement15$data$Value |> max()

supplement15$data <- supplement15$data |> tidytable::select(
  -Left, -Right
) |> tidytable::pivot_wider(
  names_from = Guild, values_from = Value
)

supplement15$dataAverage <- supplement15$data |> tidytable::group_by(
  Environment1:AffinityBins
) |> tidytable::mutate(
  Basal = mean(Basal),
  Consumer = mean(Consumer)
)

##### a: ######################################################################
####### Core Plot: ############################################################
supplement15$plotAverage <- ggplot2::ggplot(
  supplement15$dataAverage,
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = Intervention,
    shape = AffinityBins
  )
) + ggplot2::geom_point(
  show.legend = FALSE
)  + ggplot2::facet_grid(
  SpeciesPreferences + InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_log10(
) + ggplot2::scale_shape_manual(
  values = c("a", "e", "i", "o", "u", "y", "£")
) + ggplot2::labs(
  x = switch(supplement15$preferences,
             "Average Total Basal Abundance", # TOTAL Abundance
             "Average Total Basal Biomass" # TOTAL Biomass
  ),
  y = switch(supplement15$preferences,
             "Average Total Consumer Abundance", # TOTAL Abundance
             "Average Total Consumer Biomass" # TOTAL Biomass
  )
)

supplement15$plotAveragePath <- ggplot2::ggplot(
  supplement15$dataAverage |> tidytable::arrange(InterventionFinal),
  ggplot2::aes(
    x = Basal,
    y = Consumer,
    color = InterventionInitial,
    group = interaction(PoolPatchSeed, InterventionInitial,
                        AffinityBins, SpeciesPreferences),
    linetype = AffinityBins
  )
) + ggplot2::geom_path(
  alpha = 0.01
) + ggplot2::geom_path(
  data = function(x) x |> tidytable::filter(PoolPatchSeed == 1),
  alpha = 1
)  + ggplot2::facet_grid(
  SpeciesPreferences ~ .
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_y_log10(
) + ggplot2::scale_x_log10(
) + ggplot2::labs(
  x = switch(supplement15$preferences,
             "Average Total Basal Abundance", # TOTAL Abundance
             "Average Total Basal Biomass" # TOTAL Biomass
  ),
  y = switch(supplement15$preferences,
             "Average Total Consumer Abundance", # TOTAL Abundance
             "Average Total Consumer Biomass" # TOTAL Biomass
  )
) + ggplot2::guides(color = "none")

supplement15$plotConsumer <- ggplot2::ggplot(
  supplement15$data |> tidytable::arrange(Time),
  ggplot2::aes(
    x = Time,
    y = Consumer,
    color = Intervention,
    linetype = AffinityBins,
    group = interaction(PoolPatchSeed, AffinityBins,
                        Intervention, SpeciesPreferences)
  )
) + ggplot2::geom_line(
  show.legend = FALSE, alpha = 0.01
) + ggplot2::geom_line(
  data = function(x) x |> tidytable::filter(PoolPatchSeed == 1),
  show.legend = FALSE, alpha = 1
)  + ggplot2::facet_grid(
  SpeciesPreferences + InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_y_log10(
) + ggplot2::labs(
  y = switch(supplement15$preferences,
             "Average Total Consumer Abundance", # TOTAL Abundance
             "Average Total Consumer Biomass" # TOTAL Biomass
  )
)

supplement15$plotBasal <- ggplot2::ggplot(
  supplement15$data |> tidytable::arrange(Time),
  ggplot2::aes(
    x = Time,
    y = Basal,
    color = Intervention,
    linetype = AffinityBins,
    group = interaction(PoolPatchSeed, AffinityBins,
                        Intervention, SpeciesPreferences)
  )
) + ggplot2::geom_line(
  show.legend = FALSE, alpha = 0.01
) + ggplot2::geom_line(
  data = function(x) x |> tidytable::filter(PoolPatchSeed == 1),
  show.legend = FALSE, alpha = 1
)  + ggplot2::facet_grid(
  SpeciesPreferences + InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_y_log10(
) + ggplot2::labs(
    y = switch(supplement15$preferences,
               "Average Total Basal Abundance", # TOTAL Abundance
               "Average Total Basal Biomass" # TOTAL Biomass
    )
)

ggplot2::ggsave(plot = supplement15$plotAverage,
                filename = paste0("Figure_supplement15_v1A_",
                                  supplement15$preferences,
                                  ".pdf"),
                units = "cm", width = 6.5*9, height = 6.5*9)

ggplot2::ggsave(plot = supplement15$plotBasal,
                filename = paste0("Figure_supplement15_v1B_",
                                  supplement15$preferences,
                                  ".pdf"),
                units = "cm", width = 6.5*9, height = 6.5*9)

ggplot2::ggsave(plot = supplement15$plotConsumer,
                filename = paste0("Figure_supplement15_v1C_",
                                  supplement15$preferences,
                                  ".pdf"),
                units = "cm", width = 6.5*9, height = 6.5*9)
