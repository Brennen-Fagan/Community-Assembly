# Setup: ######################################################################

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsTimeBC.R")
supplement12 <- list()

### 12 Supplement: ############################################################
##### Turnover related statistics: ############################################
# Problem here: data is not sufficiently normal, so central tendency isn't quite
# catching the right information.

supplement12$data <- diversitiesTimeBC |> tidytable::filter(
  PoolPatchSeed %in% basePoolPatchSeeds,
  NicheDistance == defaultNicheDistance,
  grepl(x = Metric, pattern = "TimeBrayCurtis:", fixed = TRUE),
  is.na(Subset),
  Time > Start, Time < Stop # Not things outside of [Start, Stop]
)

supplement12$plot <- ggplot2::ggplot(
  supplement12$data,
  ggplot2::aes(x = Value, color = Intervention)
) + ggplot2::geom_histogram(
  show.legend = FALSE, binwidth = 0.001
) + ggplot2::geom_boxplot(
  data = ~ .x |> tidytable::group_by(
    # Within run average
    Intervention, PoolPatchSeed,
    SpeciesPreferences, InterventionInitial, InterventionFinal
  ) |> tidytable::summarise(
    Value = mean(Value),
    # ) |> tidytable::group_by(
    #   # Between run average.
    #   Intervention
    # ) |> tidytable::summarise(
    #   Value = mean(Value),
    y = 100
  ),
  mapping = ggplot2::aes(y = y, group = Intervention),
  color = "black",
  show.legend = FALSE
) + ggplot2::facet_grid(
  SpeciesPreferences + InterventionInitial ~ InterventionFinal
) + ggplot2::scale_color_manual(
  values = colorPalette
) + ggplot2::scale_y_log10(
)

ggplot2::ggsave(
  supplement12$plot,
  filename = "Figure_supplement12_v1.pdf",
  units = "cm", width = 20*3, height = 20*2
)
