
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

figure8 <- list(
  # pref = "100% 0",
  # pref = "50% 0, 50% 1",
  pref = "Uniform(0, 1)",
  luinitl = c("(0)", "(0.5)", "(1)"), # Land Use INITiaL
  lufinal = c("(0)", "(0.5)", "(1)") # Land Use FINAL
  # luinitl = c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  # lufinal = c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)")
)

figure8$prefstring <- switch(
  figure8$pref,
  "100% 0" = "1000",
  "50% 0, 50% 1" = "5050",
  "Uniform(0, 1)" = "Unif"
)

figure8$lustring <- paste0(
  length(figure8$luinitl),"to",length(figure8$lufinal)
)

# Main Plots: #################################################################
### Plot 9: ###################################################################
##### Data: ###################################################################

figure8$interventionTimes <- diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0",
  is.na(Subset), # Won't matter, so less data
  InterventionInitial != InterventionFinal
) |> tidytable::select(
  PoolPatch, PoolPatchSeed, Time
) |> tidytable::group_by(
  PoolPatch, PoolPatchSeed
) |> tidytable::summarise(
  Time = min(
    Time[round(Time, 6)!=round(min(Time), 6)]
  ), # Not min(Time) but the next time.
  .groups = "drop"
)

figure8$data <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == figure8$pref,
  InterventionInitial %in% figure8$luinitl,
  InterventionFinal %in% figure8$lufinal,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0" ,
  is.na(Subset)
) |> tidytable::left_join(
  figure8$interventionTimes |> tidytable::rename(
    InterventionTime = Time
  ),
  by = c("PoolPatch", "PoolPatchSeed")
) |> tidytable::filter(
  # Make sure we get from InterventionTime forwards, but note num. error.
  round(Time - InterventionTime, digits = 10) >= 0
) |> tidytable::mutate(
  Metric = factor(Metric, levels = c("Alpha Hill:0", "Alpha Abundance"),
                  labels = c("Richness", "Abundance"), ordered = TRUE),
  TimeSinceIntervention = round(Time - InterventionTime, digits = 10),
  InterventionChange = abs( # Size of Intervention, 0, 0.25, 0.5, 0.75, or 1.
    as.numeric(gsub(InterventionInitial, pattern = "[(]|[)]", replacement = ""))
    - as.numeric(gsub(InterventionFinal, pattern = "[(]|[)]", replacement = ""))
  ),
  TimeSinceIntervention = tidytable::case_when( # Create groupings for times.
    TimeSinceIntervention <= 50 ~ round(TimeSinceIntervention, 0),
    TimeSinceIntervention < 1105 ~ round(TimeSinceIntervention, -1), # Skip breaks < 5, drop.
    TimeSinceIntervention < 16350 ~ round(TimeSinceIntervention, -2),
    TRUE ~ TimeSinceIntervention
  )
) |> tidytable::group_by(
  # Average Over the now grouped times to make each sim equally weighted.
  # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
  # simple and match each other with the seeds. More complicated set-ups
  # will want to adjust the groupings here.
  Intervention, InterventionInitial, InterventionFinal, Metric,
  PoolPatchSeed, SpeciesPreferences, InterventionChange, TimeSinceIntervention
) |> tidytable::summarise(
  # Make sure there's only one value...
  Value = median(Value),
  .groups = "drop_last"
) |> tidytable::mutate(
  # Preserving most of the grouping structure
  ValueDifference = Value - Value[which(TimeSinceIntervention == 0)]
  # Note TimeSinceIntervention == 0 at round(digits = 10) scale per the filter.
)

# ggplot2::ggplot(
#   figure8$data,
#   ggplot2::aes(
#     x = TimeSinceIntervention, y = Value,
#     color = Intervention, group = interaction(Intervention, PoolPatchSeed))
# ) + ggplot2::geom_line(alpha = 0.1) + ggplot2::facet_grid(
#   InterventionInitial ~ InterventionFinal
# )

figure8$dataSummary <- figure8$data |> tidytable::group_by(
  # Summarise across PoolPatchSeed
  Intervention, InterventionInitial, InterventionFinal, Metric,
  SpeciesPreferences, InterventionChange, TimeSinceIntervention
) |> tidytable::summarise(
  # Across PoolPatchSeeds
  Total = tidytable::n(),
  Neg = sum(ValueDifference < 0),
  Zero = sum(ValueDifference == 0),
  Pos = sum(ValueDifference > 0)
) |> tidytable::filter(
  Total == length(unique(figure8$data$PoolPatchSeed)) # Num. of sims!
) |> tidytable::pivot_longer(
  cols = Neg:Pos, names_to = "Type", values_to = "Counts"
) |> tidytable::mutate(
  Percentage = Counts / Total * 100,
  Text = paste0(formatC(Percentage, digits = 1, format = "f"), "%"),
  Comp = paste0(formatC(100-Percentage, digits = 1, format = "f"), "%")
)

##### Main Plot: ##############################################################
figure8$plot <- ggplot2::ggplot(
  figure8$dataSummary |> tidytable::filter(
    Type == "Neg"
  ),
  ggplot2::aes(x = TimeSinceIntervention, color = Intervention,
               y = (100 - Percentage)/100, # for Complement
               group = interaction(Type, Intervention, SpeciesPreferences))
# ) + ggplot2::geom_vline(
#   xintercept = c(0, 10, 10000)
) + ggplot2::geom_text(
  data = data.frame(text = c("Short (t = 10)", "Long (t = 10000)"),
                    x = c(10, 10000), y = 0.25),
  ggplot2::aes(label = text, x = x, y = y), angle = 90,
  inherit.aes = FALSE
) + ggplot2::geom_line(
) + ggplot2::geom_point(
  alpha = 0.2, shape = 21 # open circles
) + ggplot2::scale_color_manual(
  values = colorPalette, aesthetics = c("color", "fill"),
  name = "Habitat Type"
) + ggplot2::theme_minimal(
# ) + ggplot2::facet_grid(
#   InterventionInitial ~ InterventionFinal
) + ggplot2::scale_x_continuous(
  transform = "log1p", breaks = c(0, 1, 10, 100, 1000, 10000), expand = c(0, 0)
) + ggplot2::scale_y_continuous(
  labels = scales::percent
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Intact Ecosystems"#expression(Delta ~ "Richness" >= 0)
)

##### Save: ###################################################################

figure8$suffix <- paste0("_", figure8$prefstring, "_", figure8$lustring)
figure8$prefix <- "Figure8"
figure8$iter <- "_Prototype1"
figure8$ext <- c(".png", ".pdf")

with(figure8,
     ggplot2::ggsave(
       plot = plot,
       filename = file.path(dirImages, paste0(
         prefix, iter, suffix, ext[1]
       )),
       units = "cm", width = 6.5*3, height = 6.5*2)
)
with(figure8,
     ggplot2::ggsave(
       plot = plot,
       filename = file.path(dirImages, paste0(
         prefix, iter, suffix, ext[2]
       )),
       units = "cm", width = 6.5*3, height = 6.5*2)
)

#
# # Take some slices so we end up with something more writable in a paper.
# supplementStatistics$STAT$shortTermLoss <-
#   # 3 times for each of the 3 preferences with each of the 3 statistics.
#   supplementStatistics$shortTermLoss df_ShortTermChanges |> tidytable::filter(
#     TimeSinceIntervention %in% c(1, 10, 20, 50)
#   ) |> tidytable::group_by(
#     TimeSinceIntervention
#   ) |> tidytable::summarise(
#     Total = sum(Total),
#     Neg = sum(Neg),
#     Zero = sum(Zero) ,
#     Pos = sum(Pos)
#   )

# # A tidytable: 5 × 5
# TimeSinceIntervention Total   Neg  Zero   Pos
# <dbl> <int> <int> <int> <int>
#   1                     1  2640  1108  1471    61
# 2                     5  2640  1963   592    85
# 3                    10  2640  2148   394    98
# 4                    20  2640  2178   281   181
# 5                    50  2640  2091   250   299
