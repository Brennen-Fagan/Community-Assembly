
source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")

figureS9 <- list(
  # pref = "100% 0",
  # pref = "50% 0, 50% 1",
  pref = "Uniform(0, 1)",
  luinitl = c("(0)", "(0.5)", "(1)"), # Land Use INITiaL
  lufinal = c("(0)", "(0.5)", "(1)") # Land Use FINAL
  # luinitl = c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  # lufinal = c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)")
)

figureS9$prefstring <- switch(
  figureS9$pref,
  "100% 0" = "1000",
  "50% 0, 50% 1" = "5050",
  "Uniform(0, 1)" = "Unif"
)

figureS9$lustring <- paste0(
  length(figureS9$luinitl),"to",length(figureS9$lufinal)
)

# Main Plots: #################################################################
### Plot 9: ###################################################################
##### Data: ###################################################################

figureS9$interventionTimes <- diversitiesRichness |> tidytable::filter(
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

figureS9$data <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences == figureS9$pref,
  InterventionInitial %in% figureS9$luinitl,
  InterventionFinal %in% figureS9$lufinal,
  NicheDistance == defaultNicheDistance,
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0" ,
  is.na(Subset)
) |> tidytable::left_join(
  figureS9$interventionTimes |> tidytable::rename(
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

figureS9$dataSummary <- figureS9$data |> tidytable::group_by(
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
  Total == length(unique(figureS9$data$PoolPatchSeed)) # Num. of sims!
) |> tidytable::pivot_longer(
  cols = Neg:Pos, names_to = "Type", values_to = "Counts"
) |> tidytable::mutate(
  Percentage = Counts / Total * 100,
  Text = paste0(formatC(Percentage, digits = 1, format = "f"), "%"),
  Comp = paste0(formatC(100-Percentage, digits = 1, format = "f"), "%")
)

##### Main Plot: ##############################################################
figureS9$plot <- ggplot2::ggplot(
  figureS9$dataSummary |> tidytable::filter(
    Type == "Neg"
  ),
  ggplot2::aes(x = TimeSinceIntervention, color = Intervention,
               y = (100 - Percentage)/100, # for Complement
               group = interaction(Type, Intervention, SpeciesPreferences))
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
) + ggplot2::scale_x_continuous(
  transform = "log1p", breaks = c(0, 1, 10, 100, 1000, 10000), expand = c(0, 0)
) + ggplot2::scale_y_continuous(
  labels = scales::percent
) + ggplot2::labs(
  x = "Time Since Intervention",
  y = "Intact Ecosystems"
)

##### Save: ###################################################################

figureS9$suffix <- paste0("_", figureS9$prefstring, "_", figureS9$lustring)
figureS9$prefix <- "FigureS9"
figureS9$iter <- "_Prototype1"
figureS9$ext <- c(".png", ".pdf")

with(figureS9,
     ggplot2::ggsave(
       plot = plot,
       filename = file.path(dirImages, paste0(
         prefix, iter, suffix, ext[1]
       )),
       units = "cm", width = 6.5*3, height = 6.5*2)
)
with(figureS9,
     ggplot2::ggsave(
       plot = plot,
       filename = file.path(dirImages, paste0(
         prefix, iter, suffix, ext[2]
       )),
       units = "cm", width = 6.5*3, height = 6.5*2)
)