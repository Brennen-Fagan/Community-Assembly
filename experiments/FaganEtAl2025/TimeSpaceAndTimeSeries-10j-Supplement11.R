# Setup: ######################################################################
# Capture the ways in which species enter and exit the basic set of systems.
# These are the "in/out" statistics, according to whether there was a
# successful/failed immigration attempt, a neutral stochastic extirpation, or an
# extirpation due to species dynamics. We keep counts and percentages and see
# how they vary as the land-use and land-use preferences vary.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")
supplement11 <- list()

### 11 Supplement: #############################################################
##### In/Out Statistics: ######################################################

supplement11$data <- tidytable::bind_rows(
  Pers |> tidytable::filter(
    # SpeciesPreferences == "100% 0",
    NicheDistance == defaultNicheDistance,
    # Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Persistence > 0,
    InType != externalNames["Dispersal"],
    In < Stop, Out > Start
  ) |> tidytable::group_by(
    SpeciesPreferences, InType, OutType, Intervention
  ) |> tidytable::summarise(
    Average = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
  ),
  ColExt |> tidytable::filter(
    # Filter what we can:
    NicheDistance == defaultNicheDistance,
    PoolPatchSeed %in% basePoolPatchSeeds,
    !Success | EventType == "Present"
  ) |> tidytable::left_join(
    # Start and Stop aren't already present in this version
    endTimes |> tidytable::select(-Times)
  ) |> tidytable::filter(
    Times > Start, Times < Stop
  ) |> tidytable::mutate(
    SpeciesPreferences =
      speciesAffinityDictionaryOrigin$SpeciesAffinities[as.numeric(SpeciesAffinity)]
  ) |> changePreferencesLevels(
    # ) |> tidytable::filter(
    #   SpeciesPreferences == "100% 0"
  ) |> tidytable::left_join(
    interventionStrings,
    by = c("PatchAffinity", "PoolPatch", "InterventionPatchType")
    # ) |> tidytable::filter(
    #   Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)")
  ) |> tidytable::mutate(
    InType = externalNames[
      ifelse(EventType == "Arrival", "Failed Arrival", "Present")
      ],
    OutType = externalNames["NA"]
  ) |> tidytable::group_by(
    SpeciesPreferences, InType, OutType, Intervention
  ) |> tidytable::summarise(
    Average = tidytable::n() / tidytable::n_distinct(PoolPatchSeed)
  )
) |> changeInterventionLevels(
)

supplement11$labels <- supplement11$data |> tidytable::group_by(
  # InType, OutType,
  SpeciesPreferences, InterventionInitial, InterventionFinal
) |> tidytable::mutate(
  Percentage = Average / sum(Average), # /44's cancel out.
  Text = paste0(formatC(Percentage*100, digits = 1, format = "f"), "%"),
  y = 400
)

supplement11$plot <- ggplot2::ggplot(
  supplement11$data,
  ggplot2::aes(x = interaction(InType,
                               OutType,
                               sep = "\n", lex.order = TRUE),#"",
               y = Average,
               fill = interaction(InType,
                                  OutType,
                                  sep = "/"))
) + ggplot2::geom_col(
  show.legend = FALSE
) + ggplot2::geom_text(
  data = supplement11$labels,
  mapping = ggplot2::aes(y = y, label = Text),
  show.legend = FALSE
) + ggplot2::facet_grid(
  SpeciesPreferences + InterventionInitial ~ InterventionFinal
) + ggplot2::theme_minimal(
) + ggplot2::labs(
  subtitle =
    "Columns = Final Land-use, Rows = Initial Land-use and Species Type",
  x = "Enter/Exit Type", y = "Average Count"#, fill = "In/Out"
)

ggplot2::ggsave(
  supplement11$plot,
  filename = "Figure_supplement11_v1.png",
  units = "cm", width = 20*3, height = 20*2
)
