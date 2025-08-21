addDistanceColumns <- function(bootstrapSamples, mindig = 1, Time = "Time") {
  bootstrapSamples %>% dplyr::mutate(
    ##### Temporal Distance: ##################################################
    TimeSinceStart = round(# Rare 1 != 1 issue.
      TimeSinceStart,
      digits = mindig
    )
  ) %>% dplyr::group_by(
    Type, Control, Bootstrap, Patch
  ) %>% dplyr::arrange(
    TimeSinceStart
  ) %>% dplyr::mutate(
    TimeGapNumber = seq_along(TimeSinceStart),
    TimeGapNumber = ifelse(
      Type == "Time series" & Control == "Control",
      rev(TimeGapNumber), TimeGapNumber
    ) # i.e.  5, 4, 3, 2, 1, ___, 1, 2, 3, 4, 5
  ) %>% dplyr::ungroup(
    ##### Spatial Distance: #################################################
  ) %>% dplyr::group_by(
    Type, Control, Bootstrap, TimeSinceStart
  ) %>% dplyr::mutate(
    DistanceFromCenter = dplyr::case_when(
      min(Patch) == 1 & max(Patch) == 10 ~ {
        ifelse(Patch < 5, Patch + 10, Patch) - median(
          ifelse(Patch < 5, Patch + 10, Patch)
        )
      },
      TRUE ~ Patch - median(Patch)
    ),
    DistanceFromCenterExpRev = ifelse(Control == "Experiment" &
                                        Type == "Space for time",
                                      -DistanceFromCenter,
                                      DistanceFromCenter)
  ) %>% dplyr::ungroup(
    ##### Define control species: ###########################################
  ) %>% dplyr::group_by(
    DistanceFromCenterExpRev, Bootstrap, TimeGapNumber, Type
  ) %>% dplyr::group_modify(
    .f = ~ computeSpeciesInControl(.x, Time = Time)
  ) %>% dplyr::ungroup(
  )
}
