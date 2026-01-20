diversitiesRichness |> tidytable::filter(
  PoolPatchSeed %in% basePoolPatchSeeds,
  NicheDistance == defaultNicheDistance,
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::ungroup(
) |> tidytable::left_join(
  InterventionTimes |> tidytable::select(
    TimeIntervention:PatchAffinitySeed
  )
) |> tidytable::filter(
  !is.na(TimeIntervention)
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal,
  PoolPatchSeed
) |> tidytable::arrange(
  Time
) |> tidytable::mutate(
  PostIntervention = Time > TimeIntervention,
  ApproxIntervention = tidytable::lead(PostIntervention) & !PostIntervention 
) |> tidytable::filter(
  ApproxIntervention | PostIntervention
) |> tidytable::group_by(
  SpeciesPreferences, Intervention, InterventionInitial, InterventionFinal,
  PoolPatchSeed
) |> tidytable::mutate(
  InterventionChange = abs(
    as.numeric(gsub(InterventionInitial, pattern = "[(]|[)]", replacement = ""))
    - as.numeric(gsub(InterventionFinal, pattern = "[(]|[)]", replacement = ""))
  ),
  TimeSinceIntervention = round(Time - Time[1], digits = 4), # Make numerically safe.
  ValueMinusIntervention = round(Value - Value[1]), # Make numerically safe.
  # !xor(TimeSinceIntervention == 0, ApproxIntervention) # should be all TRUE
) |> tidytable::filter(
  TimeSinceIntervention <= 51,
  TimeSinceIntervention == floor(TimeSinceIntervention)
) |> View()
