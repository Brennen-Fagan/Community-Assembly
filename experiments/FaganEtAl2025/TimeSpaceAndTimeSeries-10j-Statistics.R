diversitiesRichness |> tidytable::filter(
  NicheDistance == defaultNicheDistance,
  (PoolPatchSeed %in% as.character(343:386)),
  Metric == "Alpha Hill:0",
  is.na(Subset)
) |> tidytable::left_join(
  endTimes |> dplyr::select(-Times)
) |> tidytable::group_by(
  SpeciesAffinity, Intervention, PoolPatchSeed
) |> tidytable::arrange(
  Time
) |> tidytable::filter(
  InterventionInitial != InterventionFinal,
  Time == Time[1] | Time == Time[2]
) |> tidytable::summarise(
  InterventionChange = abs(
    as.numeric(gsub(InterventionInitial, pattern = "[(]|[)]", replacement = ""))
    - as.numeric(gsub(InterventionFinal, pattern = "[(]|[)]", replacement = ""))
  ),
  Time = Time[2] - Time[1],
  Value = Value[2] - Value[1],
  Method = "Temporal",
  .groups = "drop"
) |> with(table(Intervention, sign(Value), SpeciesAffinity))