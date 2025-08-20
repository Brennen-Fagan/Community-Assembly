flattenDiversity <- function(d) {
  if ("FullID" %in% names(d$Ellipsis)) {
    id <- strsplit(
      strsplit(d$Ellipsis$FullID, "_", fixed = TRUE)[[1]], # Split seeds off.
      "-", fixed = TRUE)
  } else if ("ID" %in% names(d$Ellipsis)) {
    id <- strsplit(
      strsplit(d$Ellipsis$ID, "_", fixed = TRUE)[[1]], # Split seeds off.
      "-", fixed = TRUE)
  } else {
    id <- strsplit(
      strsplit(
        strsplit(d$File, ".", fixed = TRUE)[[1]][1], # Remove .RData.
        "_", fixed = TRUE)[[1]][3:4], # Remove TSTS_Type and split seeds off.
      "-", fixed = TRUE # Separate out the id values.
    )
  }
  
  if (length(id) < 3) {
    # I.e., no intervention.
    id[[3]] <- rep(NA, 4)
    id[[4]] <- rep(NA, 2)
  }
  
  pres <- tidytable::data.table(d$Presence) %>% tidytable::mutate(
    Environment1 = Environment,
    Environment2 = NA
  ) %>% tidytable::ungroup() %>% tidytable::group_by(
    Time, Environment1, Environment2
  ) %>% tidytable::summarise(
    `Average Size:0` = mean(Size),
    `Average Size:1` = sum(Size*Abundance)/sum(Abundance),
    `St.Dev. Size:0` = sqrt(var(Size)),
    `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance))),
    `Average LSize:0` = mean(log10(Size)),
    `Average LSize:1` = sum(log10(Size)*Abundance)/sum(Abundance),
    `St.Dev. LSize:0` = sqrt(var(log10(Size))),
    `St.Dev. LSize:1` = sqrt(var(log10(Size)*Abundance/sum(Abundance))),
    `Ratio Con/Bas:0` = sum(Type == "Consumer")/sum(Type == "Basal"),
    `Ratio Con/Bas:1` = sum((Type == "Consumer") * Abundance) /
      sum((Type == "Basal") * Abundance),
    `Average Aff.:0` = mean(Affinity),
    `Average Aff.:1` = sum(Affinity*Abundance)/sum(Abundance)
  ) %>% tidytable::pivot_longer(
    cols = `Average Size:0`:`Average Aff.:1`,
    names_to = "Metric", values_to = "Value"
  ) %>% tidytable::mutate(
    Subset = NA
  )
  
  pres_subset <- tidytable::data.table(d$Presence) %>% tidytable::mutate(
    Environment1 = Environment,
    Environment2 = NA,
    Subset = paste0(Type, "_", Affinity)
  ) %>% tidytable::ungroup() %>% tidytable::group_by(
    Time, Environment1, Environment2, Subset
  ) %>% tidytable::summarise(
    `Average Size:0` = mean(Size),
    `Average Size:1` = sum(Size*Abundance)/sum(Abundance),
    `St.Dev. Size:0` = sqrt(var(Size)),
    `St.Dev. Size:1` = sqrt(var(Size*Abundance/sum(Abundance))),
    `Average LSize:0` = mean(log10(Size)),
    `Average LSize:1` = sum(log10(Size)*Abundance)/sum(Abundance),
    `St.Dev. LSize:0` = sqrt(var(log10(Size))),
    `St.Dev. LSize:1` = sqrt(var(log10(Size)*Abundance/sum(Abundance)))
  ) %>% tidytable::pivot_longer(
    cols = `Average Size:0`:`St.Dev. LSize:1`,
    names_to = "Metric", values_to = "Value"
  )
  
  tidytable::data.table(d$Diversity) %>% tidytable::bind_rows(
    pres, pres_subset
  ) %>% tidytable::mutate(
    PoolPatch = id[[1]][1],
    PoolPatchSeed = id[[2]][1],
    Interactions = id[[1]][2],
    InteractionsSeed = id[[2]][2],
    Events = id[[1]][3],
    EventsSeed = id[[2]][3],
    InitialConditions = id[[1]][4],
    InitialConditionsSeed = id[[2]][4],
    Dispersal = id[[1]][5],
    NicheDistance = id[[1]][6],
    SpeciesAffinity = id[[1]][7],
    SpeciesAffinitySeed = id[[2]][5],
    PatchAffinity = id[[1]][8],
    PatchAffinitySeed = id[[2]][6],
    InterventionPatchType = id[[3]][1],
    InterventionPatchSeed = id[[4]][1],
    InterventionTimeType = id[[3]][2],
    InterventionTimeSeed = id[[4]][2],
    InterventionDispersal = id[[3]][3],
    InterventionNicheDistance = id[[3]][4]
  )
}