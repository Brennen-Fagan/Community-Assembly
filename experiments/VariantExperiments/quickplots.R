# Quick plotting
library(ggplot2)
library(dplyr)

quickfolder <- "TSTS_Simulations_2-1_30-30_2024-10-22"
quickbase <- "TSTS_Diversity_2-1-4-2-NA-7-21_30-30-53-2-12"
# experiments
# ppDO, dynDO, eDO, icDO, dispDO, distDO, aDO,
# dir(path = quickfolder, pattern = "Diversity_2-1-4-2-NA-7-21")


load(file.path(
  quickfolder, paste0(quickbase, ".RData")
))
Diversity$Diversity %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::mutate(Value = ifelse(is.na(Value), 0, Value)) -> baseRichness

load(file.path(
  quickfolder, paste0(quickbase, "_115-1-p-p_1-1.RData")
))
Diversity$Diversity %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::mutate(Value = ifelse(is.na(Value), 0, Value)) -> Richness1

load(file.path(
  quickfolder, paste0(quickbase, "_114-1-p-p_1-1.RData")
))
Diversity$Diversity %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::mutate(Value = ifelse(is.na(Value), 0, Value)) -> Richness0.75

load(file.path(
  quickfolder, paste0(quickbase, "_113-1-p-p_1-1.RData")
))
Diversity$Diversity %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::mutate(Value = ifelse(is.na(Value), 0, Value)) -> Richness0.5

load(file.path(
  quickfolder, paste0(quickbase, "_112-1-p-p_1-1.RData")
))
Diversity$Diversity %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::mutate(Value = ifelse(is.na(Value), 0, Value)) -> Richness0.25

load(file.path(
  quickfolder, paste0(quickbase, "_111-1-p-p_1-1.RData")
))
Diversity$Diversity %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::mutate(Value = ifelse(is.na(Value), 0, Value)) -> Richness0

ggplot(
  dplyr::bind_rows(
    Richness0 %>% dplyr::mutate(patchAffinity = 0),
    Richness0.25 %>% dplyr::mutate(patchAffinity = 0.25),
    Richness0.5 %>% dplyr::mutate(patchAffinity = 0.5),
    Richness0.75 %>% dplyr::mutate(patchAffinity = 0.75),
    Richness1 %>% dplyr::mutate(patchAffinity = 1)
  ) %>% dplyr::mutate(Subset = ifelse(is.na(Subset), "Total", Subset)),
  aes(x = Time, y = Value, color = factor(patchAffinity))
) + geom_line(
  data = baseRichness %>% dplyr::mutate(
    Subset = ifelse(is.na(Subset), "Total", Subset)
  ), inherit.aes = FALSE,
  mapping = aes(x = Time, y = Value), color = "black"
) + geom_line() + ggplot2::facet_grid(
  factor(Subset, ordered = TRUE,
         levels = rev(levels(factor(Subset)))) ~ patchAffinity
) + ggplot2::labs(
  y = "Richness",
  color = "Final\nPatch\nValue"
) + ggplot2::coord_cartesian(
  expand = FALSE
)
