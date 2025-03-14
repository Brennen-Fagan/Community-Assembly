# double checks
library(dplyr)
source('/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TimeSpaceAndTimeSeries-0-Functions.R', echo=TRUE)
# load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_345-345_2025-01-21/TSTS_Intervention_142486-4929-28-1-NA-5-21_345-345-388-391-469_111-1-p-p_1-1.RData")
# load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_369-369_2025-01-23/TSTS_Intervention_142486-4929-28-1-NA-5-21_369-369-412-415-829_112-1-p-p_1-1.RData")
# load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_382-382_2025-01-24/TSTS_Intervention_142486-4929-28-1-NA-5-6_382-382-425-428-1017_113-1-p-p_1-1.RData")
# load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_342-342_2025-01-27/TSTS_Intervention_142486-4929-28-1-NA-7-28_342-342-385-388-426_115-1-p-p_1-1.RData")
load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/TSTS_Simulations_142486-4929_364-364_2025-01-22/TSTS_Simulation_142486-4929-28-1-NA-7-14_364-364-407-410-751.RData")
targets <- result$Events %>% dplyr::filter(Success) %>% dplyr::arrange(Times)

# abund <- result$Abundance[(t(outer(which(result$Abundance[, 1] %in% targets$Times), -1:1, "+"))),]
#
# stopifnot(nrow(abund)/3 == nrow(targets))

output <- calculateColExtMetrics(result)

# What are we looking for:
#   These are all the successful events, so they should all be caught properly.
# stackoverflow.com/a/22674231
targetcheck3 <- do.call(paste0, targets[, 1:3])
targetcheck5 <- do.call(paste0, targets[, 1:5])
# Not including success (trivial) or event type (remapping of types, esp.
# some Arrival -> Present).
stopifnot(all(targetcheck3 %in% do.call(paste0, output[, 1:3])))


# We should also be able to check the list of events provided and make sure
# they make sense with the actual abundance.
# Note this ignores the events that come into conflict, but should clear most.
output$Doublecheck <- NA
for (row in 1:nrow(output)) {
  outrow <- output[row, ]
  abundrownums <- which(outrow$Times == result$Abundance[, 1]) + -1:1
  abundrows <- result$Abundance[abundrownums, -1]

  if (outrow$Type == "Arrival") {
    output[row, ]$Doublecheck <- (
      abundrows[2, outrow$Species] < result$Parameters$EliminationThreshold &&
        abundrows[3, outrow$Species] > result$Parameters$EliminationThreshold &&
        do.call(paste0, outrow[, 1:5]) %in% targetcheck5
    )
  } else if (outrow$Type == "Present") {
    output[row, ]$Doublecheck <- (
      abundrows[1, outrow$Species] > result$Parameters$EliminationThreshold &&
        abundrows[2, outrow$Species] > result$Parameters$EliminationThreshold &&
        !(do.call(paste0, outrow[, 1:5]) %in% targetcheck5)
    )
  } else if (outrow$Type == "Extinct") {
    output[row, ]$Doublecheck <- (
      abundrows[2, outrow$Species] > result$Parameters$EliminationThreshold &&
        abundrows[3, outrow$Species] < result$Parameters$EliminationThreshold &&
        do.call(paste0, outrow[, 1:5]) %in% targetcheck5
    )
  } else if (outrow$Type == "Dynamic Loss") {
    output[row, ]$Doublecheck <- (
      abundrows[2, outrow$Species] > result$Parameters$EliminationThreshold &&
        abundrows[3, outrow$Species] < result$Parameters$EliminationThreshold &&
        !(do.call(paste0, outrow[, 1:5]) %in% targetcheck5)
    )
  }
}

stopifnot(all(output$Doublecheck))

output <- output %>% dplyr::mutate(Sign = dplyr::case_when(
  Type == "Arrival" ~ 1,
  Type == "Present" ~ 0,
  Type %in% c("Extinct", "Dynamic Loss") ~ -1
))

temp <-
  (result$Abundance[-1, -1] > result$Parameters$EliminationThreshold) -
  (result$Abundance[-nrow(result$Abundance), -1] > result$Parameters$EliminationThreshold)

changeCoords <- rbind(
  data.frame(
    Species = col(temp)[which(temp == 1)],
    Environment = 1,
    Times = result$Abundance[row(temp)[which(temp == 1)], 1],
    Type = 1
  ),data.frame(
    Species = col(temp)[which(temp == -1)],
    Environment = 1,
    Times = result$Abundance[row(temp)[which(temp == -1)], 1],
    Type = -1
  )
)

output %>% dplyr::filter(Sign != 0) %>% dplyr::select(Species, Environment, Times, Sign)

stopifnot(all(do.call(paste0, changeCoords) %in%
                do.call(paste0, output %>% dplyr::filter(Sign != 0) %>% dplyr::select(Species, Environment, Times, Sign))))

# To summarise: all events should be in output, and all changes should be in
#          output, and all records in output should match the abundance changes.
