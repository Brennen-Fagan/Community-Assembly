load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/RMTRCode2/experiments/RMTR-MNA/Data_2022-08-05/MNA-ExampleOutcome-Result-Env10-None.RData")
load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/RMTRCode2/experiments/RMTR-MNA/Data_2022-08-05/MNA-ExampleOutcome-PoolMats-Env10.RData")
library(RMTRCode2)
library(dplyr)
RMTRCode2::CalculateTrophicStructure(
  Pool, result$NumEnvironments, InteractionMatrices,
  result$Parameters$EliminationThreshold
) -> trophicStructureCalculator

largestOccupancy <- which.max(
  apply(result$Abundance, 1, function(x) sum(x > 0))
)

trophicStructureCalculator(result$Abundance[largestOccupancy, -1]) -> temp

sizes <- unlist(lapply(temp$EdgeVertexLists, function(x) nrow(x$Edges)))

targets <- order(sizes)[10:8]

# We need (only?) the positive flow of biomass from each node.
# EffectActual is:
#    Intrisc+: Reproductive Rate * Abundance
#    Exploit+: EffectPerUnit (Interaction Matrix Value) * ResourceAbundance
# So we need to arrive at the total amount of biomass gained from these.
# I think the best way to do so is to multiply Intrisc+ by size, but multiply
# Exploit+ by ConsumerAbundance AND Consumer size.
# Furthermore, I think Exploit+ is the contents of the matrix, while
# Intrisc+ is the contents of the import row.
# I'm going to leave out exports and respiration.

for (targ in targets) {
  targEdges <- temp$EdgeVertexLists[[targ]]$Edges %>% dplyr::filter(
    effectSign > 0
  ) %>% dplyr::ungroup() %>% dplyr::select(
    from, to, effectActual, Type
  )

  targNodes <- temp$EdgeVertexLists[[targ]]$Vertices %>% dplyr::mutate(
    Biomass = Size * N
  ) %>% dplyr::select(-ReproductionRate, -Intraspecific, -Type)

  targEdges <- dplyr::left_join(
    targEdges, targNodes,
    by = c("to" = "node")
  ) %>% dplyr::mutate(
    BiomassFlow = effectActual * ifelse(
      Type == "Intrisc+", Size, Biomass
    ),
    # Column (To) Major order
    Index =
      (as.numeric(factor(to, levels = targNodes$node, ordered = TRUE)) - 1) *
      length(targNodes$node) + as.numeric(factor(from, levels = targNodes$node, ordered = TRUE))
  )

  targFlow <- targEdges %>% dplyr::filter(Type != "Intrisc+")
  targMprt <- targEdges %>% dplyr::filter(Type == "Intrisc+")

  importRow <- matrix(0, ncol = nrow(targNodes), nrow = 1)
  importRow[
    with(targMprt,
         as.numeric(factor(from, levels = targNodes$node, ordered = TRUE)))
    ] <- targMprt$BiomassFlow

  outputMatrix <- matrix(0,
                         nrow = nrow(targNodes),
                         ncol = nrow(targNodes))
  outputMatrix[targFlow$Index] <-
    targFlow$BiomassFlow

  colnames(outputMatrix) <- targNodes$node
  rownames(outputMatrix) <- targNodes$node

  outputMatrix <- cbind(
    outputMatrix,
    "IsAlive" = 1, "Biomass" = targNodes$Biomass,
    "TrophicLevel" = Pool$Size[
      paste0("s", Pool$ID) %in% rownames(temp$TrophicLevels[[targ]])
      ]*10
    )

  outputMatrix <- rbind(
    outputMatrix,
    "Import" = c(importRow, rep(0, 3))
  )

  outputMatrix <- cbind(
    outputMatrix, "Export" = 0, "Respiration" = 0
    )

  write.table(
    tibble::rownames_to_column(data.frame(outputMatrix), var = "Names"),
    file = paste0("community-", targ, ".csv"),
    sep = ";", quote = FALSE,
    col.names = TRUE)
}

