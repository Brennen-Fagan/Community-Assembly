# Libraries: ##################################################################
library(RMTRCode2)
library(dplyr)
library(Matrix)
# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Dictionaries.R"))

load(file.path(directory,
               "CompareEliminationThresholds",
               "TSTS_PoolPatchDynamics_2-1_30-30.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "TSTS_Simulation_2-1-4-1-NA-1-1_30-30-53-1-1.RData"))

PatchAffinities <- result$Ellipsis$Affinity$PatchAffinities

# Instantiate PerCapitaDynamics: ############################################
# We've already built the "functional" interactions matrix, but we now need
# to apply the transformation r' = r rho^(sign(r)) and then combine.

# Not a function, so we don't need to on the fly.
# We'll be a bit lazy here, but hopefully readable and clear.
grid <- expand.grid(
  pool = 1:nrow(Pool), # Fastest Varying
  patch = 1:length(PatchAffinities) # Slower Varying.
)
rprime <- result$Ellipsis$Affinity$EffectiveReproductionRate

if (is.function(rprime)) {
  # Calculate rprime using Parms$Patch
  if (is.function(InteractionMatrices$Mats[[1]])) {
    # Calculate and combine interaction matrices on the fly.
    PerCapitaDynamics <- DynamicsFunction(
      rprime,
      function(t, y, parms) {
        Matrix::bdiag(lapply(
          InteractionMatrices$Mats,
          function(matfunc) {matfunc(t, y, parms)}
        ))
      },
      result$NumEnvironments
    )
  }
  else {
    # Just combine the interaction matrices.
    PerCapitaDynamics <- DynamicsFunction(
      rprime,
      Matrix::bdiag(InteractionMatrices$Mats),
      result$NumEnvironments
    )
  }
} else {
  # Treat rprime as constant and explicitly calculated.
  if (is.function(InteractionMatrices$Mats[[1]])) {
    # Calculate and combine interaction matrices on the fly.
    PerCapitaDynamics <- DynamicsFunction(
      rprime,
      function(t, y, parms) {
        Matrix::bdiag(lapply(
          InteractionMatrices$Mats,
          function(matfunc) {matfunc(t, y, parms)}
        ))
      }
    )
  }
  else {
    # Just combine the interaction matrices.
    PerCapitaDynamics <- DynamicsFunction(
      rprime,
      Matrix::bdiag(InteractionMatrices$Mats)
    )
  }
}

# Instantiate Dispersal Matrix: #############################################
# First branch shouldn't be an issue, but may need to be fixed later by
# either hard coding or extracting the appropriate dispersalDict entry.
if ( result$NumEnvironments > 1) {
  DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
    EnvironmentDistances = convertDispersalDictToDistMatrix(
      dispersalDictionary,
      nEnv =  result$NumEnvironments
    ),
    SpeciesSpeeds = Pool$Speed
  )
} else {
  DispersalMatrix <- Matrix::sparseMatrix(
    i = {}, j = {}, # From documentation
    dims = c(nrow(Pool), nrow(Pool))
  )
}

# result_0.1_0.1 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 0.1,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 0.1,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 0.1), Reduced Propagule Size (Arrival Density * 0.1)"
# )
# save(result_0.1_0.1,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_01_01.RData"))
#
#
# result_0.1_1 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 0.1,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 1,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 0.1), Reduced Propagule Size (Arrival Density * 1)"
# )
# save(result_0.1_1,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_01_1.RData"))
#
# result_0.1_10 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 0.1,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 10,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 0.1), Reduced Propagule Size (Arrival Density * 10)"
# )
# save(result_0.1_10,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_01_10.RData"))
#
# result_1_0.1 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 1,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 0.1,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 1), Reduced Propagule Size (Arrival Density * 0.1)"
# )
# save(result_1_0.1,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_1_01.RData"))
#
# result_1_1 <- result
#
# result_1_10 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 1,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 10,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 1), Reduced Propagule Size (Arrival Density * 10)"
# )
# save(result_1_10,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_1_10.RData"))
#
# result_10_0.1 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 10,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 0.1,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 10), Reduced Propagule Size (Arrival Density * 0.1)"
# )
# save(result_10_0.1,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_10_01.RData"))
#
# result_10_1 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 10,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 1,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 10), Reduced Propagule Size (Arrival Density * 1)"
# )
# save(result_10_1,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_10_1.RData"))
#
# result_10_10 <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
#   Pool = Pool,
#   PopulationInitial = result$Abundance[1, -1],
#   NumEnvironments = result$NumEnvironments,
#   CharacteristicRate = CharacteristicRate,
#   Events = list(Events = result$Events %>% dplyr::mutate(Success = NA)),
#   PerCapitaDynamics = PerCapitaDynamics,
#   DispersalMatrix = DispersalMatrix,
#   EliminationThreshold = result$Parameters$EliminationThreshold * 10,
#   ArrivalDensity = result$Parameters$ArrivalDensity * 10,
#   ExtinctionProportion = result$Parameters$ExtinctionProportion,
#   MaximumTimeStep = result$Parameters$MaximumTimeStep,
#   BetweenEventSteps = result$Parameters$BetweenEventSteps,
#   Verbose = FALSE,
#   # Using the ellipsis pass through feature:
#   Timescale = "Simulation",
#   ID = result$Ellipsis$ID,
#   Affinity = result$Ellipsis$Affinity,
#   Notes = "Reduced Elim. Thresh. (* 10), Reduced Propagule Size (Arrival Density * 10)"
# )
# save(result_10_10,
#      file = file.path(directory,
#                       "CompareEliminationThresholds",
#                       "ComparisonSim_10_10.RData"))

result_1_1 <- result
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_01_01.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_01_1.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_01_10.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_1_01.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_1_10.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_10_01.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_10_1.RData"))
load(file.path(directory,
               "CompareEliminationThresholds",
               "ComparisonSim_10_10.RData"))

# ####
# We want species that would persist under a different parameter set.
testEarlyElimination <- function(abundance) {
  traj <- apply(abundance[, -1], MARGIN = 1, PerCapitaDynamics, t = 0)
  for (row in (nrow(abundance) - 1):1) {
    abundance[row, -1] <-
      ifelse(abundance[row, -1] < 1e-3 &
               traj[[row]] < 0 &
               abundance[row+1, -1] == 0,
             0,
             abundance[row, -1])
  }
  return(abundance)
}

do.call(
  gridExtra::grid.arrange,
  lapply(list(result_0.1_0.1, result_0.1_1, result_0.1_10, # by row.
              result_1_0.1, result_1_1, result_1_10,
              result_10_0.1, result_10_1, result_10_10),
         FUN = function(r) {
           # Get comparators
           rows <- floor(seq(from = 1, to = nrow(r$Abundance), length.out = 10000))
           temp <- r$Abundance[rows, ]
           base <- result_1_1$Abundance[rows, -1]

           # Species about to be eliminated? Toss it out early to prevent noise,
           # since we're looking for different community configurations.
           # These are essentially false positives.
           # Work from largest time,
           #    if below 1e-3 & negative slope & future is 0, then eliminate.
           # (Ruling out stochastic extirpation and non-trivial presence.)
           temp <- testEarlyElimination(temp)
           base <- testEarlyElimination(cbind(temp[, 1], base))[, -1]

           # Species present in both places? We don't care then.
           elims <- temp[, -1] > 1e-3 & base > 1e-3
           temp[, -1] <- ifelse(elims, 0, temp[, -1])
           base <- ifelse(elims, 0, base)

           # Canberra distances, note that NaN <=> base + temp == 0 <=>
           # base == 0 & temp == 0, but we don't care about double absence.
           temp[, -1] <- (temp[, -1] - base) / (temp[, -1] + base)
           temp[, -1] <- ifelse(is.na(temp[,-1]), 0, temp[,-1])
           LawMorton1996_PlotAbundance(
             temp, guides = FALSE
           # ) + ggplot2::scale_y_log10(
           # ) + ggplot2::coord_cartesian(
           #   ylim = c(10^-5, 1)
           ) + ggplot2::ggtitle(paste(r$Parameters$EliminationThreshold,
                                      r$Parameters$ArrivalDensity))
         }
  )
)

load("CompareEliminationThresholds/Diversity_TSTS_Simulation_2-1-4-1-NA-1-1_30-30-53-1-1.RData")
Diversity_1_1 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_01_01.RData")
Diversity_0.1_0.1 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_01_1.RData")
Diversity_0.1_1 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_01_10.RData")
Diversity_0.1_10 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_1_01.RData")
Diversity_1_0.1 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_1_10.RData")
Diversity_1_10 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_10_01.RData")
Diversity_10_0.1 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_10_1.RData")
Diversity_10_1 <- Diversity
load("CompareEliminationThresholds/Diversity_ComparisonSim_10_10.RData")
Diversity_10_10 <- Diversity

do.call(
  gridExtra::grid.arrange,
  lapply(list(Diversity_0.1_0.1, Diversity_0.1_1, Diversity_0.1_10, # by row.
              Diversity_1_0.1, Diversity_1_1, Diversity_1_10,
              Diversity_10_0.1, Diversity_10_1, Diversity_10_10),
         FUN = function(d) {
           ggplot2::ggplot(d$Presence,
                           ggplot2::aes(x = Time, y = Species,
                                        color = interaction(Type, Affinity))
                           ) + ggplot2::geom_point(
                             shape = '.', show.legend = FALSE
                             )
         }
  )
)

do.call(
  gridExtra::grid.arrange,
  lapply(list(result_0.1_0.1, result_0.1_1, result_0.1_10, # by row.
              result_1_0.1, result_1_1, result_1_10,
              result_10_0.1, result_10_1, result_10_10),
         FUN = function(r) {
           # Get comparators
           rows <- floor(seq(from = 1, to = nrow(r$Abundance), length.out = 10000))
           temp <- r$Abundance[rows, ]
           base <- result_1_1$Abundance[rows, -1]

           # Species about to be eliminated? Toss it out early to prevent noise,
           # since we're looking for different community configurations.
           # These are essentially false positives.
           # Work from largest time,
           #    if below 1e-3 & negative slope & future is 0, then eliminate.
           # (Ruling out stochastic extirpation and non-trivial presence.)
           temp <- testEarlyElimination(temp)
           base <- testEarlyElimination(cbind(temp[, 1], base))[, -1]

           # Species present in both places? We don't care then.
           elims <- temp[, -1] > 1e-3 & base > 1e-3
           temp[, -1] <- ifelse(elims, 0, temp[, -1])
           base <- ifelse(elims, 0, base)

           temp[, -1] <- (temp[, -1] - base) / (temp[, -1] + base)
           temp[, -1] <- ifelse(is.na(temp[,-1]), 0, temp[,-1])

           r$Abundance <- temp
           presence <- RMTRCode2::Calculate_Species(r, bintimes = FALSE)

           r$Abundance[, -1] <- -temp[,-1]
           presence2 <- RMTRCode2::Calculate_Species(r, bintimes = FALSE)

           if(nrow(presence) > 0){
             presence2$Abundance <- -presence2$Abundance
             ggplot2::ggplot(dplyr::bind_rows(presence, presence2),
                             ggplot2::aes(x = Time, y = Species,
                                          color = Abundance)
             ) + ggplot2::geom_point(
               shape = '.'
             ) + ggplot2::scale_color_viridis_c(
               direction = -1,
               breaks = c(-1, -0.5, 0, 0.5, 1),
               limits = c(-1, 1)
             ) + ggplot2::ggtitle(paste(r$Parameters$EliminationThreshold,
                                        r$Parameters$ArrivalDensity))
           } else {
             ggplot2::ggplot()
           }
         }
  )
)


do.call(
  gridExtra::grid.arrange,
  lapply(list(Diversity_0.1_0.1, Diversity_0.1_1, Diversity_0.1_10, # by row.
              Diversity_1_0.1, Diversity_1_1, Diversity_1_10,
              Diversity_10_0.1, Diversity_10_1, Diversity_10_10),
         FUN = function(d) {
           ggplot2::ggplot(d$Diversity %>% dplyr::filter(Metric == "Alpha Hill:0", is.na(Subset)),
                           ggplot2::aes(x = Time, y = Value)
           ) + ggplot2::geom_point(
             shape = '.', show.legend = FALSE
           ) + ggplot2::scale_y_continuous(
             limits = c(0, 13)
           )
         }
  )
)

dplyr::bind_rows(
  Diversity_0.1_0.1$Diversity %>% dplyr::mutate(Vers = "0.1_0.1"),
  Diversity_0.1_1$Diversity %>% dplyr::mutate(Vers = "0.1_1"),
  Diversity_0.1_10$Diversity %>% dplyr::mutate(Vers = "0.1_10"),
  Diversity_1_0.1$Diversity %>% dplyr::mutate(Vers = "1_0.1"),
  Diversity_1_1$Diversity %>% dplyr::mutate(Vers = "1_1"),
  Diversity_1_10$Diversity %>% dplyr::mutate(Vers = "1_10"),
  Diversity_10_0.1$Diversity %>% dplyr::mutate(Vers = "10_0.1"),
  Diversity_10_1$Diversity %>% dplyr::mutate(Vers = "10_1"),
  Diversity_10_10$Diversity %>% dplyr::mutate(Vers = "10_10")
) %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% ggplot2::ggplot(
  ggplot2::aes(x = Time, y = Value, color = Vers)
) + ggplot2::geom_line(
) + ggplot2::coord_cartesian(
  ylim = c(0, 15)
) + ggplot2::facet_wrap(
  .~Subset
)

dplyr::bind_rows(
  Diversity_0.1_0.1$Diversity %>% dplyr::mutate(Vers = "0.1_0.1"),
  Diversity_0.1_1$Diversity %>% dplyr::mutate(Vers = "0.1_1"),
  Diversity_0.1_10$Diversity %>% dplyr::mutate(Vers = "0.1_10"),
  Diversity_1_0.1$Diversity %>% dplyr::mutate(Vers = "1_0.1"),
  Diversity_1_1$Diversity %>% dplyr::mutate(Vers = "1_1"),
  Diversity_1_10$Diversity %>% dplyr::mutate(Vers = "1_10"),
  Diversity_10_0.1$Diversity %>% dplyr::mutate(Vers = "10_0.1"),
  Diversity_10_1$Diversity %>% dplyr::mutate(Vers = "10_1"),
  Diversity_10_10$Diversity %>% dplyr::mutate(Vers = "10_10")
) %>% dplyr::filter(
  Metric == "Alpha Hill:0"
) %>% dplyr::mutate(
  Time = floor(Time / 10000) * 10000,
  Subset = ifelse(is.na(Subset), "NA", Subset)
) %>% ggplot2::ggplot(
  ggplot2::aes(x = Time, y = Value, color = Vers,
               group = interaction(Time, Vers, Subset))
) + ggplot2::geom_boxplot(
) + ggplot2::coord_cartesian(
  ylim = c(0, 15)
) + ggplot2::facet_wrap(
  .~Subset
)
