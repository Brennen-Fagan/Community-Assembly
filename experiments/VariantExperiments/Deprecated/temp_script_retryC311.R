load("MNA-MasterC-Setup-3-1-1.RData")
# If trying affinity. If not, comment out down to the result run.
Affinities <- rep(0, nrow(Pool)); PatchAffinities <- rep(0.5, 10);

rhofunction <- function(
  base = 2, offset = 0, multiplier = 1, metric = "euclidean"
) {
  force(base);force(offset);force(multiplier)
  function(m, n) {
    base ^ (offset - multiplier * dist(
      matrix(c(m, n), byrow = TRUE, nrow = 2), method = metric)
    )
  }
}

rho.noop <- function(m, n) {1}
rho.10.1.2.euclidean <- rhofunction(10, 1, 2)

grid <- expand.grid(
  pool = 1:nrow(Pool), # Fastest Varying
  patch = 1:length(PatchAffinities) # Slower Varying.
)
rprime <-
  rep(Pool$ReproductionRate, NumberEnvironments) *
  mapply(
    grid$pool,
    grid$patch,
    FUN = function(i, j) {
      rho.10.1.2.euclidean(
        Affinities[i],
        PatchAffinities[j]
      )[1]
    }
  ) ^ sign(rep(Pool$ReproductionRate, NumberEnvironments))
stopifnot(rprime[1:100] == Pool$ReproductionRate)

PerCapitaDynamics <- RMTRCode2::PerCapitaDynamics_Type1(
  rprime, Matrix::bdiag(InteractionMatrices$Mats))

result <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
  Pool = Pool,
  PopulationInitial = popInitial,
  NumEnvironments = NumberEnvironments,
  CharacteristicRate = CharacteristicRate,
  Events = Events,
  PerCapitaDynamics = PerCapitaDynamics,
  DispersalMatrix = DispersalMatrix,
  EliminationThreshold = params$EliminationThreshold,
  ArrivalDensity = params$ArrivalDensity,
  ExtinctionProportion = ExtirpationProportion,
  MaximumTimeStep = params$MaximumTimeStep,
  BetweenEventSteps = params$BetweenEventSteps,
  Verbose = FALSE,
  # Using the ellipsis pass through feature:
  Timescale = "Simulation",
  Affinity = list(
    SpeciesAffinities = Affinities,
    PatchAffinities = PatchAffinities,
    EffectiveReproductionRate = rprime
  )
)
save(result, file = "Rerun-AffHigh0vHalf-3-1.RData")



# Plots for Thoughts: #########################################################
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# options(bitmapType = "cairo")

## Pools: #####################################################################
# # I don't think the relationships look different for the evensplit case:
# tempenv <- new.env()
# load("MNA-MasterC-Setup-3-1-1.RData", envir = tempenv)
# tempenv2 <- new.env()
# load(file.path("TSTS_Simulations_2-1_30-30_2024-10-11",
#                "TSTS_PoolPatchDynamics_2-1_30-30.RData"), envir = tempenv2)
# ggplot2::ggplot(
#   rbind(tempenv$Pool %>% dplyr::mutate(Affinity = NA, Pool = "Original"),
#         tempenv2$Pool %>% dplyr::mutate(Affinity = ID%%2, Pool = "New")),
#   ggplot2::aes(x = Size, y = ReproductionRate, color = Type)
# ) + ggplot2::geom_hex(
# ) + ggplot2::scale_fill_viridis_c(
# ) + ggplot2::facet_wrap(
#   .~Affinity+Pool
# ) + ggplot2::scale_x_log10(
# )
# # Reassuringly, ks.tests don't reject it
# ks.test(tempenv$Pool$Size, tempenv2$Pool$Size)
# ks.test(tempenv$Pool$ReproductionRate, tempenv2$Pool$ReproductionRate)

## Interaction Matrices: ######################################################
# Do the MNA-311 matrices look different than the TSTS_2-1_30-30_2024-10-11's?
# If there is randomness in the interaction types, Chi-square on the
# interaction types some form of histogram for interaction effects seem like an
# obvious bare-minimum check. For now, we don't have random interactions like
# that (because every consumer eats everything, and there's essentially a
# fixed number of basals and consumers each).

# convertIntMatToDoubleDF <- function(mat, pool = NULL) {
#   # Return dataframe with To, From, ToType, FromType, ToChange, FromChange.
#
#   # Put in single direction format.
#   rownames(mat) <- 1:nrow(mat)
#   colnames(mat) <- 1:ncol(mat)
#   sing <- reshape2::melt(mat)
#   colnames(sing) <- c("To", "From", "Change")
#   sing <- sing %>% dplyr::filter(Change != 0)
#
#   # Now convert to the desired double direction format.
#   # In the matrix, its the effect TO the row FROM the column.
#   # As an initial convention, we'll use row < column (upper triangular).
#   doub <- rbind(dplyr::full_join(
#     sing %>% dplyr::filter(To < From) %>% dplyr::rename(ToChange = Change),
#     sing %>% dplyr::filter(To > From) %>% dplyr::rename(
#       tidyr::all_of(c(To = "From", From = "To", FromChange = "Change"))
#     ),
#     by = c("From", "To")
#   ),
#   sing %>% dplyr::filter(To == From) %>% dplyr::rename(
#     ToChange = Change) %>% dplyr::mutate(FromChange = NA)
#   )
#
#   # And assign types.
#   if (!is.null(pool))
#     doub <- dplyr::left_join(
#       dplyr::left_join(
#         doub,
#         pool %>% dplyr::select(ID, Type), by = c("To" = "ID")
#       ) %>% dplyr::rename(ToType = Type),
#       pool %>% dplyr::select(ID, Type), by = c("From" = "ID")
#     ) %>% dplyr::rename(FromType = Type)
#
#   return(doub)
# }
#
# tempenv <- new.env()
# load("MNA-MasterC-Setup-3-1-1.RData", envir = tempenv)
# tempenv2 <- new.env()
# load(file.path("TSTS_Simulations_2-1_30-30_2024-10-11",
#                "TSTS_PoolPatchDynamics_2-1_30-30.RData"), envir = tempenv2)
#
# tempdat <- dplyr::bind_rows(
#   lapply(seq_along(tempenv$InteractionMatrices$Mats), function(i, mats) {
#     convertIntMatToDoubleDF(mats[[i]], pool = tempenv$Pool) %>% dplyr::mutate(
#       Run = "Original", Env = i
#     )
#   }, mats = tempenv$InteractionMatrices$Mats),
#   lapply(seq_along(tempenv2$InteractionMatrices$Mats), function(i, mats) {
#     convertIntMatToDoubleDF(mats[[i]], pool = tempenv$Pool) %>% dplyr::mutate(
#       Run = "New", Env = i
#     )
#   }, mats = tempenv$InteractionMatrices$Mats)
# )
#
# tempdat <- dplyr::bind_rows(
#   tempdat %>% dplyr::filter(ToChange <= 0 | FromChange >= 0),
#   tempdat %>% dplyr::filter(ToChange > 0, FromChange < 0) %>% dplyr::rename(
#     tidyr::all_of(c(To = "From", From = "To",
#                     ToChange = "FromChange", FromChange = "ToChange",
#                     ToType = "FromType", FromType = "ToType"))
#   )
# )
#
# 2-D Histogram
# ggplot2::ggplot(
#   tempdat %>% dplyr::mutate(
#     FromChange = ifelse(is.na(FromChange), 0, FromChange)
#   ),
#   ggplot2::aes(x = ToChange,
#                y = FromChange#,
#                #color = interaction(ToType, FromType)
#                )
# ) + ggplot2::geom_hex(
# ) + ggplot2::scale_fill_viridis_c(trans = "log10"
# ) + ggplot2::facet_wrap(
#   .~Run+Env + interaction(ToType, FromType)
# )
#
# ks.test((tempdat %>% dplyr::filter(Run == "Original"))$FromChange,
#         (tempdat %>% dplyr::filter(Run == "New"))$ToChange)
# ks.test((tempdat %>% dplyr::filter(Run == "Original"))$ToChange,
#         (tempdat %>% dplyr::filter(Run == "New"))$ToChange)

## Events: ####################################################################
# Trivial to check whether the reruns give the same results.
# load("MNA-MasterC-Cases-Prepared-3-1.RData")
# resultenv <- new.env()
# load("Rerun-NoChanges-3-1.RData", envir = resultenv)
# all.equal(resultenv$result$Events, results[[1]]$Events)

# # What about the distribution of events between runs. Do those look the same?
# tempenv <- new.env()
# load("MNA-MasterC-Cases-Prepared-3-1.RData", envir = tempenv)
# tempenv2 <- new.env()
# load(file.path("TSTS_Simulations_2-1_30-30_2024-10-11",
#                "TSTS_Simulation_2-1-4-1-NA-1-1_30-30-53-1-1.RData"),
#      envir = tempenv2)
# tempenv$results[[1]]$Events <- tempenv$results[[1]]$Events %>% dplyr::mutate(
#   Times = Times / tempenv$results[[1]]$ReactionTime
# )
# tempenv2$result$Events <- tempenv2$result$Events %>% dplyr::mutate(
#   Times = Times / tempenv2$result$ReactionTime
# )
#
# # We should expect a difference in number of species.
# # What's odder to me is the difference in times. A factor of 6?
# ggplot2::ggplot(
#   rbind(
#     tempenv$results[[1]]$Events %>% dplyr::filter(
#       Environment == 1
#     ) %>% dplyr::mutate(Run = "Original"),
#     tempenv2$result$Events %>% dplyr::mutate(Run = "New")),
#   ggplot2::aes(x = Times, y = Species, color = Type, alpha = Success)
# ) + ggplot2::geom_point(
# ) + ggplot2::facet_wrap(
#   .~Run
# )
#
# # Checking the Number of Events over the characteristic time period yields
# # essentially the same rates:
# tempenv$results[[1]]$Events %>% dplyr::filter(
#   Environment == 1
# ) %>% with(table(Type)/max(Times))
# tempenv2$result$Events %>% with(table(Type)/max(Times))
#
# # But the number of events I'm running over within a simulation is drastically
# # fewer. This might be the change I made to the coupon collecting problem,
# # but also suggests that we might feel comfortable about making the changeover
# # much earlier than halfway.
# dim(tempenv$results[[1]]$Events %>% dplyr::filter(
#   Environment == 1
# ))
# dim(tempenv2$result$Events)
#
# # Of course, we already know we're getting ridiculously different results.
# # The instability of the systems we're assembling likely explains the large
# # number of successes in the "new" compared to the old amongst basals.
# # This should be verifiable with testing:
# chisq.test(
#   # Convert to frequencies
#   rbind(with(tempenv2$result$Events, table(Type, Success))[1:4],
#         with(tempenv$results[[1]]$Events %>% dplyr::filter(Environment == 1), table(Type, Success))[1:4]
#   )
# )
# # Which rejects the hypothesis of independence/same distribution easily,
# # but doesn't do so for the types alone.
# # And, of course, the time gaps should be similar:
# ks.test(diff(tempenv2$result$Events$Times),
#         diff(tempenv$results[[1]]$Events %>%
#                dplyr::filter(Environment == 1) %>% dplyr::pull(Times)))
# Contrast the species-type event distributions:
# (NOTING THE DIFFERING NUMBERS OF SPECIES, SO TWO NON-INDEPENDENT BATCHES)
# chisq.test(
#   # Convert to frequencies
#   rbind(with(tempenv2$result$Events, table(Type, Species))[1:200],
#         with(tempenv$results[[1]]$Events %>% dplyr::filter(Environment == 1), table(Type, Species))[1:200]
#   )
# )
# chisq.test(
#   # Convert to frequencies
#   rbind(with(tempenv2$result$Events, table(Type, Species))[201:400],
#         with(tempenv$results[[1]]$Events %>% dplyr::filter(Environment == 1), table(Type, Species))[1:200]
#   )
# )

## Abundance Matrices: ########################################################
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# options(bitmapType = "cairo")
# load("MNA-MasterC-Cases-Prepared-3-1.RData")
# resultenv <- new.env()
# load("Rerun-NoChanges-3-1.RData", envir = resultenv)
# temp <- results[[1]]$Abundance
# temp[,-1][temp[,-1] < results[[1]]$Parameters$EliminationThreshold] <- 0
# temp2 <- resultenv$result$Abundance
# temp2[,-1][temp2[,-1] < resultenv$result$Parameters$EliminationThreshold] <- 0
# temp3 <- 2*(temp - temp2) / (temp + temp2)
# temp3[is.na(temp3)] <- 0 # <=> temp + temp2 = 0, but => temp - temp2 = 0.
# obj <- gridExtra::grid.arrange(
#   RMTRCode2::LawMorton1996_PlotAbundance(
#     cbind(temp2[, 1], abs(temp3)), guides = FALSE
#   ) + ggplot2::scale_y_log10(limits = c(1e-4, NA)),
#   RMTRCode2::LawMorton1996_PlotAbundance(
#     cbind(temp2[, 1], abs(temp - temp2)), guides = FALSE
#   ) + ggplot2::scale_y_log10(limits = c(1e-4, NA)),
#   ncol = 2
# )
# ggplot2::ggsave(obj, filename = "rerun-abundance-equivalence.png")
# quit(save = "no")

# Can we ask if there are sustained differences?
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# options(bitmapType = "cairo")
# load("MNA-MasterC-Cases-Prepared-3-1.RData")
# resultenv <- new.env()
# load("Rerun-NoChanges-3-1.RData", envir = resultenv)
# temp <- results[[1]]$Abundance
# temp[,-1][temp[,-1] < results[[1]]$Parameters$EliminationThreshold] <- 0
# temp2 <- resultenv$result$Abundance
# temp2[,-1][temp2[,-1] < resultenv$result$Parameters$EliminationThreshold] <- 0
# temp3 <- 2*(temp - temp2) / (temp + temp2)
# temp3[is.na(temp3)] <- 0 # <=> temp + temp2 = 0, but => temp - temp2 = 0.
# mapRelvsAbs <- dplyr::left_join(
# reshape2::melt(temp3) %>% dplyr::rename(Relative = value),
# reshape2::melt(temp - temp2) %>% dplyr::rename(Absolute = value),
# by = c("Var1", "Var2"))
# mapRelvsAbs <- mapRelvsAbs %>% dplyr::filter(Relative != 0 | Absolute != 0)
# # Already a decent plot, although some care needed for interpretation:
# ggplot2::ggplot(mapRelvsAbs, ggplot2::aes(x = Relative, y = Absolute, group = Var2)) + ggplot2::geom_path()
# ggplot2::ggplot(mapRelvsAbs, ggplot2::aes(
#  x = abs(Relative), y = abs(Absolute), group = Var2)
# ) + ggplot2::geom_path() + ggplot2::scale_y_log10() + ggplot2::scale_x_log10()
#
# mapRelvsAbs <- mapRelvsAbs %>% dplyr::group_by(Var2) %>% dplyr::mutate(
#   NotContiguous = (Var1 - lag(Var1)) - 1,
#   NotContiguous = ifelse(is.na(NotContiguous), 0, NotContiguous),
#   Grouping = cumsum(NotContiguous != 0))
# # Then we can inspect how many of them are sustained above some threshold
# mapRelvsAbs %>% dplyr::group_by(Var2, Grouping) %>% dplyr::summarise(
#   RelThresh = any(abs(Relative) > 0.1),
#   RelAboveThresh = sum(abs(Relative) > 0.1),
#   Eliminated = any(abs(Relative) == 2),
#   TimeSpan = max(Time) - min(Time),
#   .groups = "drop"
# ) %>% dplyr::filter(RelThresh)
# # Only one looks like it deviates significantly and isn't on track to
# # elimination, but it also doesn't remain deviated too much.
