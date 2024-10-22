# Compare between old, new, and hybrid?
library(dplyr)
library(tidyr)
library(ggplot2)
options(bitmapType = "cairo")

# Loading: ####################################################################
oldenv <- new.env() # OLD -- Behaves well
newenv1 <- new.env() # NEW -- Behaves poorly
newenv2 <- new.env() # HYBRID -- Behaves...?
newenvEvents <- new.env()
newenvPool <- new.env()
newenvMats <- new.env()
newenvRerun <- new.env()
newenvPoolMat <- new.env()

load("MNA-MasterC-Setup-3-1-1.RData",
     envir = oldenv)
load("MNA-MasterC-Cases-Prepared-3-1.RData",
     envir = oldenv)

load(file.path("TSTS_Simulations_2-1_30-30_2024-10-11",
               "TSTS_PoolPatchDynamics_2-1_30-30.RData"),
     envir = newenv1)
load(file.path("TSTS_Simulations_2-1_30-30_2024-10-11",
               "TSTS_Simulation_2-1-4-1-NA-1-1_30-30-53-1-1.RData"),
     envir = newenv1)

load(file.path("TSTS_Simulations_7-1_34-34_2024-10-17",
               "TSTS_PoolPatchDynamics_7-1_34-34.RData"),
     envir = newenv2)
load(file.path("TSTS_Simulations_7-1_34-34_2024-10-17",
               "TSTS_Simulation_7-1-18-1-NA-7-15_34-34-57-81-13.RData"),
     envir = newenv2)

load(file.path("TSTS_Simulations_7-1_34-34-InjectEvents_2024-10-18",
               "TSTS_PoolPatchDynamics_7-1_34-34-InjectEvents.RData"),
     envir = newenvEvents)
load(file.path("TSTS_Simulations_7-1_34-34-InjectEvents_2024-10-18",
               "TSTS_Simulation_7-1-18-1-NA-7-15_34-34-57-81-13.RData"),
     envir = newenvEvents)

load(file.path("TSTS_Simulations_7-1_34-34-InjectPool_2024-10-18",
               "TSTS_PoolPatchDynamics_7-1_34-34-InjectPool.RData"),
     envir = newenvPool)
load(file.path("TSTS_Simulations_7-1_34-34-InjectPool_2024-10-18",
               "TSTS_Simulation_7-1-18-1-NA-7-15_34-34-57-81-13.RData"),
     envir = newenvPool)

load(file.path("TSTS_Simulations_7-1_34-34-SizeAndMatrices_2024-10-18",
               "TSTS_PoolPatchDynamics_7-1_34-34-SizeAndMatrices.RData"),
     envir = newenvMats)
load(file.path("TSTS_Simulations_7-1_34-34-SizeAndMatrices_2024-10-18",
               "TSTS_Simulation_7-1-18-1-NA-7-15_34-34-57-81-13.RData"),
     envir = newenvMats)

load("MNA-MasterC-Setup-3-1-1.RData",
     envir = newenvRerun)
load("Rerun-AffHigh0vHalf-3-1.RData",
     envir = newenvRerun)

load(file.path("TSTS_Simulations_7-1_34-34-InjectPoolMat_2024-10-21",
               "TSTS_PoolPatchDynamics_7-1_34-34-InjectPoolMat.RData"),
     envir = newenvPoolMat)
load(file.path("TSTS_Simulations_7-1_34-34-InjectPoolMat_2024-10-21",
               "TSTS_Simulation_7-1-18-1-NA-7-15_34-34-57-81-13.RData"),
     envir = newenvPoolMat)

# Compare the effective reproduction rates: equal to base rates? ##############
# Trivially true for old.
all(newenv1$Pool$ReproductionRate ==
      newenv1$result$Ellipsis$Affinity$EffectiveReproductionRate)
all(newenv2$Pool$ReproductionRate ==
      newenv2$result$Ellipsis$Affinity$EffectiveReproductionRate[1:100])
all(newenvEvents$Pool$ReproductionRate ==
      newenvEvents$result$Ellipsis$Affinity$EffectiveReproductionRate[1:100])
all(newenvPool$Pool$ReproductionRate ==
      newenvPool$result$Ellipsis$Affinity$EffectiveReproductionRate[1:100])
all(newenvMats$Pool$ReproductionRate ==
      newenvMats$result$Ellipsis$Affinity$EffectiveReproductionRate[1:100])
all(newenvRerun$Pool$ReproductionRate ==
      newenvRerun$result$Ellipsis$Affinity$EffectiveReproductionRate[1:100])
all(newenvPoolMat$Pool$ReproductionRate ==
      newenvPoolMat$result$Ellipsis$Affinity$EffectiveReproductionRate[1:100])
# So affinity is recorded as doing nothing!

# Compare the pools: are they obviously different, outside of affinity?
ggplot2::ggplot(
  rbind(oldenv$Pool %>% dplyr::mutate(Affinity = NA, Pool = "Original"),
        newenv1$Pool %>% dplyr::mutate(
          Affinity = newenv1$result$Ellipsis$Affinity$SpeciesAffinities,
          Pool = "New"
        ),
        newenv2$Pool %>% dplyr::mutate(
          Affinity = newenv2$result$Ellipsis$Affinity$SpeciesAffinities,
          Pool = "Hybrid"
        ),
        newenvEvents$Pool %>% dplyr::mutate(
          Affinity = newenvEvents$result$Ellipsis$Affinity$SpeciesAffinities,
          Pool = "InjectEvents"
        ),
        newenvPool$Pool %>% dplyr::mutate(
          Affinity = newenvPool$result$Ellipsis$Affinity$SpeciesAffinities,
          Pool = "InjectPool"
        ),
        newenvMats$Pool %>% dplyr::mutate(
          Affinity = newenvMats$result$Ellipsis$Affinity$SpeciesAffinities,
          Pool = "InjectSize&Mats"
        ),
        newenvRerun$Pool %>% dplyr::mutate(
          Affinity = newenvRerun$result$Ellipsis$Affinity$SpeciesAffinities,
          Pool = "RerunOriginal"
        ),
        newenvPoolMat$Pool %>% dplyr::mutate(
          Affinity = newenvPoolMat$result$Ellipsis$Affinity$SpeciesAffinities,
          Pool = "InjectPoolMat"
        )
  ),
  ggplot2::aes(x = Size, y = ReproductionRate, color = Type)
) + ggplot2::geom_hex(
) + ggplot2::scale_fill_viridis_c(
) + ggplot2::facet_wrap(
  .~Affinity+Pool
) + ggplot2::scale_x_log10(
)

# Reassuringly, ks.tests don't reject a shared population distribution.
ks.test(oldenv$Pool$Size, newenv1$Pool$Size)
ks.test(oldenv$Pool$Size, newenv2$Pool$Size)
ks.test(newenv1$Pool$Size, newenv2$Pool$Size)
ks.test(oldenv$Pool$ReproductionRate, newenv1$Pool$ReproductionRate)
ks.test(oldenv$Pool$ReproductionRate, newenv2$Pool$ReproductionRate)
ks.test(newenv1$Pool$ReproductionRate, newenv2$Pool$ReproductionRate)

# How about the Interaction Matrices? #########################################
# If there is randomness in the interaction types, Chi-square on the
# interaction types some form of histogram for interaction effects seem like an
# obvious bare-minimum check. For now, we don't have random interactions like
# that (because every consumer eats everything, and there's essentially a
# fixed number of basals and consumers each).

# Helper function to look at interactions instead of one-ways.
convertIntMatToDoubleDF <- function(mat, pool = NULL) {
  # Return dataframe with To, From, ToType, FromType, ToChange, FromChange.

  # Put in single direction format.
  rownames(mat) <- 1:nrow(mat)
  colnames(mat) <- 1:ncol(mat)
  sing <- reshape2::melt(mat)
  colnames(sing) <- c("To", "From", "Change")
  sing <- sing %>% dplyr::filter(Change != 0)

  # Now convert to the desired double direction format.
  # In the matrix, its the effect TO the row FROM the column.
  # As an initial convention, we'll use row < column (upper triangular).
  doub <- rbind(dplyr::full_join(
    sing %>% dplyr::filter(To < From) %>% dplyr::rename(ToChange = Change),
    sing %>% dplyr::filter(To > From) %>% dplyr::rename(
      tidyr::all_of(c(To = "From", From = "To", FromChange = "Change"))
    ),
    by = c("From", "To")
  ),
  sing %>% dplyr::filter(To == From) %>% dplyr::rename(
    ToChange = Change) %>% dplyr::mutate(FromChange = NA)
  )

  # And assign types.
  if (!is.null(pool))
    doub <- dplyr::left_join(
      dplyr::left_join(
        doub,
        pool %>% dplyr::select(ID, Type), by = c("To" = "ID")
      ) %>% dplyr::rename(ToType = Type),
      pool %>% dplyr::select(ID, Type), by = c("From" = "ID")
    ) %>% dplyr::rename(FromType = Type)

  return(doub)
}

tempdat <- dplyr::bind_rows(
  lapply(seq_along(oldenv$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = oldenv$Pool) %>% dplyr::mutate(
      Run = "Original", Env = i
    )
  }, mats = oldenv$InteractionMatrices$Mats),
  lapply(seq_along(newenv1$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = newenv1$Pool) %>% dplyr::mutate(
      Run = "New", Env = i
    )
  }, mats = newenv1$InteractionMatrices$Mats),
  lapply(seq_along(newenv2$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = newenv2$Pool) %>% dplyr::mutate(
      Run = "Hybrid", Env = i
    )
  }, mats = newenv2$InteractionMatrices$Mats),
  lapply(seq_along(newenvEvents$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = newenvEvents$Pool) %>% dplyr::mutate(
      Run = "InjectEvents", Env = i
    )
  }, mats = newenvEvents$InteractionMatrices$Mats),
  lapply(seq_along(newenvPool$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = newenvPool$Pool) %>% dplyr::mutate(
      Run = "InjectPool", Env = i
    )
  }, mats = newenvPool$InteractionMatrices$Mats),
  lapply(seq_along(newenvMats$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = newenvMats$Pool) %>% dplyr::mutate(
      Run = "InjectSize&Mats", Env = i
    )
  }, mats = newenvMats$InteractionMatrices$Mats),
  lapply(seq_along(newenvRerun$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = newenvRerun$Pool) %>% dplyr::mutate(
      Run = "RerunOriginal", Env = i
    )
  }, mats = newenvRerun$InteractionMatrices$Mats),
  lapply(seq_along(newenvPoolMat$InteractionMatrices$Mats), function(i, mats) {
    convertIntMatToDoubleDF(mats[[i]], pool = newenvPoolMat$Pool) %>% dplyr::mutate(
      Run = "InjectPoolMat", Env = i
    )
  }, mats = newenvPoolMat$InteractionMatrices$Mats)
)

# Rearrange so as to push predator prey into a single quadrant.
tempdat <- dplyr::bind_rows(
  tempdat %>% dplyr::filter(ToChange <= 0 | FromChange >= 0),
  tempdat %>% dplyr::filter(ToChange > 0, FromChange < 0) %>% dplyr::rename(
    tidyr::all_of(c(To = "From", From = "To",
                    ToChange = "FromChange", FromChange = "ToChange",
                    ToType = "FromType", FromType = "ToType"))
  )
)

# 2-D Histogram. Remember that New has more species and less environments!
# Attempting to average greatly reduces the perceived variance as well.
ggplot2::ggplot(
  tempdat %>% dplyr::mutate(
    FromChange = ifelse(is.na(FromChange), 0, FromChange)
    # ) %>% dplyr::group_by(To, From, ToType, FromType, Run) %>% dplyr::summarise(
    #   # Over Env:
    #   ToChange = mean(ToChange), FromChange = mean(FromChange)
  ) %>% dplyr::filter(
    Env == 1, From < 101
  ),
  ggplot2::aes(x = ToChange,
               y = FromChange#,
               #color = interaction(ToType, FromType)
  )
) + ggplot2::geom_hex(
) + ggplot2::scale_fill_viridis_c(trans = "log10"
) + ggplot2::facet_wrap(
  .~Run + interaction(ToType, FromType) # + Env
)

# Significance! (But not if we compare for envs for the same pool.)
ks.test(c((tempdat %>% dplyr::filter(Run == "Original"))$ToChange,
          (tempdat %>% dplyr::filter(Run == "Original"))$FromChange),
        c((tempdat %>% dplyr::filter(Run == "New"))$ToChange,
          (tempdat %>% dplyr::filter(Run == "New"))$FromChange))

ks.test(c((tempdat %>% dplyr::filter(Run == "Original"))$ToChange,
          (tempdat %>% dplyr::filter(Run == "Original"))$FromChange),
        c((tempdat %>% dplyr::filter(Run == "Hybrid"))$ToChange,
          (tempdat %>% dplyr::filter(Run == "Hybrid"))$FromChange))

ks.test(c((tempdat %>% dplyr::filter(Run == "New"))$ToChange,
          (tempdat %>% dplyr::filter(Run == "New"))$FromChange),
        c((tempdat %>% dplyr::filter(Run == "Hybrid"))$ToChange,
          (tempdat %>% dplyr::filter(Run == "Hybrid"))$FromChange))

# But the ecdfs are quite close, so this might be just too much data?
plot(ecdf(c((tempdat %>% dplyr::filter(Run == "Original"))$ToChange,
            (tempdat %>% dplyr::filter(Run == "Original"))$FromChange)))
plot(ecdf(c((tempdat %>% dplyr::filter(Run == "New"))$ToChange,
            (tempdat %>% dplyr::filter(Run == "New"))$FromChange)),
     add = TRUE, col = "red", lty = "dashed")
plot(ecdf(c((tempdat %>% dplyr::filter(Run == "Hybrid"))$ToChange,
            (tempdat %>% dplyr::filter(Run == "Hybrid"))$FromChange)),
     add = TRUE, col = "blue", lty = "dotted")

# This might be the avenue to follow up on?

# Events: #####################################################################
oldenv$results <- lapply(oldenv$results, function(r) {
  r$Events <- r$Events %>% dplyr::mutate(
    Times = Times / r$ReactionTime
  )
  return(r)
})
newenv1$result$Events <- newenv1$result$Events %>% dplyr::mutate(
  Times = Times / newenv1$result$ReactionTime
)
newenv2$result$Events <- newenv2$result$Events %>% dplyr::mutate(
  Times = Times / newenv2$result$ReactionTime
)
newenvEvents$result$Events <- newenvEvents$result$Events %>% dplyr::mutate(
  Times = Times / newenvEvents$result$ReactionTime
)
newenvPool$result$Events <- newenvPool$result$Events %>% dplyr::mutate(
  Times = Times / newenvPool$result$ReactionTime
)
newenvMats$result$Events <- newenvMats$result$Events %>% dplyr::mutate(
  Times = Times / newenvMats$result$ReactionTime
)
newenvRerun$result$Events <- newenvRerun$result$Events %>% dplyr::mutate(
  Times = Times / newenvRerun$result$ReactionTime
)
newenvPoolMat$result$Events <- newenvPoolMat$result$Events %>% dplyr::mutate(
  Times = Times / newenvPoolMat$result$ReactionTime
)

# We should expect a difference in number of species, and we've messed with
# the length of time of the simulation (but not with the event rates).
# The problem in the plot is that the hybrid and new have a visually distinct
# pattern of increased arrivals of basals compared to the old versions.
ggplot2::ggplot(
  dplyr::bind_rows(
    dplyr::bind_rows(lapply(oldenv$results, function(r) {
      r$Events %>% dplyr::filter(Environment == 1) %>% dplyr::mutate(
        Run = paste0("Old: ", r$HistorySeed), Species = Species / max(Species)
      )
    })),
    newenv1$result$Events %>% dplyr::mutate(
      Run = "New", Species = Species / max(Species)
      ),
    newenv2$result$Events %>% dplyr::filter(
      Environment == 1) %>% dplyr::mutate(
        Run = "Hybrid", Species = Species / max(Species)),
    newenvEvents$result$Events %>% dplyr::filter(
      Environment == 1) %>% dplyr::mutate(
        Run = "InjectEvents", Species = Species / max(Species)),
    newenvPool$result$Events %>% dplyr::filter(
      Environment == 1) %>% dplyr::mutate(
        Run = "InjectPool", Species = Species / max(Species)),
    newenvMats$result$Events %>% dplyr::filter(
      Environment == 1) %>% dplyr::mutate(
        Run = "InjectSize&Mats", Species = Species / max(Species)),
    newenvRerun$result$Events %>% dplyr::filter(
      Environment == 1) %>% dplyr::mutate(
        Run = "RerunOriginal", Species = Species / max(Species)),
    newenvPoolMat$result$Events %>% dplyr::filter(
      Environment == 1) %>% dplyr::mutate(
        Run = "InjectPoolMat", Species = Species / max(Species))
  ),
  ggplot2::aes(x = Times, y = Species, color = Type, alpha = Success)
) + ggplot2::geom_hline(
  yintercept = 1/3
) + ggplot2::geom_vline(
  xintercept = 4500
) + ggplot2::geom_point(
) + ggplot2::facet_wrap(
  .~Run
)

# Checking the Number of Events over the characteristic time period yields
# essentially the same rates over the patches:
oldenv$results[[1]]$Events %>% dplyr::filter(
  Environment == 1
) %>% with(table(Type)/max(Times))
newenv1$result$Events %>% with(table(Type)/max(Times))
newenv2$result$Events %>% dplyr::filter(
  Environment == 1
) %>% with(table(Type)/max(Times))

# Same distribution of events effectively.
chisq.test(
  # Convert to frequencies
  rbind(with(oldenv$results[[1]]$Events %>% dplyr::filter(Environment == 1), table(Type))[1:2],
        with(newenv1$result$Events %>% dplyr::filter(Environment == 1), table(Type))[1:2],
        with(newenv2$result$Events, table(Type))[1:2]
  )
)

# Even amongst species, accounting for the different numbers of course.
chisq.test(
  # Convert to frequencies
  rbind(with(oldenv$results[[1]]$Events %>% dplyr::filter(Environment == 1), table(Type, Species))[1:200],
        with(newenv1$result$Events %>% dplyr::filter(Environment == 1), table(Type, Species))[1:200],
        with(newenv2$result$Events, table(Type, Species))[1:200]
  ),
  simulate.p.value = TRUE, B = 1000
)

# Checking time gaps:
ks.test(
  diff(oldenv$results[[1]]$Events %>%
         dplyr::filter(Environment == 1) %>% dplyr::pull(Times)),
  diff(newenv1$result$Events$Times)
  )
ks.test(
  diff(oldenv$results[[1]]$Events %>%
         dplyr::filter(Environment == 1) %>% dplyr::pull(Times)),
  diff(newenv2$result$Events %>%
         dplyr::filter(Environment == 1) %>% dplyr::pull(Times))
)
ks.test(
  diff(newenv1$result$Events$Times),
  diff(newenv2$result$Events %>%
         dplyr::filter(Environment == 1) %>% dplyr::pull(Times))
)
