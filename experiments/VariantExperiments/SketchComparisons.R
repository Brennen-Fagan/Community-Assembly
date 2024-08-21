# Preparation: ################################################################

### Libraries: ################################################################
library(dplyr)
library(tidyr)
library(ggplot2)

### Files: ####################################################################
datfolder <-
  "TSTS_Simulations_18-1_9-9_2024-08-06"
results <-
  lapply(dir(datfolder, "Simulation", full.names = T),
         function(x) {nm <- load(x); return(get(nm))})
interventions <-
  lapply(dir(datfolder, "Intervention", full.names = T),
         function(x) {nm <- load(x); return(get(nm))})

### Adjust for connections: ###################################################
# interventions store the ParentRun, which we need to assign to each result.
results <- lapply(seq_along(results), function(i, r, n) {
  r[[i]]$RunName <- n[i]
  return(r[[i]])
}, r = results, n = dir(datfolder, "Simulation", full.names = T))
names(results) <- dir(datfolder, "Simulation", full.names = T)
names(interventions) <- dir(datfolder, "Intervention", full.names = T)

numberAttempts <- 10000 # 30,000 simulation scale time units (simts) TS~8 simts
# For non-"suitability", this is the number of equispaced chronological points
# used for determining how the relationship develops through time.
# For "suitability", this is the number of random samples we take to determine
# how well the control is at replicating the unperturbed system.

# Note:
stopifnot(
  all(results[[1]]$ReactionTime ==
        unlist(lapply(c(results, interventions), function(r) r$ReactionTime)))
)
# So it doesn't matter which timescale we use, but it does matter start and end
# times.
# We have unperturbed-before, unperturbed-after, and perturbed-after.
unperturbedRange <-
  c(500, min(unlist(lapply(interventions, function(r) r$Abundance[1, 1]))))

# Work in TimeSincePerturbation units.
interventions <- lapply(
  interventions,
  function(r) {
    r$Abundance[, 1] <-
      (r$Abundance[, 1] - r$Ellipsis$Affinity$TimeIntervention)
    return(r)
  }
)
perturbedRange <-
  c(0, min(unlist(lapply(
    interventions,
    function(r) {
      (r$Events$Times[length(r$Events$Times)]
       - r$Ellipsis$Affinity$TimeIntervention)
    }))))
# REMINDER: MAKE SURE TO CONSIDER TEMPORAL PROXIMITY TO PERTURBATION.

# Schema: #####################################################################
# 1 aaaaaaaaaaaaaaaaaaaaaaaaaa/bbbbbbbbbbbbbbbbbbbbbbb Perturbation
# 2 cccccccccccccccccccccccccc/ccccccccccccccccccccccc Spatial Control (Theory)
# 3 aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa Unperturbed
# 4 cccccccccccccccccccccccccccccccccccccccccccccccccc Spatial Control
# Comparing 1 after / to 3 after /  -> True Perturbation Effect on focal patch.
# Comparing 2 after / to 4 after /  -> True Perturbation Effect on non-focal.
# Comparing 1 before / to 3 after / -> Temporal Control Suitability.
# Comparing 2 after / to 3 after /  -> Spatial Control Suitability.
# Comparing 1 before / to 1 after / -> Temporal Control.
# Comparing 2 after / to 1 after /  -> Spatial Control.

# Tools: ######################################################################
extractAbund <- function(abundance, nEnv, epsilon = NULL) {
  nSpec <- (ncol(abundance) - 1)/nEnv
  stopifnot(nSpec == floor(nSpec))
  function(env, time) {
    timerow <- which.max(time < abundance[, 1]) - 1
    stopifnot(env >= 1, env <= nEnv, timerow >= 1, timerow <= nrow(abundance))
    retval <-
      abundance[timerow, (env - 1) *nSpec + 1:nSpec + 1] # + 1 for Time Column
    if (is.null(epsilon)) {
      return(retval)
    }
    retval[retval < epsilon] <- 0
    return(retval)
  }
}

extractJaccDissimilarity <- function(ab1, ab2) {
  ab1 <- ab1 > 0; ab2 <- ab2 > 0
  if (sum(ab1) == 0) {
    if (sum(ab2) == 0) return(NaN) # 0 / 0
    return(0) # 0 / nonzero
  } else if (sum(ab2) == 0) return(0) # 0 / nonzero
  vegan::vegdist(method = "jaccard", x = rbind(ab1, ab2))[1]
}

extractBrayDissimilarity <- function(ab1, ab2) {
  if (sum(ab1) == 0) {
    if (sum(ab2) == 0) return(NaN) # 0 / 0
    return(0) # 0 / nonzero
  } else if (sum(ab2) == 0) return(0) # 0 / nonzero
  vegan::vegdist(method = "bray", x = rbind(ab1, ab2))[1]
}

# E.g.
# unperturbed1 <- extractAbund(
#   results[[1]]$Abundance,
#   results[[1]]$NumEnvironments,
#   epsilon = results[[1]]$Parameters$EliminationThreshold)
# unperturbed3 <- extractAbund(
#   results[[3]]$Abundance,
#   results[[3]]$NumEnvironments,
#   epsilon = results[[3]]$Parameters$EliminationThreshold)

# extractJaccDissimilarity(unperturbed1(1, 10000), unperturbed3(1, 10000))

# Comparisons: ################################################################
### Final Setup: ##############################################################
# TODO CHANGE FOR OTHER CONFIGURATIONS.
interventionDetails <- data.frame(
  InterventionRun = names(interventions),
  InterventionTime = unlist(lapply(
    interventions, function(r) r$Ellipsis$Affinity$TimeIntervention
  )),
  EnvIntervention = unlist(lapply(
    interventions, function(r) r$Ellipsis$Affinity$PatchInterventions
  )),
  EnvControl = unlist(lapply(
    interventions,
    function(r) setdiff(1:r$NumEnvironments,
                        r$Ellipsis$Affinity$PatchInterventions)
  )),
  stringsAsFactors = FALSE
)

# From above, but 2 + 4 -> ControlEnv, 1 + 3 InterventionEnv,
#                 1 + 2 after is InterventionRun
# Comparing 1 after / to 3 after /  -> True Focal.
# Comparing 2 after / to 4 after /  -> True Nonfocal.
# Comparing 1 before / to 3 after / -> Temporal Suitability
# Comparing 2 after / to 3 after /  -> Spatial Suitability.
# Comparing 1 before / to 1 after / -> Temporal Control.
# Comparing 2 after / to 1 after /  -> Spatial Control.
comparisonTypes <- c("True Focal",
                     "True Nonfocal",
                     "Temporal Suitability",
                     "Spatial Suitability",
                     "Temporal Control",
                     "Spatial Control")

abundances <- lapply(c(results, interventions), function(r) {
  extractAbund(
    r$Abundance,
    nEnv = r$NumEnvironments,
    epsilon = r$Parameters$EliminationThreshold)
})

# Expectations (developed after having seen some (corrupted) data mind you!):
# Temporal variations are 2-D due to varying start times.
# Other variations are 1-D due to paired nature.
# The Controls should aim to replicate the True Focal.
# The temporal control should not matter what comparison time is used.
#   (I.e., moving along x-axis should have about the same beta diversity.)
#   Hence the action should be in the vertical movement.
# Stripes should correspond to community configurations.
# In 1-D a sudden spike corresponds to a change in either comparator.
# In the bottom right of the plots, the differences should be decreased.
# After looking at the data more, a vertical stripe in temporal indicates an
#   insufficient burn-in.
# If suitability is good, the corresponding control should match the True Focal
#   well.
# One "nuisance" (that's a big problem in principle I guess?) is that the
#   non-perturbed focal patch could change configuration, and the controls
#   cannot pick up on that I think.
#

### Configuration of Comparisons: #############################################
# Note: we pair times for the 2-D exploration in temporal control and
#       temporal suitability.
temporalpairs <- matrix(c(
  runif(n = numberAttempts,
        min = unperturbedRange[1],
        max = unperturbedRange[2]),# col1: before
  runif(n = numberAttempts,
        min = perturbedRange[1],
        max = perturbedRange[2])# col2: after
  ), ncol = 2)
spatialpairs <- matrix(c(
  runif(n = numberAttempts,
        min = perturbedRange[1],
        max = perturbedRange[2]),# col1: before
  runif(n = numberAttempts,
        min = perturbedRange[1],
        max = perturbedRange[2])# col2: after
), ncol = 2)

attempts <- expand.grid(
  Replicate = 1:numberAttempts,
  ParentRun = names(results),
  InterventionType = c("12-1-p-2_1-1", "12-1-p-3_1-1"),
  ComparisonType = comparisonTypes,
  stringsAsFactors = FALSE
) %>% dplyr::mutate(
  InterventionRun =
    paste0(gsub(pattern = "Simulation_",
                replacement = "Intervention_",
                x = tools::file_path_sans_ext(ParentRun)),
           "_", InterventionType, ".RData")
) %>% dplyr::left_join(
  interventionDetails, by = "InterventionRun"
) %>% dplyr::mutate(
  Env1 = dplyr::case_when(
    ComparisonType == comparisonTypes[1] ~ EnvIntervention,
    ComparisonType == comparisonTypes[2] ~ EnvControl,
    ComparisonType == comparisonTypes[3] ~ EnvIntervention,
    ComparisonType == comparisonTypes[4] ~ EnvControl,
    ComparisonType == comparisonTypes[5] ~ EnvIntervention,
    ComparisonType == comparisonTypes[6] ~ EnvControl,
    TRUE ~ NA_integer_),
  Env2 = dplyr::case_when(
    ComparisonType == comparisonTypes[1] ~ EnvIntervention,
    ComparisonType == comparisonTypes[2] ~ EnvControl,
    ComparisonType == comparisonTypes[3] ~ EnvIntervention,
    ComparisonType == comparisonTypes[4] ~ EnvIntervention,
    ComparisonType == comparisonTypes[5] ~ EnvIntervention,
    ComparisonType == comparisonTypes[6] ~ EnvIntervention,
    TRUE ~ NA_integer_),
  Env1File = dplyr::case_when(
    ComparisonType == comparisonTypes[1] ~ InterventionRun,
    ComparisonType == comparisonTypes[2] ~ InterventionRun,
    ComparisonType == comparisonTypes[3] ~ ParentRun,
    ComparisonType == comparisonTypes[4] ~ InterventionRun,
    ComparisonType == comparisonTypes[5] ~ ParentRun,
    ComparisonType == comparisonTypes[6] ~ InterventionRun,
    TRUE ~ NA_character_),
  Env2File = dplyr::case_when(
    ComparisonType == comparisonTypes[1] ~ ParentRun,
    ComparisonType == comparisonTypes[2] ~ ParentRun,
    ComparisonType == comparisonTypes[3] ~ ParentRun,
    ComparisonType == comparisonTypes[4] ~ ParentRun,
    ComparisonType == comparisonTypes[5] ~ InterventionRun,
    ComparisonType == comparisonTypes[6] ~ InterventionRun,
    TRUE ~ NA_character_),
  Env1Abundance = unlist(lapply(Env1File,
                                function(f) which(f == names(abundances)))),
  Env2Abundance = unlist(lapply(Env2File,
                                function(f) which(f == names(abundances)))),
  Env1Timespan = dplyr::case_when(
    ComparisonType == comparisonTypes[1] ~ "After",
    ComparisonType == comparisonTypes[2] ~ "After",
    ComparisonType == comparisonTypes[3] ~ "Before",
    ComparisonType == comparisonTypes[4] ~ "After",
    ComparisonType == comparisonTypes[5] ~ "Before",
    ComparisonType == comparisonTypes[6] ~ "After",
    TRUE ~ NA_character_),
  Env2Timespan = dplyr::case_when(
    ComparisonType == comparisonTypes[1] ~ "After",
    ComparisonType == comparisonTypes[2] ~ "After",
    ComparisonType == comparisonTypes[3] ~ "After",
    ComparisonType == comparisonTypes[4] ~ "After",
    ComparisonType == comparisonTypes[5] ~ "After",
    ComparisonType == comparisonTypes[6] ~ "After",
    TRUE ~ NA_character_),
  Env1Time = dplyr::case_when(
    ComparisonType %in% comparisonTypes[c(3, 5)] &
      Env1Timespan == "Before" ~ temporalpairs[Replicate, 1],
    ComparisonType %in% comparisonTypes[c(4, 6)] &
      Env1Timespan == "After" ~ spatialpairs[Replicate, 1],

    Env1Timespan == Env2Timespan & Env1Timespan == "Before" ~
      unperturbedRange[1] +
      diff(unperturbedRange)/(numberAttempts - 1) * (Replicate - 1),
    Env1Timespan == Env2Timespan & Env1Timespan == "After" ~
      perturbedRange[1] +
      diff(perturbedRange)/(numberAttempts - 1) * (Replicate - 1),

    Env1Timespan == "Before" ~ runif(n = length(Env1Timespan),
                                     min = unperturbedRange[1],
                                     max = unperturbedRange[2]),
    Env1Timespan == "After" ~ runif(n = length(Env1Timespan),
                                    min = perturbedRange[1],
                                    max = perturbedRange[2]),
    TRUE ~ NA_real_
  ),
  Env2Time = dplyr::case_when(
    ComparisonType %in% comparisonTypes[c(3, 5)] &
      Env2Timespan == "After" ~ temporalpairs[Replicate, 2],
    ComparisonType %in% comparisonTypes[c(4, 6)] &
      Env2Timespan == "After" ~ spatialpairs[Replicate, 2],

    Env1Timespan == Env2Timespan & Env2Timespan == "Before" ~
      unperturbedRange[1] +
      diff(unperturbedRange)/(numberAttempts - 1) * (Replicate - 1),
    Env1Timespan == Env2Timespan & Env2Timespan == "After" ~
      perturbedRange[1] +
      diff(perturbedRange)/(numberAttempts - 1) * (Replicate - 1),

    Env2Timespan == "Before" ~ runif(n = length(Env2Timespan),
                                     min = unperturbedRange[1],
                                     max = unperturbedRange[2]),
    Env2Timespan == "After" ~  runif(n = length(Env2Timespan),
                                     min = perturbedRange[1],
                                     max = perturbedRange[2]),
    TRUE ~ NA_real_
  ),
  Env1TimeEval = dplyr::case_when(
    Env1Timespan == "After" &
      Env1File == ParentRun ~ Env1Time + InterventionTime,
    TRUE ~ Env1Time
  ),
  Env2TimeEval = dplyr::case_when(
    Env2Timespan == "After" &
      Env2File == ParentRun ~ Env2Time + InterventionTime,
    TRUE ~ Env2Time
  )
  #,
  # Jaccard = NA,
  # Bray = NA
)

### Evaluate Comparisons: #####################################################
attempts <- attempts %>% dplyr::group_by(# Rowwise
  Replicate, ParentRun, InterventionType, ComparisonType
) %>% dplyr::group_modify(
  .f = function(.x, .y) {# non-group cols, group cols.
    ab1 <- abundances[[.x$Env1Abundance]](.x$Env1, .x$Env1TimeEval)
    ab2 <- abundances[[.x$Env2Abundance]](.x$Env2, .x$Env2TimeEval)
    jacc <- extractJaccDissimilarity(ab1, ab2)
    bray <- extractBrayDissimilarity(ab1, ab2)
    return(.x %>% dplyr::mutate(
      Jaccard = jacc, Bray = bray
    ))
  }
) %>% dplyr::ungroup()

### Plot Comparisons: #########################################################
ggplot2::ggplot(
  attempts,
  ggplot2::aes(x = Env1Time, y = Env2Time, color = Jaccard)
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c()

ggplot2::ggplot(
  attempts,
  ggplot2::aes(x = Env1Time, y = Env2Time, color = Bray)
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c()

### Temporal: Distance Decay?: ################################################

attemptsTime <- attempts %>% dplyr::filter(
  ComparisonType %in% comparisonTypes[c(3, 5)]
)

# Generally, the temporal metrics become more dissimilar if the temporal
# control is further away (temporally!) from the comparison time, but the story
# appears to be somewhat complicated by varying intervention strengths.
attemptsTime %>% dplyr::mutate(
  TimeDistance =
    abs(
      (Env1TimeEval + ifelse(Env1File == InterventionRun,
                             InterventionTime,
                             0)) -
        (Env2TimeEval + ifelse(Env2File == InterventionRun,
                               InterventionTime,
                               0))
    )
) %>% ggplot2::ggplot(
  ggplot2::aes(x = TimeDistance, y = Jaccard, color = Env1Time)
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::geom_smooth(color = "red"
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c()

attemptsTime %>% dplyr::mutate(
  TimeDistance =
    abs(
      (Env1TimeEval + ifelse(Env1File == InterventionRun,
                             InterventionTime,
                             0)) -
        (Env2TimeEval + ifelse(Env2File == InterventionRun,
                               InterventionTime,
                               0))
    )
) %>% ggplot2::ggplot(
  ggplot2::aes(x = TimeDistance, y = Bray, color = Env1Time)
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::geom_smooth(color = "red"
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c()


### Mismatch with True Focal: #################################################
# Quite a few things we want to compare around this.
#   Difference between beta calculated from control and from unperturbed.
#     Time1, Time2, Color = True Focal(Time2) - Control
#   Predictability from the suitability calculation.
#     Suitability, True Focal(Time2) - Control, Color = Time2
#   Systemic deviations between control and unperturbed.
#     Control, True Focal, Color = Time2
# Be careful with temporals to keep any pairs.

##### Spatial: ################################################################
attemptsTRUE <- attempts %>% dplyr::filter(
  ComparisonType == comparisonTypes[1]
)

# Pivot Suitability to be in the same rows as Control. This keeps pairs.
attemptsSpace <- attempts %>% dplyr::filter(
  ComparisonType %in% comparisonTypes[c(4, 6)]
) %>% tidyr::pivot_wider(
  id_cols = c(# Shared Observation Columns (differs from Time version!)
    Replicate, ParentRun, InterventionRun, InterventionType, InterventionTime,
    EnvIntervention, EnvControl
  ),
  names_from = ComparisonType,
  values_from = c(Env1:Bray)
)

attemptsSpaceControlTime2Col <- attemptsSpace$`Env2Time_Spatial Control`

# Remove Duplicated Columns. Note that this is aggressive!
attemptsSpace <- attemptsSpace[!duplicated(as.list(attemptsSpace))]

# Need to match Env1Time of True Focal to Env2Time of Control.
# Retrieve the new index.
attemptsSpaceControlTime2ColIndex <-
  which(unlist(lapply(as.list(attemptsSpace),
                      function(x) all(x == attemptsSpaceControlTime2Col))))

# attemptsSpace <- dplyr::right_join(
#   attemptsTRUE, attemptsSpace,
#   by = c(# Shared Observation Columns (differs from Time version!)
#     "Replicate", "ParentRun", "InterventionRun",
#     "InterventionType", "InterventionTime",
#     "EnvIntervention", "EnvControl",
#     "Env1Time" = colnames(attemptsSpace)[attemptsSpaceControlTime2ColIndex]
#   ),
#   suffix = c("True", "")
# )
# Instead, we need to match the ControlTime2s to the "nearest" Env1Times.
# Easiest way to match up would be stepfun.
# Note first: all(diff(unique(attemptsTRUE$Env1Time)) > 0)
timeBins <- stepfun(unique(attemptsTRUE$Env1Time)[-1], # y 1 longer, take off 0
                    unique(attemptsTRUE$Env1Time))
attemptsSpace$TimeBin <-
  timeBins(attemptsSpace[, attemptsSpaceControlTime2ColIndex][[1]])
# such that: all(attemptsSpace$TimeBin %in% attemptsTRUE$Env1Time).
# Hence:
attemptsSpace <- dplyr::right_join(
  attemptsTRUE,
  attemptsSpace,
  by = c(# Shared Observation Columns (differs from Time version!)
    #"Replicate",
    "ParentRun", "InterventionRun",
    "InterventionType", "InterventionTime",
    "EnvIntervention", "EnvControl",
    "Env1Time" = "TimeBin"
  ),
  suffix = c("True", "")
)

attemptsSpace <- attemptsSpace[!duplicated(as.list(attemptsSpace))]

#   Difference between beta calculated from control and from unperturbed.
ggplot2::ggplot(
  attemptsSpace, ggplot2::aes(
    x = Env2TimeEval, y = Env1Time, color = `Jaccard_Spatial Control` - Jaccard
  )
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ ParentRun + InterventionType
) + ggplot2::scale_color_gradient2(
) + ggplot2::labs(
  color = "J(Pert. Focal, Pert. Control)\n-J(Pert. Focal, Unpert. Focal)"
)
ggplot2::ggplot(
  attemptsSpace, ggplot2::aes(
    x = Env2TimeEval, y = Env1Time, color = `Bray_Spatial Control` - Bray
  )
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_gradient2(
) + ggplot2::labs(
  color = "B(Pert. Focal, Pert. Control)\n-B(Pert. Focal, Unpert. Focal)"
)

#   Predictability from the suitability calculation.
ggplot2::ggplot(
  attemptsSpace, ggplot2::aes(
    x = `Jaccard_Spatial Suitability`,
    y = `Jaccard_Spatial Control` - Jaccard,
    color = Env2TimeEval
  )
) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "J(Unpert. Focal, Pert. Control) Dissimilarity",
  y = "J(Pert. Focal, Pert. Control) - J(Pert. Focal, Unpert. Focal)"
)

ggplot2::ggplot(
  attemptsSpace, ggplot2::aes(
    x = `Bray_Spatial Suitability`,
    y = `Bray_Spatial Control` - Bray,
    color = Env2TimeEval
  )
) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "B(Unpert. Focal, Pert. Control) Dissimilarity",
  y = "B(Pert. Focal, Pert. Control) - B(Pert. Focal, Unpert. Focal)"
)

#   Systemic deviations between control and unperturbed.
ggplot2::ggplot(
  attemptsSpace, ggplot2::aes(
    x = `Jaccard_Spatial Control`,
    y =  Jaccard,
    color = Env2TimeEval
  )
) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "J(Pert. Focal, Pert. Control) Dissimilarity",
  y = "J(Pert. Focal, Unpert. Focal)"
)
ggplot2::ggplot(
  attemptsSpace, ggplot2::aes(
    x = `Bray_Spatial Control`,
    y =  Bray,
    color = Env2TimeEval
  )
) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "B(Pert. Focal, Pert. Control) Dissimilarity",
  y = "B(Pert. Focal, Unpert. Focal)",
  color = "Time Since\nPerturbation"
)

# Density of differences
ggplot2::ggplot(
  attemptsSpace,
  ggplot2::aes(x = `Jaccard_Spatial Control` - Jaccard)
  # ) + ggplot2::geom_freqpoly(
  #   bins = 10000,
  #   ggplot2::aes(group = round(Env2TimeEval/10)*10,
  #                color = round(Env2TimeEval/10)*10),
  #   alpha = 0.1
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::geom_density(
) + ggplot2::geom_rug(
  ggplot2::aes(color = Env2TimeEval),
  alpha = 0.01
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "J(Pert. Focal, Pert. Control) - J(Pert. Focal, Unpert. Focal)",
  color = "Time Since\nPerturbation"
) + ggplot2::facet_grid(
  . ~ basename(ParentRun) + InterventionType
)
# But notice!
ggplot2::ggplot(
  attemptsSpace %>% dplyr::filter(
    abs(Env1Time - `Env1Time_Spatial Control`) / Env1Time < 0.05
  ),
  ggplot2::aes(x = `Jaccard_Spatial Control` - Jaccard)
  # ) + ggplot2::geom_freqpoly(
  #   bins = 10000,
  #   ggplot2::aes(group = round(Env2TimeEval/10)*10,
  #                color = round(Env2TimeEval/10)*10),
  #   alpha = 0.1
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::geom_density(
) + ggplot2::geom_rug(
  ggplot2::aes(color = Env2TimeEval),
  alpha = 0.01
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "J(Pert. Focal, Pert. Control) - J(Pert. Focal, Unpert. Focal)",
  color = "Time Since\nPerturbation"
) + ggplot2::facet_grid(
  . ~ basename(ParentRun) + InterventionType
)

ggplot2::ggplot(
  attemptsSpace,
  ggplot2::aes(x = `Bray_Spatial Control` - Bray)
# ) + ggplot2::geom_freqpoly(
#   bins = 10000,
#   ggplot2::aes(group = round(Env2TimeEval/100)*100,
#                color = round(Env2TimeEval/100)*100),
#   alpha = 0.1
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::geom_density(
) + ggplot2::geom_rug(
  ggplot2::aes(color = Env2TimeEval),
  alpha = 0.01
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "B(Pert. Focal, Pert. Control) - B(Pert. Focal, Unpert. Focal)",
  color = "Time Since\nPerturbation"
) + ggplot2::facet_grid(
  . ~ basename(ParentRun) + InterventionType
)

##### Temporal: ###############################################################
# Reminder:
# Quite a few things we want to compare around this.
#   Difference between beta calculated from control and from unperturbed.
#     Time1, Time2, Color = True Focal(Time2) - Control
#   Predictability from the suitability calculation.
#     Suitability, True Focal(Time2) - Control, Color = Time2
#   Systemic deviations between control and unperturbed.
#     Control, True Focal, Color = Time2
# Be careful with temporals to keep any pairs.

# Pivot Suitability to be in the same rows as Control. This keeps pairs.
attemptsTime <- attempts %>% dplyr::filter(
  ComparisonType %in% comparisonTypes[c(3, 5)]
) %>% tidyr::pivot_wider(
  id_cols = c(# Shared Observation Columns (differs from Time version!)
    Replicate, ParentRun, InterventionRun, InterventionType, InterventionTime,
    EnvIntervention, EnvControl
  ),
  names_from = ComparisonType,
  values_from = c(Env1:Bray)
)

attemptsTimeControlTime2Col <- attemptsTime$`Env2Time_Temporal Control`

# Remove Duplicated Columns. Note that this is aggressive!
attemptsTime <- attemptsTime[!duplicated(as.list(attemptsTime))]

# Need to match Env1Time of True Focal to Env2Time of Control.
# Retrieve the new index.
attemptsTimeControlTime2ColIndex <-
  which(unlist(lapply(as.list(attemptsTime),
                      function(x) all(x == attemptsTimeControlTime2Col))))

# If all of the times lined up, we could do.
# attemptsTime <- dplyr::right_join(
#   attemptsTRUE, attemptsTime,
#   by = c(# Shared Observation Columns (differs from Time version!)
#     "Replicate", "ParentRun", "InterventionRun",
#     "InterventionType", "InterventionTime",
#     "EnvIntervention", "EnvControl",
#     "Env1Time" = colnames(attemptsTime)[attemptsTimeControlTime2ColIndex]
#   ),
#   suffix = c("True", "")
# )
# Instead, we need to match the ControlTime2s to the "nearest" Env1Times.
# Easiest way to match up would be stepfun.
# Note first: all(diff(unique(attemptsTRUE$Env1Time)) > 0)
timeBins <- stepfun(unique(attemptsTRUE$Env1Time)[-1], # y 1 longer, take off 0
                    unique(attemptsTRUE$Env1Time))
attemptsTime$TimeBin <-
  timeBins(attemptsTime[, attemptsTimeControlTime2ColIndex][[1]])
# such that: all(attemptsTime$TimeBin %in% attemptsTRUE$Env1Time).
# Hence:
attemptsTime <- dplyr::right_join(
  attemptsTRUE,
  attemptsTime,
  by = c(# Shared Observation Columns (differs from Time version!)
    #"Replicate",
    "ParentRun", "InterventionRun",
    "InterventionType", "InterventionTime",
    "EnvIntervention", "EnvControl",
    "Env1Time" = "TimeBin"
  ),
  suffix = c("True", "")
)

attemptsTime <- attemptsTime[!duplicated(as.list(attemptsTime))]

#   Difference between beta calculated from control and from unperturbed.
ggplot2::ggplot(
  attemptsTime, ggplot2::aes(
    x = Env2TimeEval, y = Env1Time, color = `Jaccard_Temporal Control` - Jaccard
  )
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ ParentRun + InterventionType
) + ggplot2::scale_color_gradient2(
) + ggplot2::labs(
  color = "J(Pert. Focal, Hist. Focal)\n-J(Pert. Focal, Unpert. Focal)"
)
ggplot2::ggplot(
  attemptsTime, ggplot2::aes(
    x = Env2TimeEval, y = Env1Time, color = `Bray_Temporal Control` - Bray
  )
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_gradient2(
) + ggplot2::labs(
  color = "B(Pert. Focal, Hist. Focal)\n-B(Pert. Focal, Unpert. Focal)"
)

#   Predictability from the suitability calculation.
ggplot2::ggplot(
  attemptsTime, ggplot2::aes(
    x = `Jaccard_Temporal Suitability`,
    y = `Jaccard_Temporal Control` - Jaccard,
    color = Env2TimeEval
  )
) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "J(Unpert. Focal, Hist. Focal) Dissimilarity",
  y = "J(Pert. Focal, Hist. Focal) - J(Pert. Focal, Unpert. Focal)"
)

ggplot2::ggplot(
  attemptsTime, ggplot2::aes(
    x = `Bray_Temporal Suitability`,
    y = `Bray_Temporal Control` - Bray,
    color = Env2TimeEval
  )
) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "B(Unpert. Focal, Hist. Focal) Dissimilarity",
  y = "B(Pert. Focal, Hist. Focal) - B(Pert. Focal, Unpert. Focal)"
)

#   Systemic deviations between control and unperturbed.
ggplot2::ggplot(
  attemptsTime, ggplot2::aes(
    x = `Jaccard_Temporal Control`,
    y =  Jaccard,
    color = Env2TimeEval
  )
) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "J(Pert. Focal, Hist. Focal) Dissimilarity",
  y = "J(Pert. Focal, Unpert. Focal)"
)
ggplot2::ggplot(
  attemptsTime, ggplot2::aes(
    x = `Bray_Temporal Control`,
    y =  Bray,
    color = Env2TimeEval
  )
) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed"
) + ggplot2::geom_point(alpha = 0.1
) + ggplot2::facet_grid(
  ComparisonType ~ basename(ParentRun) + InterventionType
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "B(Pert. Focal, Hist. Focal) Dissimilarity",
  y = "B(Pert. Focal, Unpert. Focal)",
  color = "Time Since\nPerturbation"
)

# Density of differences
ggplot2::ggplot(
  attemptsTime,
  ggplot2::aes(x = `Jaccard_Temporal Control` - Jaccard)
  # ) + ggplot2::geom_freqpoly(
  #   bins = 10000,
  #   ggplot2::aes(group = round(Env2TimeEval/10)*10,
  #                color = round(Env2TimeEval/10)*10),
  #   alpha = 0.1
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::geom_density(
) + ggplot2::geom_rug(
  ggplot2::aes(color = Env2TimeEval),
  alpha = 0.01
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "J(Pert. Focal, Hist. Focal) - J(Pert. Focal, Unpert. Focal)",
  color = "Time Since\nPerturbation"
) + ggplot2::facet_grid(
  . ~ basename(ParentRun) + InterventionType
)
ggplot2::ggplot(
  attemptsTime,
  ggplot2::aes(x = `Bray_Temporal Control` - Bray)
  # ) + ggplot2::geom_freqpoly(
  #   bins = 10000,
  #   ggplot2::aes(group = round(Env2TimeEval/100)*100,
  #                color = round(Env2TimeEval/100)*100),
  #   alpha = 0.1
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::geom_density(
) + ggplot2::geom_rug(
  ggplot2::aes(color = Env2TimeEval),
  alpha = 0.01
) + ggplot2::scale_color_viridis_c(
) + ggplot2::labs(
  x = "B(Pert. Focal, Hist. Focal) - B(Pert. Focal, Unpert. Focal)",
  color = "Time Since\nPerturbation"
) + ggplot2::facet_grid(
  . ~ basename(ParentRun) + InterventionType
)
