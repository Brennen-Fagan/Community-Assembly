# source("standardizeGapSize.R") # maybe? Thinking I might try splinefun.
library(dplyr)

datfolder <-
  "TSTS_Simulations_18-1_9-9_2024-08-06"
results <-
  lapply(dir(datfolder, "Simulation", full.names = T),
         function(x) {nm <- load(x); return(get(nm))})
interventions <-
  lapply(dir(datfolder, "Intervention", full.names = T),
         function(x) {nm <- load(x); return(get(nm))})

# interventions store the ParentRun, which we need to assign to each result.
results <- lapply(seq_along(results), function(i, r, n) {
  r[[i]]$RunName <- n[i]
  return(r[[i]])
}, r = results, n = dir(datfolder, "Simulation", full.names = T))
names(results) <- dir(datfolder, "Simulation", full.names = T) 
names(interventions) <- dir(datfolder, "Intervention", full.names = T) 

numberAttempts <- 100#0
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

# Schema:
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

# Tool:
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

# TODO CHANGE FOR OTHER CONFIGURATIONS.
interventionDetails <- data.frame(
  InterventionRun = names(interventions),
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
comparisonTypes <- c("True Focal", "True Nonfocal", 
                     "Temporal Suitability", "Spatial Suitability",
                     "Temporal Control", "Spatial Control")

abundances <- lapply(c(results, interventions), function(r) {
  extractAbund(
    r$Abundance,
    nEnv = r$NumEnvironments,
    epsilon = r$Parameters$EliminationThreshold)
})

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
    Env1Timespan == Env2Timespan & Env2Timespan == "Before" ~ 
      unperturbedRange[1] + 
      diff(unperturbedRange)/(numberAttempts - 1) * (Replicate - 1),
    Env1Timespan == Env2Timespan & Env2Timespan == "After" ~ 
      perturbedRange[1] + 
      diff(perturbedRange)/(numberAttempts - 1) * (Replicate - 1),
    Env2Timespan == "Before" ~ runif(n = length(Env2Timespan), 
                                     min = unperturbedRange[1],
                                     max = unperturbedRange[2]),
    Env2Timespan == "After" ~ runif(n = length(Env2Timespan), 
                                    min = perturbedRange[1],
                                    max = perturbedRange[2]),
    TRUE ~ NA_real_
  )#,
  # Jaccard = NA,
  # Bray = NA
)

attempts <- attempts %>% dplyr::group_by(# Rowwise
  Replicate, ParentRun, InterventionType, ComparisonType
) %>% dplyr::group_modify(
  .f = function(.x, .y) {# non-group cols, group cols.
    jacc <- extractJaccDissimilarity(
      abundances[[.x$Env1Abundance]](.x$Env1, .x$Env1Time),
      abundances[[.x$Env2Abundance]](.x$Env2, .x$Env2Time)
      )
    bray <- extractBrayDissimilarity(
      abundances[[.x$Env1Abundance]](.x$Env1, .x$Env1Time),
      abundances[[.x$Env2Abundance]](.x$Env2, .x$Env2Time)
    )
    return(.x %>% dplyr::mutate(
      Jaccard = jacc, Bray = bray
    ))
  }
) %>% dplyr::ungroup()
