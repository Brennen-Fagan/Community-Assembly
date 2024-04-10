# Tests-calculateAbundanceMetrics.R

# Standard Opening: ###########################################################
library(RMTRCode2)
library(dplyr)

source("TimeSpaceAndTimeSeries-0-Functions.R")

print("Script stops early if a test fails.")

# Trivial Tests: ##############################################################
# Stochastic Tests: ###########################################################
# Deterministic Tests: ########################################################
environments <- 4; nspecies <- c(2, 1)
configurations <- list(
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    0:5, 5:0, rep(0, 6)
  )),
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    5:0, rep(0, 6), 0:5 
  )),
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    c(1, 1, 2, 2, 3, 3), c(1, 1, 2, 2, 3, 3), c(1, 1, 2, 2, 3, 3)
  )),
  matrix(byrow = FALSE, ncol = sum(nspecies), c(
    c(0, 0, 1, 2, 3, 3), rep(0, 6), rep(0, 6)
  ))
)
times <- matrix(1:nrow(configurations[[1]]), ncol = 1)
configurations[[1]] <- cbind(times, configurations[[1]])
ecosystem <- list(
  Abundance = do.call(cbind, configurations),
  Parameters = list(
    EliminationThreshold = 0.1
  ),
  NumEnvironments = environments
)

# Spatial, so compares between environments
result <- calculateAbundanceMetrics(
  abundance = ecosystem$Abundance, 
  nspecies = nspecies, 
  nenvironments = environments
)

# Reminder: only based on Species Identity, not abundance.
#           Addtionally breaks up by basal (1:2) and consumer (3).
#           No species: NaN.

expected <- list(
  alpha = expand.grid(
    Time = times,
    Environment = 1:environments
  ) %>% dplyr::group_by(Time, Environment) %>% dplyr::group_modify(
    .f = function(.x, .y) {
      env <- configurations[[.y$Environment]][.y$Time,, drop = FALSE]
      if(.y$Environment == 1) env <- env[, -1, drop = FALSE]
      env <- env/sum(env)
      envBasal <- env[, 1:nspecies[1], drop = FALSE]/
        sum(env[, 1:nspecies[1], drop = FALSE])
      envConsumer <- 
        env[, nspecies[1]+(1:nspecies[2]), drop = FALSE]/
        sum(env[, nspecies[1]+(1:nspecies[2]), drop = FALSE])
      
      retvalAll <- rbind.data.frame(
        c(lapply(c(0, 2^c((-2:6))), function(a) {
          if (sum(env != 0, na.rm = TRUE) == 0) {
            NA
          } else if (a != 1) {
            exp(1/(1 - a) * log(sum(env[env != 0]^a)))
          } else {
            exp(sum(- env[env != 0] * 
                      log(env[env != 0])))
          }
        }), 1 / max(env))
      )
      
      names(retvalAll) <- paste0(c(0, 2^c((-2:6)), Inf), ", All")
      
      retvalBasal <- rbind.data.frame(
        c(lapply(c(0, 2^c((-2:6))), function(a) {
          if (sum(envBasal != 0, na.rm = TRUE) == 0) {
            NA
          } else if (a != 1) {
            exp(1/(1 - a) * log(sum(envBasal[envBasal != 0]^a)))
          } else {
            exp(sum(- envBasal[envBasal != 0] * 
                      log(envBasal[envBasal != 0])))
          }
        }), 1 / max(envBasal)))
      
      names(retvalBasal) <- paste0(c(0, 2^c((-2:6)), Inf), ", Basal")
      
      retvalConsumer <- rbind.data.frame(
        c(lapply(c(0, 2^c((-2:6))), function(a) {
          if (sum(envConsumer != 0, na.rm = TRUE) == 0) {
            NA
          } else if (a != 1) {
            exp(1/(1 - a) * log(sum(envConsumer[envConsumer != 0]^a)))
          } else {
            exp(sum(- envConsumer[envConsumer != 0] * 
                      log(envConsumer[envConsumer != 0])))
          }
        }), 1 / max(envConsumer)))
      
      names(retvalConsumer) <- paste0(c(0, 2^c((-2:6)), Inf), ", Consumer")
      
      cbind(retvalAll, retvalBasal, retvalConsumer)
    }
  ) %>% dplyr::select(
    Time, `0, All`:`Inf, Consumer`, Environment
  ) %>% dplyr::arrange(Environment, Time) %>% as.data.frame,
  
  beta = expand.grid(
    Env1 = 1:environments,
    Env2 = 1:environments,
    Time = times
  ) %>% dplyr::filter(
    Env1 < Env2
  ) %>% dplyr::group_by(Env1, Env2, Time) %>% dplyr::group_modify(
    .f = function(.x, .y) {
      env1 <- configurations[[.y$Env1]][.y$Time,, drop = FALSE]
      if(.y$Env1 == 1) env1 <- env1[, -1, drop = FALSE]
      env2 <- configurations[[.y$Env2]][.y$Time,, drop = FALSE]
      if(.y$Env2 == 1) env2 <- env2[, -1, drop = FALSE]
      
      envs <- rbind(env1, env2)
      data.frame(
        BrayCurtis = 1 - 2 * sum(apply(envs, 2, min)) / sum(colSums(envs))
        )
    }
  ) %>% dplyr::arrange(Time, Env1, Env2) %>% as.data.frame,
  
  gamma = expand.grid(
    Time = times,
    Environment = NA
  ) %>% dplyr::group_by(Time, Environment) %>% dplyr::group_modify(
    .f = function(.x, .y) {
      env <- configurations[[1]][.y$Time, -1, drop = FALSE]
      for(config in configurations[-1]) {
        env <- env + config[.y$Time,, drop = FALSE]
      }
      
      env <- env/sum(env)
      envBasal <- env[, 1:nspecies[1], drop = FALSE]/
        sum(env[, 1:nspecies[1], drop = FALSE])
      envConsumer <- 
        env[, nspecies[1]+(1:nspecies[2]), drop = FALSE]/
        sum(env[, nspecies[1]+(1:nspecies[2]), drop = FALSE])
      
      retvalAll <- rbind.data.frame(
        c(lapply(c(0, 2^c((-2:6))), function(a) {
          if (sum(env != 0, na.rm = TRUE) == 0) {
            NA
          } else if (a != 1) {
            exp(1/(1 - a) * log(sum(env[env != 0]^a)))
          } else {
            exp(sum(- env[env != 0] * 
                      log(env[env != 0])))
          }
        }), 1 / max(env))
      )
      
      names(retvalAll) <- paste0(c(0, 2^c((-2:6)), Inf), ", All")
      
      retvalBasal <- rbind.data.frame(
        c(lapply(c(0, 2^c((-2:6))), function(a) {
          if (sum(envBasal != 0, na.rm = TRUE) == 0) {
            NA
          } else if (a != 1) {
            exp(1/(1 - a) * log(sum(envBasal[envBasal != 0]^a)))
          } else {
            exp(sum(- envBasal[envBasal != 0] * 
                      log(envBasal[envBasal != 0])))
          }
        }), 1 / max(envBasal)))
      
      names(retvalBasal) <- paste0(c(0, 2^c((-2:6)), Inf), ", Basal")
      
      retvalConsumer <- rbind.data.frame(
        c(lapply(c(0, 2^c((-2:6))), function(a) {
          if (sum(envConsumer != 0, na.rm = TRUE) == 0) {
            NA
          } else if (a != 1) {
            exp(1/(1 - a) * log(sum(envConsumer[envConsumer != 0]^a)))
          } else {
            exp(sum(- envConsumer[envConsumer != 0] * 
                      log(envConsumer[envConsumer != 0])))
          }
        }), 1 / max(envConsumer)))
      
      names(retvalConsumer) <- paste0(c(0, 2^c((-2:6)), Inf), ", Consumer")
      
      cbind(retvalAll, retvalBasal, retvalConsumer)
    }
  ) %>% dplyr::select(
    Time, `0, All`:`Inf, Consumer`, Environment
  ) %>% dplyr::arrange(Environment, Time) %>% as.data.frame
)

stopifnot(
  dim(expected$alpha) == dim(result$alpha),
  outer(
    1:nrow(expected$alpha), 
    1:ncol(expected$alpha), 
    FUN = function(i, j) {
      (is.na(expected$alpha) & is.na(result$alpha)) | 
        expected$alpha - result$alpha < 0.01
    }),
  dim(expected$beta) == dim(result$beta),
  outer(
    1:nrow(expected$beta), 
    1:ncol(expected$beta), 
    FUN = function(i, j) {
      (is.na(expected$beta) & is.na(result$beta)) | 
        expected$beta - result$beta < 0.01
    }),
  dim(expected$gamma) == dim(result$gamma),
  outer(
    1:nrow(expected$gamma), 
    1:ncol(expected$gamma), 
    FUN = function(i, j) {
      (is.na(expected$gamma) & is.na(result$gamma)) | 
        expected$gamma - result$gamma < 0.01
    })
)

# Standard Closing: ###########################################################
print("Success.")
