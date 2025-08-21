# Like Calculating Diversities, but not for binary data only.
calculateAbundanceMetrics <- function(abundance, nspecies, nenvironments) {
  stopifnot(nenvironments >= 1)

  envs <- lapply(
    1:nenvironments,
    function(i, abund, numSpecies) {
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
    },
    abund = abundance,
    numSpecies = sum(nspecies)
  )

  time <- abundance[, 1]

  ### Alpha Diversity: ##################################################
  diversity_alpha <- lapply(
    1:nenvironments,
    function(i, abund, numSpecies) {
      env <- envs[[i]]
      env_basal <- env[, 1:nspecies[1], drop = FALSE]
      env_consumer <- env[, nspecies[1] + 1:nspecies[2], drop = FALSE]
      # fix a bug in vegan for single species ecosystems...
      if(ncol(env_consumer) == 1) env_consumer <- cbind(env_consumer, 0)
      metrics <- vegan::renyi(env, hill = TRUE)
      metrics_basal <- vegan::renyi(env_basal, hill = TRUE)
      metrics_consumer <- vegan::renyi(env_consumer, hill = TRUE)
      names(metrics) <- paste0(names(metrics), ", All")
      names(metrics_basal) <- paste0(names(metrics_basal), ", Basal")
      names(metrics_consumer) <- paste0(names(metrics_consumer), ", Consumer")
      cbind(Time = time,
            metrics, metrics_basal, metrics_consumer,
            Environment = i,
            stringsAsFactors = FALSE)
    },
    abund = abundance,
    numSpecies = sum(nspecies)
  )

  diversity_alpha <- dplyr::bind_rows(diversity_alpha)

  print("alpha")
  ### Gamma Diversity: ##################################################
  # In contrast, need to combine across environments.
  envgamma <- envs[[1]]
  if (nenvironments > 1) {
    for (i in 2:nenvironments) {
      envgamma <- envgamma + envs[[i]]
    }
  }

  env_basal <- envgamma[, 1:nspecies[1], drop = FALSE]
  env_consumer <- envgamma[, nspecies[1] + 1:nspecies[2], drop = FALSE]
  # fix a bug in vegan for single species ecosystems...
  if(ncol(env_consumer) == 1) env_consumer <- cbind(env_consumer, 0)

  metrics <- vegan::renyi(envgamma, hill = TRUE)
  metrics_basal <- vegan::renyi(env_basal, hill = TRUE)
  metrics_consumer <- vegan::renyi(env_consumer, hill = TRUE)
  names(metrics) <- paste0(names(metrics), ", All")
  names(metrics_basal) <- paste0(names(metrics_basal), ", Basal")
  names(metrics_consumer) <- paste0(names(metrics_consumer), ", Consumer")
  diversity_gamma <- cbind(Time = time,
                           metrics, metrics_basal, metrics_consumer,
                           Environment = NA,
                           stringsAsFactors = FALSE)

  print("gamma")
  ### Beta Diversity (Jaccard, Space): ##################################
  diversity_beta <- apply(
    abundance,
    MARGIN = 1, # Rows
    function(row, envs) {
      thistime <- row[1]

      # Vegan complains about rows with all 0's.
      # The warning is generic, so we cannot silence it specifically.
      dists <- suppressWarnings(vegan::vegdist(
        method = "bray",
        x = matrix(row[-1], nrow = envs, byrow = TRUE)
      ))

      dataf <- expand.grid(
        Env1 = 1:envs,
        Env2 = 1:envs
      ) %>% dplyr::filter(
        Env1 < Env2
      ) %>% dplyr::arrange(
        Env1, Env2
      ) %>% dplyr::mutate(
        Time = thistime,
        BrayCurtis = as.vector(dists)
      )

      return(dataf)
    },
    envs = nenvironments
  ) %>% dplyr::bind_rows()

  print("beta")
  ### Return Diversities: ###############################################
  return(list(
    alpha = diversity_alpha,
    beta = diversity_beta,
    gamma = diversity_gamma
  ))
}
