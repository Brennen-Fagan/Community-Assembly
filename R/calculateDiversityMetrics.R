# Based on RMTRCode2::Calculate_Diversity and calculateAbundanceMetrics as well
# as the betapart procedure.
# The main difference is to try to provide a uniform layout from the beginning
calculateDiversityMetrics <- function(
  abundance, nspecies, nenvironments, sizes
) {
  stopifnot(nenvironments >= 1)

  envs <- lapply(
    1:nenvironments,
    function(i, abund, numSpecies) {
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      env <- matrix(env, nrow = nrow(abund))
      return(env)
    },
    abund = abundance,
    numSpecies = sum(nspecies)
  )

  time <- abundance[, 1]

  # Format: data.frame: Time, Environment 1, Environment 2, Metric, Value
  # If an environment is not appropriate, we instead use NA.

  # Alpha: ####################################################################
  diversityAlpha <- dplyr::bind_rows(lapply(
    1:nenvironments,
    function(i) {
      vals <- hillWrapper(env = envs[[i]], i = i)
      abundance <-
        data.frame(
          Time = time,
          Environment1 = i,
          Environment2 = NA,
          Metric = "Abundance",
          Value = apply(envs[[i]], MARGIN = 1, function(x) sum(x))
        )
      # ~ individuals, can get average abundance by abundance/hill:0
      # similarly for biomass, can also get gamma by summing.
      if (!is.null(sizes)) {
        biomass <-
          data.frame(
            Time = time,
            Environment1 = i,
            Environment2 = NA,
            Metric = "Biomass",
            Value = apply(envs[[i]], MARGIN = 1, function(x) sum(x * sizes))
          )
        abundance <- rbind(abundance, biomass)
      }
      rbind(
        cbind(Time = rep(time, each = nrow(vals)/length(time)), vals),
        abundance
      )
    }
  )) %>% dplyr::mutate(
    Metric = paste0("Alpha ", Metric)
  )

  # Gamma: ####################################################################
  if (nenvironments > 1) {
    envgamma <- envs[[1]]
    if (nenvironments > 1) {
      for (i in 2:nenvironments) {
        envgamma <- envgamma + envs[[i]]
      }
    }

    diversityGamma <-
      hillWrapper(env = envgamma, i = NA) %>% dplyr::mutate(
        Metric = paste0("Gamma ", Metric)
      )
    diversityGamma <- cbind(
      Time = rep(time, each = nrow(diversityGamma)/length(time)),
      diversityGamma
    )
  }

  # Gamma Temporal: ###########################################################
  diversityGammaTemporal <- dplyr::bind_rows(lapply(
    1:nenvironments,
    function(i) {
      # Gamma (spatial) adds across environments; this adds across times.
      vals <- hillWrapper(env = colSums(envs[[i]]), i = i)
      cbind(Time = NA, vals)
    }
  )) %>% dplyr::mutate(
    Metric = paste0("GammaTemporal ", Metric)
  )

  # Beta Spatial: #############################################################
  if (nenvironments > 1) {
    diversityBetaSpace <- apply(
      abundance,
      MARGIN = 1, # Rows
      function(row, nenv) {
        thistime <- row[1]

        # List with three components:
        #   balance/turnover, gradient/nestedness, and total.
        distsBC <- betapart::beta.pair.abund(
          x = matrix(row[-1], nrow = nenv, byrow = TRUE),
          index.family = "bray"
        )
        distsJ <- betapart::beta.pair(
          x = matrix((row[-1] > 0) + 0, nrow = nenv, byrow = TRUE),
          index.family = "jaccard"
        )

        dataf <- expand.grid(
          Environment1 = 1:nenv,
          Environment2 = 1:nenv
        ) %>% dplyr::filter(
          Environment1 < Environment2
        ) %>% dplyr::arrange(
          Environment1, Environment2
        ) %>% dplyr::mutate(
          Time = thistime,
          SpaceBrayCurtisBalance = as.vector(distsBC$beta.bray.bal),
          SpaceBrayCurtisGradient = as.vector(distsBC$beta.bray.gra),
          SpaceBrayCurtis = as.vector(distsBC$beta.bray),
          SpaceJaccardTurnover = as.vector(distsJ$beta.jtu),
          SpaceJaccardNestedness = as.vector(distsJ$beta.jne),
          SpaceJaccard = as.vector(distsJ$beta.jac)
        )

        return(dataf)
      },
      nenv = nenvironments
    ) %>% dplyr::bind_rows() %>% tidyr::pivot_longer(
      cols = SpaceBrayCurtisBalance:SpaceJaccard,
      names_to = "Metric", values_to = "Value"
    )
  }

  # Beta Temporal: ############################################################
  # Picking most common 4 time differences, but assuming sorted already.
  # Note, we do it at this step to make sure we don't have a really common
  # primary difference causing more 2+-ary differences than other primaries.
  candidates <- sort(as.numeric(names( # Diff subtracts 1, hence 3.
    tail(sort(table(diff(time)), partial = length(time) - 3), n = 4)
  )))

  # Look for the candidate differences in the times.
  timediffs <- outer(time, time, `-`)
  # Arrange so that the candidate differences are what are analysed.
  timesets <- lapply(candidates, function(cand) {
    ro <- row(timediffs)[cand == timediffs]
    co <- col(timediffs)[cand == timediffs]
    data.frame(row = ro, col = co)
  })

  diversityBetaTime <- dplyr::bind_rows(lapply(
    1:nenvironments, function(i) {
      env <- envs[[i]]
      lapply(1:length(timesets), function(j) {
        rcs <- timesets[[j]]
        lapply(1:nrow(rcs), function(r) {
          target <- rbind(env[rcs$col[r], ], env[rcs$row[r], ])
          distsBC <- betapart::beta.pair.abund(
            x = target, index.family = "bray")
          distsJ <- betapart::beta.pair(
            x = (target > 0) + 0, index.family = "jaccard")
          retval <- expand.grid(
            Environment1 = i,
            Environment2 = NA
          ) %>% dplyr::mutate(
            Ti = time[rcs$col[r]], # Will fix in a moment, don't worry!
            TimeBrayCurtisBalance = as.vector(distsBC$beta.bray.bal),
            TimeBrayCurtisGradient = as.vector(distsBC$beta.bray.gra),
            TimeBrayCurtis = as.vector(distsBC$beta.bray),
            TimeJaccardTurnover = as.vector(distsJ$beta.jtu),
            TimeJaccardNestedness = as.vector(distsJ$beta.jne),
            TimeJaccard = as.vector(distsJ$beta.jac)
          )
          names(retval)[4:9] <- paste0(names(retval)[4:9], ": ", candidates[j])
          retval <- retval %>% tidyr::pivot_longer(
            cols = tidyr::starts_with("Time"),
            names_to = "Metric", values_to = "Value"
          ) %>% dplyr::rename(
            Time = Ti
          )
          return(retval)
        }) %>% dplyr::bind_rows()
      })  %>% dplyr::bind_rows()
    }
  ))

  return(dplyr::bind_rows(
    diversityAlpha,
    if (nenvironments > 1) diversityGamma,
    if (nenvironments > 1) diversityBetaSpace,
    diversityBetaTime
  ))
}
