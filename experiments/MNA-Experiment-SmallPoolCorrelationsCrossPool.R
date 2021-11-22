# Post-MNA-Experiment-SmallPoolCorrelationsFragged.R
# We check for correlations across pools.

targets <- c(2, 4)

temp <- Events[[targets[1]]]
temp$Events$Type[temp$Events$Environment != 1] <- "Dummy"
temp2 <- Events[[targets[2]]]
temp2$Events$Type[temp2$Events$Environment != 2] <- "Dummy"
tempevents <- list(Events = rbind(temp$Events, temp2$Events),
                   Seed = unique(temp$Seed, temp2$Seed))
tempevents$Events <- tempevents$Events %>% dplyr::arrange(Times)

tempcharrate <- min(CharacteristicRates[[targets[1]]],
                    CharacteristicRates[[targets[2]]])

tempouts <- lapply(
  1, function(
    i, pools, mats, rates, events, dispmats,
    elm, arrdens, maxt, bett, targs = targets
  ) {
    print(Sys.time())

    # Divide up.
    tempresults <- lapply(
      1:length(targs),
      function(j, pl, Di, mat, rate, events, dynamics) {
        dispmat <- Di[[j]][1:nrow(pl[[j]]) + nrow(pl[[j]]) * (j - 1),
                           1:nrow(pl[[j]]) + nrow(pl[[j]]) * (j - 1)]
        dyns <- PerCapitaDynamics_Type1(pl[[j]]$ReproductionRate,
                                        mat[[j]]$Mats[[j]],
                                        NumEnvironments = 1)
        eventset <- events
        eventset$Events$Type[eventset$Events$Environment != j] <- "Dummy"
        eventset$Events$Environment <- 1

        MultipleNumericalAssembly_Dispersal(
          Pool = pl[[j]],
          NumEnvironments = 1,
          CharacteristicRate = rate,
          Events = eventset,
          PerCapitaDynamics = dyns,
          DispersalMatrix = dispmat,
          EliminationThreshold = elm,
          ArrivalDensity = arrdens,
          MaximumTimeStep = maxt,
          BetweenEventSteps = bett,
          Verbose = FALSE
        )
      },
      pl = pools,
      mat = mats,
      rate = rates,
      Di = dispmats,
      events = events
    )

    # Rebuild:
    result <- list()
    result$Events <- do.call(
      rbind,
      lapply(seq_along(tempresults), function(i, dat) {
        temp <- dat[[i]]$Events[!is.na(dat[[i]]$Events$Success), ]
        temp$Environment <- i
        temp
      }, dat = tempresults)
    )
    result$Events <- result$Events[order(result$Events$Times), ]

    result$Abundance <- cbind(
      time = tempresults[[1]]$Abundance[, 1],
      do.call(
        cbind,
        lapply(seq_along(tempresults), function(i, x, pl, targ = targs) {
          x[[i]]$Abundance[
            , 1:nrow(pl[[i]]) + 1
          ]
        }, x = tempresults, pl = pools)
      )
    )

    result$NumEnvironments <- length(tempresults)
    result$ReactionTime <- unique(unlist(lapply(
      tempresults, function(x) x$ReactionTime
    )))
    result$HistorySeed <- unique(unlist(lapply(
      tempresults, function(x) x$HistorySeed
    )))
    result$Parameters <- unique(unlist(lapply(
      tempresults, function(x) x$Parameters
    ), recursive = FALSE))
    names(result$Parameters) <- names(tempresults[[1]]$Parameters)
    result$Ellipsis <- unique(unlist(lapply(
      tempresults, function(x) x$Ellipsis
    ), recursive = FALSE))

    return(result)
  },
  pools = Pools[targets],
  rates = tempcharrate,
  mats = InteractionMatrices[targets],
  events = tempevents,
  dispmats = DispersalMatrices[targets],
  elm = EliminationThreshold,
  arrdens = ArrivalDensity,
  maxt = MaximumTimeStep,
  bett = BetweenEventSteps
)

tempdiversity <- lapply(tempouts, Calculate_Diversity)

temprich <- lapply(
  seq_along(tempdiversity),
  function(i, div, charrate) {
    div[[i]] %>% dplyr::filter(
      Measurement == "Richness",
      Environment != "Mean",
      Environment != "Gamma",
      Time > 100 / charrate
    ) %>% dplyr::mutate(
      Time = floor(Time / 10) * 10
    ) %>% tidyr::pivot_wider(
      names_from = Environment,
      values_from = Value,
      values_fn = median
      # Steps where a species is removed or added seem to
      # be resulting in two entries.
    )
  },
  div = tempdiversity,
  charrate = tempcharrate
)

tempcorr <- lapply(
  temprich,
  function(x) testcorr::rcorr.test(x %>% dplyr::select(`1`:`2`))
)
