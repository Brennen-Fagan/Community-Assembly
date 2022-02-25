# Presence Absence Matrix
tempPresAbs <- data.frame(cbind(
  Time = Results$`MNA-FirstAttempt-Result-Env10-None.RData`$Abundance[, 1],
  (Results$`MNA-FirstAttempt-Result-Env10-None.RData`$Abundance[, -1]) > 0
  ))

# Binning for equivalent time steps.
tempPresAbs[, 1] <- floor(tempPresAbs[, 1] * 100) / 100

tempPresAbs <- tempPresAbs %>% dplyr::group_by(
  Time
) %>% dplyr::summarise_all(
  ~ .x[which.min(
    tail(dist(matrix(c(.x, median(.x)))), length(.x))
  )]
)

tempPresAbs <- tempPresAbs %>% tidyr::pivot_longer(
  cols = `X1`:`X1000`,
  names_to = "EnvSpecies",
  values_to = "Present"
) %>% dplyr::mutate(
  EnvSpecies = as.numeric(substring(EnvSpecies, 2)),
  Environment = ((EnvSpecies - 1) %/% 100 + 1),
  Species = ((EnvSpecies - 1) %% 100) + 1
) %>% dplyr::select(
  -EnvSpecies
)

tempPresAbs %>% dplyr::group_by(Species) %>% dplyr::group_map(
  ~ .x %>% tidyr::pivot_wider(
    names_from = "Environment", values_from = "Present"
  ) %>% dplyr::select(`1`:`10`) %>% testcorr::rcorr.test()
)

tempPresAbsADF <- tempPresAbs %>% dplyr::group_by(Species) %>% dplyr::group_map(
  ~ .x %>% tidyr::pivot_wider(
    names_from = "Environment",
    values_from = "Present"
  ) %>% dplyr::select(
    `1`:`10`
  ) %>% lapply(function(x) tseries::adf.test(x)$p.value)
) %>% unlist() %>% matrix(nrow = 10)

image(1:ncol(tempPresAbsADFs), 1:nrow(tempPresAbsADFs),
      t(log10(tempPresAbsADFs)),
      col = terrain.colors(60),
      axes = FALSE, main = "ADF Significance")

tempPresAbsKPSS <- tempPresAbs %>% dplyr::group_by(Species) %>% dplyr::group_map(
  ~ .x %>% tidyr::pivot_wider(
    names_from = "Environment",
    values_from = "Present"
  ) %>% dplyr::select(
    `1`:`10`
  ) %>% lapply(function(x) tseries::kpss.test(x)$p.value)
) %>% unlist() %>% matrix(nrow = 10)

image(1:ncol(tempPresAbsKPSS), 1:nrow(tempPresAbsKPSS),
      t(log10(tempPresAbsKPSS)),
      col = terrain.colors(60),
      axes = FALSE, main = "KPSS Significance")


# Autocorrelations
tempPresAbsacf <- tempPresAbs %>% dplyr::group_by(Species) %>% dplyr::group_map(
  ~ .x %>% tidyr::pivot_wider(
    names_from = "Environment",
    values_from = "Present"
  ) %>% dplyr::select(
    `1`:`10`
  ) %>% lapply(function(x) if(any(diff(x) > 0)) acf(x) else NA)
)

# Partial Autocorrelations
tempPresAbspacf <- tempPresAbs %>% dplyr::group_by(Species) %>% dplyr::group_map(
  ~ .x %>% tidyr::pivot_wider(
    names_from = "Environment",
    values_from = "Present"
  ) %>% dplyr::select(
    `1`:`10`
  ) %>% lapply(function(x) if(any(diff(x) > 0)) pacf(x) else NA)
)

tempPresAbsMean <- tempPresAbs %>% dplyr::group_by(Species) %>% dplyr::group_map(
  ~ .x %>% tidyr::pivot_wider(
    names_from = "Environment",
    values_from = "Present"
  ) %>% dplyr::select(
    `1`:`10`
  ) %>% lapply(changepoint::cpt.mean, method = "PELT")
)

tempPresAbsMeanMat <- unlist(lapply(tempPresAbsMean, function(lst) {
  lapply(lst, function(s4) {length(s4@param.est$mean)})
})) %>% matrix(nrow = 10)

image(1:ncol(tempPresAbsMeanMat), 1:nrow(tempPresAbsMeanMat),
      t(tempPresAbsMeanMat - 1),
      col = terrain.colors(60),
      axes = FALSE, main = "Changepoints")

tempPresAbsTimeSlices <- tempPresAbs %>% dplyr::group_by(Species) %>% dplyr::group_map(
  ~ do.call(cbind, lapply(4:7, dat = .x %>% tidyr::pivot_wider(
    names_from = "Environment",
    values_from = "Present"
  ), function(t, dat) {
    dat %>% dplyr::filter(
      t < Time, Time <= t + 1
    ) %>% dplyr::select(
      `1`:`10`
    ) %>% dplyr::rename_all(
      ~ paste0("TF", t, t + 1, " ", .)
    )
  }))
)

tempPresAbsTimeSlicescorrs <- lapply(tempPresAbsTimeSlices, testcorr::rcorr.test,
                                     plot = FALSE, table = FALSE)

tempPresAbsTimeSlicescorrsplottable <-
  dplyr::bind_rows(lapply(seq_along(tempPresAbsTimeSlicescorrs), function(i, tcs) {
    data.frame(
      Correlation = tcs[[i]]$pc[lower.tri(tcs[[i]]$pc)],
      PValue = tcs[[i]]$pv[lower.tri(tcs[[i]]$pv)],
      row = row(tcs[[i]]$pc)[lower.tri(tcs[[i]]$pc)],
      col = col(tcs[[i]]$pc)[lower.tri(tcs[[i]]$pc)],
      Species = i,
      Space = "None",
      Distance = "1e+00",
      Modifier = "Result",
      ModIntensity = "NA",
      stringsAsFactors = FALSE
    )
  }, tcs = tempPresAbsTimeSlicescorrs))

tempgrob1 <- ggplot2::ggplot(
  tempPresAbsTimeSlicescorrsplottable,
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "Species Pres.-Abs. Time Slice Corr."
)

tempgrob2 <- ggplot2::ggplot(
  tempPresAbsTimeSlicescorrsplottable %>% dplyr::filter(row == col + 10),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "Same Environment, Lag 10,000"
)

tempgrob3 <- ggplot2::ggplot(
  tempPresAbsTimeSlicescorrsplottable %>% dplyr::filter(row == col + 20),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "Same Environment, Lag 20,000"
)

tempgrob4 <- ggplot2::ggplot(
  tempPresAbsTimeSlicescorrsplottable %>% dplyr::filter(
    (1 <= row & row <= 10 & 1 <= col & col <= 10) |
      (11 <= row & row  <= 20 & 11 <= col & col <= 20) |
      (21 <= row & row  <= 30 & 21 <= col & col <= 30)
  ),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "Different Environments, No Lag"
)

tempgrob5 <- ggplot2::ggplot(
  tempPresAbsTimeSlicescorrsplottable %>% dplyr::filter(
    Species < 34
  ),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "No Consumers, All Lags"
)

tempgrob6 <- ggplot2::ggplot(
  tempPresAbsTimeSlicescorrsplottable %>% dplyr::filter(
    Species > 34
  ),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "No Basals, All Lags"
)

print(
  do.call(gridExtra::grid.arrange, list(tempgrob1, tempgrob2,
                                        tempgrob3, tempgrob4,
                                        tempgrob5, tempgrob6))
)
