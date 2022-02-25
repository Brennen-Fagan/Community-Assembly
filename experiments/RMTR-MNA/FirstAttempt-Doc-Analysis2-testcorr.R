Diversity %>% dplyr::group_by(
  Simulation
) %>% dplyr::arrange(
  Time
) %>% dplyr::filter(
  Measurement == "Richness",
  Environment != "Mean",
  Environment != "Gamma",
  Time > 3, Time <= 6
) %>% dplyr::mutate(
  TimeFloor = floor(Time * 100) / 100,
  Environment = as.numeric(Environment)
  # Will be converted back to character, but this makes the order natural.
) %>% dplyr::group_by(
  Simulation,
  TimeFloor, Environment,
  Modifier, ModIntensity, Space, Distance
) %>% dplyr::summarise(
  Value = Value[which.min(
    tail(dist(matrix(c(Value, median(Value)))), length(Value))
  )]
) %>% tidyr::pivot_wider(
  names_from = Environment,
  values_from = Value,
  values_fn = list(ceiling)
  # Steps where a species is removed or added seem to
  # be resulting in two entries.
) %>% dplyr::ungroup() %>% dplyr::group_by(
  Simulation
) %>% dplyr::group_map(
  ~ .x
) -> temp

names(temp) <- unlist(lapply(temp,
                             function(x) toString(unique(c(x$Space, x$Distance, x$Modifier, x$ModIntensity)))))

tempDeMean <- lapply(temp, function(x) {
  x %>% dplyr::mutate_at(
    as.character(1:10), ~ .x - mean(.x, na.rm = TRUE)
  )
})

tempDeAR1 <- lapply(tempDeMean, function(x) {
  x %>% dplyr::mutate_at(
    as.character(1:10), ~ residuals(arima(.x, order = c(2,0,0)))
  )
})

tempcorrs <- lapply(temp, function(x) testcorr::rcorr.test(x %>% dplyr::select(`1`:`10`)))

#lapply(tempcorrs, function(tc) {hist(tc$pc[lower.tri(tc$pc)])})

tempcorrsplottable <- dplyr::bind_rows(lapply(seq_along(tempcorrs), function(i, tcs, nms) {
  nvals <- strsplit(nms[[i]], split = ", ", fixed = TRUE)[[1]]
  data.frame(
    Correlation = tcs[[i]]$pc[lower.tri(tcs[[i]]$pc)],
    PValue = tcs[[i]]$pv[lower.tri(tcs[[i]]$pv)],
    row = row(tcs[[i]]$pc)[lower.tri(tcs[[i]]$pc)],
    col = col(tcs[[i]]$pc)[lower.tri(tcs[[i]]$pc)],
    Space = nvals[1],
    Distance = nvals[2],
    Modifier = nvals[3],
    ModIntensity = nvals[4],
    stringsAsFactors = FALSE
  )
}, tcs = tempcorrs, nms = names(tempcorrs)))

print(ggplot2::ggplot(
  tempcorrsplottable %>% dplyr::filter(Space != "Ring"),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::facet_grid(
  Modifier + ModIntensity ~ Space + Distance
))

# Test: Alternative is Stationary
tempadfs <- lapply(temp, function(df) {
  lapply(df %>% dplyr::select(`1`:`10`), tseries::adf.test)
})

# Test: Null is Stationary
tempkpss <- lapply(temp, function(df) {
  lapply(df %>% dplyr::select(`1`:`10`), tseries::kpss.test)
})

# Autocorrelations
tempacf <- lapply(temp, function(df) {
  lapply(df %>% dplyr::select(`1`:`10`), acf)
})

# Partial Autocorrelations
temppacf <- lapply(temp, function(df) {
  lapply(df %>% dplyr::select(`1`:`10`), pacf)
})

# Compare [7,8] with [4,5] within and between environments
# in the no relationship case.
tempTimeSlices <- cbind(
  temp$`None, 1e+00, Result, NA` %>% dplyr::filter(
    TimeFloor >= 4, TimeFloor < 5
  ) %>% dplyr::select(`1`:`10`) %>% dplyr::rename_all(
    ~ paste("TF45", .)),
  temp$`None, 1e+00, Result, NA` %>% dplyr::filter(
    TimeFloor >= 5, TimeFloor < 6
  ) %>% dplyr::select(`1`:`10`) %>% dplyr::rename_all(
    ~ paste("TF56", .)),
  temp$`None, 1e+00, Result, NA` %>% dplyr::filter(
    TimeFloor >= 6, TimeFloor < 7
  ) %>% dplyr::select(`1`:`10`) %>% dplyr::rename_all(
    ~ paste("TF67", .)),
  temp$`None, 1e+00, Result, NA` %>% dplyr::filter(
    TimeFloor >= 7, TimeFloor < 8
  ) %>% dplyr::select(`1`:`10`) %>% dplyr::rename_all(
    ~ paste("TF78", .))
)

tempTimeSlicescorrs <- testcorr::rcorr.test(tempTimeSlices)

tempTimeSlicescorrsplottable <-
  data.frame(
    Correlation = tempTimeSlicescorrs$pc[lower.tri(tempTimeSlicescorrs$pc)],
    PValue = tempTimeSlicescorrs$pv[lower.tri(tempTimeSlicescorrs$pv)],
    row = row(tempTimeSlicescorrs$pc)[lower.tri(tempTimeSlicescorrs$pc)],
    col = col(tempTimeSlicescorrs$pc)[lower.tri(tempTimeSlicescorrs$pc)],
    Space = "None",
    Distance = "1e+00",
    Modifier = "Result",
    ModIntensity = "NA",
    stringsAsFactors = FALSE
  )

tempgrob1 <- ggplot2::ggplot(
  tempTimeSlicescorrsplottable,
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "Independent Islands' Time Slice Correlations"
)

tempgrob2 <- ggplot2::ggplot(
  tempTimeSlicescorrsplottable %>% dplyr::filter(row == col + 10),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "Same Environment, Lag 10,000"
)

tempgrob3 <- ggplot2::ggplot(
  tempTimeSlicescorrsplottable %>% dplyr::filter(row == col + 20),
  ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
  binwidth = 0.05
) + ggplot2::geom_vline(
  xintercept = 0, linetype = "dashed"
) + ggplot2::ggtitle(
  "Same Environment, Lag 20,000"
)

tempgrob4 <- ggplot2::ggplot(
  tempTimeSlicescorrsplottable %>% dplyr::filter(
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

print(
  do.call(gridExtra::grid.arrange, list(tempgrob1, tempgrob2, tempgrob3, tempgrob4))
)
