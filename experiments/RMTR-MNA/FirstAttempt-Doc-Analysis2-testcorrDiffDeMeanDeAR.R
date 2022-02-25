load("~/Random-Matrix-Theory/RMTRCode2/experiments/RMTR-MNA/Viking_SaveDiversity_2021-12-28_2022-01-15/MNA-MasterB-Cases-Diversity-2-1.RData")

dplyr::bind_rows(lapply(1:10, function(i) Diversity[[i]]$alpha %>% dplyr::mutate(
Simulation = i, Space = "Ring", Distance = Diversity$SpaceMod, Modifier = "NA", ModIntensity = toString(Diversity$NeutralMod), Measurement = "Richness", Value = Richness)
)) -> Diversity

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

tempDeDiff <-  lapply(temp, function(x) {
  x %>% dplyr::mutate_at(
    as.character(1:10), ~ c(NA, diff(.x))
  )
})

tempDeMean <- lapply(tempDeDiff, function(x) {
x %>% dplyr::mutate_at(
as.character(1:10), ~ .x - mean(.x, na.rm = TRUE)
)
})

tempDeAR1 <- lapply(tempDeMean, function(x) {
x %>% dplyr::mutate_at(
as.character(1:10), ~ residuals(arima(.x, order = c(2,0,0)))
)
})

tempcorrs <- lapply(tempDeAR1, function(x) testcorr::rcorr.test(x[-1, ] %>% dplyr::select(`1`:`10`)))

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
tempcorrsplottable,
ggplot2::aes(x = Correlation, fill = PValue, group = PValue)
) + ggplot2::geom_histogram(
binwidth = 0.05
) + ggplot2::geom_vline(
xintercept = 0, linetype = "dashed"
) + ggplot2::facet_grid(
Modifier + ModIntensity ~ Space + Distance
))
