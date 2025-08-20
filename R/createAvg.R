createAvg <- function(d, target, ps = ribbonQuantiles, cE = centralEst) {
  d |> dplyr::mutate(
    Time = floor(Time * time_averaging_size)/time_averaging_size
  ) |>  dplyr::group_by(
    Time, Pool, Noise, Neutral, Space
  ) |> dplyr::summarise(
    Low = unlist(dplyr::across(dplyr::any_of(target),
                               .fns = ~ quantile(.x, p = ps[1], na.rm = TRUE))),
    High = unlist(dplyr::across(dplyr::any_of(target),
                                .fns = ~ quantile(.x, p = ps[2], na.rm = TRUE))),
    Central = unlist(dplyr::across(dplyr::any_of(target),
                                   .fns = ~ cE(.x))),
    .groups = "drop"
  )
}