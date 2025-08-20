addNeutral <- function(loaded, presenceChanges, divide_time_by) {
  # Make Consistent
  events <- loaded$Events |> dplyr::filter(
    Times > burn_in
  ) |> dplyr::mutate(
    Times = Times / divide_time_by,
    Type = dplyr::case_when(
      Type == "Extinct" ~ "Extirpation",
      Type == "Arrival" ~ "Immigration",
      TRUE ~ "OOPS"
    ),
    Environment = as.character(Environment)
  ) |> dplyr::rename(
    Time = Times
  )
  
  events |> dplyr::full_join(
    presenceChanges, by = c("Time", "Species", "Environment", "Type")
  ) |> dplyr::mutate(
    Neutral = !is.na(Success),
    Success = ifelse(is.na(Success), TRUE, Success),
    Regional = ifelse(is.na(Regional), FALSE, Regional)
  )
}