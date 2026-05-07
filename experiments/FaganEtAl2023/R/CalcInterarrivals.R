CalcInterarrivals <- function(presenceChanges) {
  presenceChanges |> dplyr::group_by(
    Simulation
  ) |> dplyr::arrange(
    Time
  ) |> dplyr::mutate(
    LagTime = dplyr::lag(Time),
    WaitTime = Time - LagTime,
    LagType = dplyr::lag(Type)
  )
  
}