ArrivalFUN_Example2 <- function(Events, Rate) {
  if (Rate == 0 || Events == 0) {
    # Waiting finite time for infinite wait times.
    # Or waiting for no events to occur.
    # Return no events.
    return(numeric(0))
  }
  Duration <- Events / Rate # events / (events / s)
  # Note: Sum_i(Xi) ~ Gamma(shape = length(Xi), rate), Xi ~ Exp(rate)
  # mean = length(Xi) * rate
  # For Events "small" (possibly < 1), * 2 not large enough.
  # For Events "large", it appears to be far more than enough.
  # The + 1 is empirical; there are a few times where the system seems one event
  # away from meeting our quota.
  WaitingTimes <- cumsum(rexp(max(floor(Events * 2), 1) + 1, Rate))
  LastEvent <- which.max(WaitingTimes > Duration) - 1
  if (LastEvent == 0 && all(WaitingTimes < Duration)) {
    return(WaitingTimes)
  } else {
    WaitingTimes <- WaitingTimes[
      c(rep(TRUE, LastEvent),
        rep(FALSE, length(WaitingTimes) - LastEvent))
      ]
    return(WaitingTimes)
  }
}
ExtinctFUN_Example2 <- ArrivalFUN_Example2
