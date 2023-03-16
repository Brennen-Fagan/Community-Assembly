GeneralisedLotkaVolterra <- function(
  t, y, parms
) {
  # t is the current time point in the integration,
  # y is the current estimate of the variables in the ODE system.
  # parms is a vector or list of parameters
  # The return value of func should be a list, whose first element is a vector
  # containing the derivatives of y with respect to time, and whose next
  # elements are global values that are required at each point in times.
  # The derivatives must be specified in the same order as... y.
  # -- DESolve::rk documentation
  #
  # Our parms contains interaction matrix a and reproduction rates r.

  with(as.list(parms), {
    list(y * (r + a %*% y))
  })
}