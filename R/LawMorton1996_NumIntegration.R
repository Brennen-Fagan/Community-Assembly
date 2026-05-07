LawMorton1996_NumIntegration <- function(
  A, R, X,
  OuterTimeStepSize = 1000,
  InnerTimeStepSize = 1,
  Tolerance = 0
) {
  # A = matrix of aij's (community matrix)
  # R = vector of ri's (basic reproduction rate)
  # X = vector of initial abundances.
  # OuterTimeStepSize is the length of the solution,
  # InnerTimeStepSize is the time in between records of the solution.
  
  ts <- seq(from = 0,
            to = OuterTimeStepSize,
            by = InnerTimeStepSize)
  
  deSolve::ode(
    X,
    times = ts,
    func = GeneralisedLotkaVolterra,
    parms = list(a = A, r = R, epsilon = Tolerance),
    events = list(func = function(t, y, parms) {
      y[y < parms$epsilon] <- 0
      y
    }, time = ts)
  )
}
