# Test case:
# Purely Mutualistic TF2010
temptemp <- PerCapitaDynamics_Mutualistic2(
  c(1, 2, -1, -2),
  matrix(c(
     -1,   0, 0.5,  0, rep(0, 4),
      0,  -1, 0.5,  0, rep(0, 4),
    0.5, 0.5,  -1,  0, rep(0, 4),
      0,   0,   0, -1, rep(0, 4),
    rep(0, 4),   -1,    0, 0.25,  0,
    rep(0, 4),    0,   -1, 0.75,  0,
    rep(0, 4), 0.25, 0.75,   -1,  0,
    rep(0, 4),    0,    0,    0, -1
  ), byrow = TRUE, nrow = 4 * 2, ncol = 4 * 2),
  2, c(2, 2),
  matrix(c(
    Inf,   0, 0.1,   0, rep(0, 4),
      0, Inf, 0.1,   0, rep(0, 4),
    0.2, 0.2, Inf,   0, rep(0, 4),
      0,   0,   0, Inf, rep(0, 4),
    rep(0, 4),  Inf,   0,  0.3,   0,
    rep(0, 4),    0, Inf,  0.3,   0,
    rep(0, 4),  0.4, 0.4,  Inf,   0,
    rep(0, 4),    0,   0,    0, Inf
  ), byrow = TRUE, nrow = 4 * 2, ncol = 4 * 2)
)


temptemp2 <- PerCapitaDynamics_Mutualistic2(
  c(1, 2, -1, -2),
  matrix(c(
     -1,   0, -0.5,  0, rep(0, 4),
      0,  -1, -0.5,  0, rep(0, 4),
    0.5, 0.5,   -1,  0, rep(0, 4),
      0,   0,    0, -1, rep(0, 4),
    rep(0, 4),   -1,    0, -0.25,  0,
    rep(0, 4),    0,   -1, -0.75,  0,
    rep(0, 4), 0.25, 0.75,    -1,  0,
    rep(0, 4),    0,    0,     0, -1
  ), byrow = TRUE, nrow = 4 * 2, ncol = 4 * 2),
  2, c(2, 2),
  matrix(c(
    Inf,   0, 0.1,   0, rep(0, 4),
    0, Inf, 0.1,   0, rep(0, 4),
    0.2, 0.2, Inf,   0, rep(0, 4),
    0,   0,   0, Inf, rep(0, 4),
    rep(0, 4),  Inf,   0,  0.3,   0,
    rep(0, 4),    0, Inf,  0.3,   0,
    rep(0, 4),  0.4, 0.4,  Inf,   0,
    rep(0, 4),    0,   0,    0, Inf
  ), byrow = TRUE, nrow = 4 * 2, ncol = 4 * 2)
)
