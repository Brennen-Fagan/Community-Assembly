evensplit <- function(values = c(0, 1)) {
  force(values)
  function(n) {
    c(rep(values, times = floor(n / length(values))),
      if (n %% length(values) != 0) {
        values[1:(n %% length(values))]
      })
  }
}
evensplit_01 <- evensplit()
evensplit_0.51 <- evensplit(c(0.5, 1))
