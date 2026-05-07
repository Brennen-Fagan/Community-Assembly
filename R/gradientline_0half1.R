gradientline_0half1 <- function(n) {
  left <- rep(0, floor(n / 3)); right <- rep(1, floor(n/3))
  return(c(left, rep(0.5, n - length(left) - length(right)), right))
}
