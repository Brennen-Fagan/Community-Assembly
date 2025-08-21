# Why? so we can have single argument functions with partials.
repFixed <- function(value = 0.5) {
  force(value)
  function(n) {rep(value, n)}
}
rep_0 <- repFixed(0)
rep_0.25 <- repFixed(0.25)
rep_0.5 <- repFixed()
rep_0.75 <- repFixed(0.75)
rep_1 <- repFixed(1)
