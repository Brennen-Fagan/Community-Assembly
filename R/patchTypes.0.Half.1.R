# Easy ring gradients.
patchTypes.0.Half.1 <- function(n) {
  toString(
    c(0,
      rep(0, floor((n - 2)/3)),
      rep(0.5, ceiling((n - 2)/6)),
      1,
      rep(1, floor((n - 2)/3)),
      rep(0.5, ceiling((n - 2)/6))
    )
  )
}
