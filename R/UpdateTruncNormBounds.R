# Mathematica:
# Quantile[
#   TruncatedDistribution[{a, b},
#                         NormalDistribution[\[Mu], \[Sigma]]], {q}]
# > \mu - \sqrt{2} \sigma \text{erfc}^{-1}\left(
#       2 \left(
#           q \left(
#               \frac{1}{2} \text{erfc} \left(
#                   \frac{\mu - b}{\sqrt{2} \sigma}
#               \right) - \frac{1}{2} \text{erfc} \left(
#                   \frac{\mu - a}{\sqrt{2} \sigma}
#               \right)
#           \right) + \frac{1}{2} \text{erfc} \left(
#               \frac{\mu - a}{\sqrt{2} \sigma}
#           \right)
#       \right)
#   \right)
# where erfc, complementary error function erfc(z) = 1 - erf(z),
#       erf, error function erf(z) = 2/\sqrt{\pi} \int_0^z \exp{-t^2} dt

# Okay, so from the documentation of ?pnorm:
## if you want the so-called 'error function'
# erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## (see Abramowitz and Stegun 29.2.29)
## and the so-called 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## and the inverses
# erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)

# From which we can assemble:
UpdateTruncNormBounds <- function(ps, a, b, mu, sigma) {
  stopifnot(all(ps >= 0 & ps <= 1))
  mu - sqrt(2) * sigma * erfcinv(
    2 * (
      ps * (
        1/2 * erfc(
          (mu - b) / (sqrt(2) * sigma)
        ) - 1/2 * erfc(
          (mu - a) / (sqrt(2) * sigma)
        )
      ) + 1/2 * erfc(
        (mu - a) / (sqrt(2) * sigma)
      )
    )
  )
}
