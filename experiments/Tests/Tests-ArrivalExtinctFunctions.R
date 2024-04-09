# Tests-ArrivalExtinctFunctions.R

# Standard Opening: ###########################################################
library(RMTRCode2)

if (!exists("testseed")) {
  testseed <- runif(1) * 1E8
  if (exists(".Random.seed")) {
    old.seed <- .Random.seed
  }
  set.seed(testseed)
  print(testseed)
}

# Trivial Tests: ##############################################################
stopifnot(identical(ArrivalFUN_Example, ExtinctFUN_Example))
stopifnot(identical(ArrivalFUN_Example2, ExtinctFUN_Example2))

# Stochastic Tests: ###########################################################
## 1: ArrivalFun_Example: #####################################################
testrate <- 10^runif(1, -3, 3)
testlength <- 100000
testtimes <- ArrivalFUN_Example(testlength, testrate)

stopifnot(testlength == length(testtimes))

estrate <- MASS::fitdistr(diff(testtimes), "exponential")

stopifnot(abs((estrate$estimate - testrate)/testrate) < 0.05)

## 2: ArrivalFun_Example2: ####################################################
# This version aims to get more then enough events, then stop just shy of the
# implied duration.
testtimes2 <- ArrivalFUN_Example2(testlength, testrate)

# Stops short.
stopifnot(testtimes2[length(testtimes2)] < testlength / testrate)
# Is not unusually far away from the end.
stopifnot(
  pexp(testlength / testrate - testtimes2[length(testtimes2)], testrate) > 
    1 - 1 / testlength
)

# Has about the right rate.
estrate2 <- MASS::fitdistr(diff(testtimes2), "exponential")
stopifnot(abs((estrate2$estimate - testrate)/testrate) < 0.05)

# Standard Closing: ###########################################################
print("Success.")
print("Testseed:")
print(testseed)

if (exists("old.seed")) set.seed(old.seed)