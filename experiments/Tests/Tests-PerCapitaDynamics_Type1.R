print("test t")
with(
  list(temp = function(t, ...) {t}),
  with(
    list(tempfunc = RMTRCode2::PerCapitaDynamics_Type1(temp, temp, 3)),
    stopifnot(
      tempfunc(1, c(1, 1, 1)) == matrix(2, nrow = 1, ncol = 3),
      tempfunc(4, c(1, 1, 1)) == matrix(8, nrow = 1, ncol = 3),
      tempfunc(1, c(1, 2, 3)) == matrix(c(2, 3, 4), nrow = 1, ncol = 3),
      tempfunc(4, c(1, 2, 3), parms = list(a = 1)
      ) == matrix(c(8, 12, 16), nrow = 1, ncol = 3))
  ))

print("test y")
with(
  list(temp = function(y, ...) {y}),
  with(
    list(tempfunc = RMTRCode2::PerCapitaDynamics_Type1(1, temp, 1)),
    stopifnot(
      tempfunc(1, c(1, 1, 1)) == matrix(4, nrow = 1, ncol = 1),
      tempfunc(4, c(1, 1, 1)) == matrix(4, nrow = 1, ncol = 1),
      tempfunc(1, c(1, 2, 3)) == matrix(15, nrow = 1, ncol = 1),
      tempfunc(4, c(1, 2, 3), parms = list(a = 1)
      ) == matrix(15, nrow = 1, ncol = 1))
  ))


print("test parms")
with(
  list(temp1 = function(parms, ...) {parms$a},
       temp2 = function(parms, ...) {parms$a * parms$b}),
  with(
    list(tempfunc1 = RMTRCode2::PerCapitaDynamics_Type1(temp1, temp1, 3),
         tempfunc2 = RMTRCode2::PerCapitaDynamics_Type1(temp2, temp2, 3),
         parms1 = list(a = 1, b = 2, c = 3),
         parms2 = list(a = 2, b = 2, c = 0)),
    stopifnot(
      tempfunc1(1, c(1, 1, 1), parms1) == matrix(2, nrow = 1, ncol = 3),
      tempfunc1(4, c(1, 1, 1), parms2) == matrix(4, nrow = 1, ncol = 3),
      tempfunc1(1, c(1, 2, 3), parms1) == matrix(c(2, 3, 4), nrow = 1, ncol = 3),
      tempfunc1(4, c(1, 2, 3), parms2) == matrix(c(4, 6, 8), nrow = 1, ncol = 3),
      tempfunc2(1, c(1, 1, 1), parms1) == matrix(4, nrow = 1, ncol = 3),
      tempfunc2(4, c(1, 1, 1), parms2) == matrix(8, nrow = 1, ncol = 3),
      tempfunc2(1, c(1, 2, 3), parms1) == matrix(c(4, 6, 8), nrow = 1, ncol = 3),
      tempfunc2(4, c(1, 2, 3), parms2) == matrix(c(8, 12, 16), nrow = 1, ncol = 3))
  ))

print("test all")
with(
  list(temp = function(t, y, parms) {
    print(paste("y is here too:", y, collapse = "!"))
    t/parms$a + parms$b
  }),
  with(
    list(tempfunc = RMTRCode2::PerCapitaDynamics_Type1(temp, temp, 3),
         parms1 = list(a = 1, b = 2, c = 3),
         parms2 = list(a = 2, b = 2, c = 0),
         captures = list(),
         returns = list()),
    {
      captures[[1]] <- capture.output(returns[[1]] <-
        tempfunc(1, c(1, 1, 1), parms1) == matrix(6, nrow = 1, ncol = 3)
      )
      captures[[2]] <- capture.output(returns[[2]] <-
        tempfunc(4, c(1, 1, 1), parms1) == matrix(12, nrow = 1, ncol = 3)
      )
      captures[[3]] <- capture.output(returns[[3]] <-
        tempfunc(1, c(1, 2, 3), parms1) == matrix(c(6, 9, 12), nrow = 1, ncol = 3)
      )
      captures[[4]] <- capture.output(returns[[4]] <-
        tempfunc(4, c(1, 2, 3), parms1) == matrix(c(12, 18, 24), nrow = 1, ncol = 3)
      )
      captures[[5]] <- capture.output(returns[[5]] <-
        tempfunc(1, c(1, 1, 1), parms2) == matrix(5, nrow = 1, ncol = 3)
      )
      captures[[6]] <- capture.output(returns[[6]] <-
        tempfunc(4, c(1, 1, 1), parms2) == matrix(8, nrow = 1, ncol = 3)
      )
      captures[[7]] <- capture.output(returns[[7]] <-
        tempfunc(1, c(1, 2, 3), parms2) == matrix(c(5, 7.5, 10), nrow = 1, ncol = 3)
      )
      captures[[8]] <- capture.output(returns[[8]] <-
        tempfunc(4, c(1, 2, 3), parms2) == matrix(c(8, 12, 16), nrow = 1, ncol = 3)
      )
      stopifnot(all(unlist(returns)),
                unlist(captures) == rep(rep(
                  c("[1] \"y is here too: 1!y is here too: 1!y is here too: 1\"",
                    "[1] \"y is here too: 1!y is here too: 2!y is here too: 3\""),
                  each = 4
                ), times = 2))
    }
  ))

print("test static")
with(
  list(temp = function(...) {exp(1)}),
  with(
    list(tempfunc = RMTRCode2::PerCapitaDynamics_Type1(temp, temp, 3)),
    stopifnot(
      tempfunc(1, c(1, 1, 1)) == matrix(2*exp(1), nrow = 1, ncol = 3),
      tempfunc(4, c(1, 1, 1)) == matrix(2*exp(1), nrow = 1, ncol = 3),
      tempfunc(1, c(1, 2, 3)) == matrix(c(2, 3, 4)*exp(1), nrow = 1, ncol = 3),
      tempfunc(4, c(1, 2, 3)) == matrix(c(2, 3, 4)*exp(1), nrow = 1, ncol = 3))
  ))


print("test misspecified")
with(
  list(temp1 = function() {1},
       temp2 = function(t) {1},
       temp3 = function(y) {1},
       temp4 = function(t, y) {1},
       temp5 = function(t, parms) {1},
       temp6 = function(y, parms) {1}),
  {
    tools::assertError(RMTRCode2::PerCapitaDynamics_Type1(temp1, 1, 3))
    #"isvalid(ReproductionRate) is not TRUE",
    tools::assertError(RMTRCode2::PerCapitaDynamics_Type1(1, temp2, 3))
    #"isvalid(InteractionMatrices) is not TRUE",
    tools::assertError(RMTRCode2::PerCapitaDynamics_Type1(temp3, 1, 3))
    #"isvalid(ReproductionRate) is not TRUE",
    tools::assertError(RMTRCode2::PerCapitaDynamics_Type1(temp4, 1, 3))
    #"isvalid(ReproductionRate) is not TRUE",
    tools::assertError(RMTRCode2::PerCapitaDynamics_Type1(temp5, 1, 3))
    #"isvalid(ReproductionRate) is not TRUE",
    tools::assertError(RMTRCode2::PerCapitaDynamics_Type1(temp6, 1, 3))
    #"isvalid(ReproductionRate) is not TRUE"
  }
)
