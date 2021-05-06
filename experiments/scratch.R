candidateData$CommunityAbund <- ""
for (r in 1:nrow(candidateData)) {
  candidateData$CommunityAbund[r] <- toString(with(
    candidateData[r, ], {
      invasions <- CommunitySeq[[1]]$Result.Addition
      abundance <- rep(0, Basals + Consumers)
      for (i in invasions) {
        abundance[i] <- abundance[i] + 1
        abundance <- RMTRCode2::quiet(
          deSolve::ode(
            #rootSolve::steady(
            y = abundance,
            times = c(0:10000),
            func = RMTRCode2::GeneralisedLotkaVolterra,
            parms = list(a = mats[[CombnNum]],
                         r = pools[[CombnNum]]$ReproductionRate),
            #method = "runsteady"
          )[10001, -1])
      }
      print(paste("Expected", Communities))
      print(paste("Calculated", toString(which(abundance > 0))))
      print(paste("> Epsilon:", toString(which(abundance > 1E-6))))
      if (Communities == toString(which(abundance > 1E-6))) {
        abundance[abundance > 1E-6]
      } else {
        "Failure"
      }
    }
  ))
}
