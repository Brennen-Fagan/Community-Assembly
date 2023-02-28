CalculateLocalInvasibles_BareBones <- function(
    Abundance, # Not Time, All Environments
    PerCapitaDynamics, # function(t, y, parms)
    Environments,
    ArrivalDensity,
    DispersalMatrix,
    Tolerance = sqrt(.Machine$double.eps),
    EliminationThreshold = 1E-4
) {
  stopifnot(
    (n <- length(Abundance)) == ncol(DispersalMatrix),
    n == nrow(DispersalMatrix)
  )

  # Who are you?
  candidates <- which(
    abs(Abundance - rep(0, n)) < Tolerance
  )

  if (length(candidates) == 0) {return(NULL)}

  # Will the patch accept you? Escape early if it accepts no-one.
  pcd <- PerCapitaDynamics(0, Abundance, list())

  if (!any(pcd[candidates] > 0)) {
    retval <- pcd[candidates] > 0
    names(retval) <- candidates
    return(retval)
  }

  # Prepare for mini-assembly process.
  tempDynamics <- compiler::cmpfun(function(t, y, parms) {
    list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
      as.numeric(
        y * PerCapitaDynamics(t, y, parms)
        # Transport: Dispersal means movement of abundance between nodes.
        + DispersalMatrix %*% y
      )
    )
  })

  # Time steps: more = longer calculation but smoother.
  tempTimings <- (0:10)/1

  invadable <- unlist(lapply(candidates, function(i) {
    if (pcd[i] < 0) return(FALSE)

    testAbund <- Abundance
    testAbund[i] <- testAbund[i] + ArrivalDensity

    retAbund <- deSolve::lsoda(
      y = testAbund,
      times = tempTimings,
      func = tempDynamics
    )

    retAbund[, -1] <- ifelse(retAbund[, -1] <= EliminationThreshold,
                             0, retAbund[, -1])

    # Must grow at least once over the timescale and not die out by the end.
    value <- any(diff(retAbund[, i + 1]) >= 0) &
      all(tail(retAbund[, i + 1]) > 0)

    names(value) <- NULL
    return(value)
  }))

  return(invadable)
}

CalculateLocalInvadables_BareBones <- CalculateLocalInvasibles_BareBones

CalculateLocalInvasibles_KnockOn <- function(
    Abundance, # Not Time, All Environments
    PerCapitaDynamics, # function(t, y, parms)
    Environments, # Deprecated.
    ArrivalDensity,
    DispersalMatrix,
    TimeScale = 2, # Size of Time Scale
    TimeIterations = 5, # Number of Time Scales to Run
    TimeSteps = 10, # Number of Check-In points.
    Tolerance = sqrt(.Machine$double.eps),
    EliminationThreshold = 1E-4
) {
  stopifnot(
    (n <- length(Abundance)) == ncol(DispersalMatrix),
    n == nrow(DispersalMatrix)
  )

  # Who are you?
  # Why abs(x - 0)? Switching to the more obvious Abund - Thresh for testing.
  candidates <- which(
    #abs
    (Abundance - EliminationThreshold
       #rep(0, n)
     ) < Tolerance
  )

  if (length(candidates) == 0) {return(NULL)}

  # Will the patch accept you? Escape early if it accepts no-one.
  pcd <- PerCapitaDynamics(0, Abundance, list())

  if (!any(pcd[candidates] > 0)) {
    outputs <- list(
      candidates = candidates,
      invadable = rep(FALSE, length(candidates)),
      effectAbundance =
        Matrix::Matrix(as.numeric(NA), nrow = length(candidates), ncol = n),
      effectRichness = {
        temp <- rep(as.numeric(NA), n)
        temp[Abundance > EliminationThreshold] <- 0
        do.call(rbind, replicate(
          length(candidates),
          Matrix::Matrix(temp, nrow = 1, ncol = n)
        ))
      },
      effectEstablish = {
        temp <- rep(as.numeric(NA), n)
        do.call(rbind, replicate(
          length(candidates),
          Matrix::Matrix(temp, nrow = 1, ncol = n)
        ))
      }
    )

    return(outputs)
  }

  # Prepare for mini-assembly process.
  tempDynamics <- compiler::cmpfun(function(t, y, parms) {
    list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
      as.numeric(
        y * PerCapitaDynamics(t, y, parms)
        # Transport: Dispersal means movement of abundance between nodes.
        + DispersalMatrix %*% y
      )
    )
  })

  # Time steps: more = longer calculation but smoother.
  tempTimings <- seq(from = 0,
                     to = TimeScale * TimeIterations,
                     length.out = TimeSteps + 1)

  outputs <- lapply(candidates, function(i) {
    if (pcd[i] < 0) return(list(
      candidate = i,
      invadable = FALSE
      ))

    testAbund <- Abundance
    testAbund[i] <- testAbund[i] + ArrivalDensity

    retAbund <- deSolve::lsoda(
      y = testAbund,
      times = tempTimings,
      func = tempDynamics
    )

    retAbund[, -1] <- ifelse(retAbund[, -1] <= EliminationThreshold,
                             0, retAbund[, -1])

    # # Must grow at least once over the timescale and not die out by the end.
    # value <- any(diff(retAbund[, i + 1]) >= 0) &
    #   all(tail(retAbund[, i + 1]) > 0)
    # While I agree with the above in principal, this can be discriminatory
    # against the largest consumers.
    value <- all(tail(retAbund[, i + 1]) > 0)

    names(value) <- NULL

    speciesBase <- which(Abundance > EliminationThreshold)
    speciesNow <- which(retAbund[nrow(retAbund), -1] > EliminationThreshold)
    speciesGain <- speciesNow[!speciesNow %in% speciesBase]
    speciesPers <- intersect(speciesBase, speciesNow)
    speciesLost <- speciesBase[!speciesBase %in% speciesNow]
    species <- rep(as.numeric(NA), n)
    species[speciesGain] <- 1
    species[speciesPers] <- 0
    species[speciesLost] <- -1

    establishBase <- which(
      pcd[,] > 0 &
      Abundance < EliminationThreshold
    )
    establishNow <- which(
      PerCapitaDynamics(0, retAbund[nrow(retAbund), -1], list())[,] > 0 &
      retAbund[nrow(retAbund), -1] < EliminationThreshold
    )

    # Can Now Establish but couldn't before = 1
    # Can Now Establish and Could Before = 0
    # Can't Now Establish but could before = -1 (note includes focal)
    # Already Present or Can't and Couldn't = NA
    # You can differentiate using species that persist, are gained or are lost.
    establish <- rep(as.numeric(NA), n)
    establish[establishNow[!establishNow %in% establishBase]] <- 1
    establish[intersect(establishBase, establishNow)] <- 0
    establish[establishBase[!establishBase %in% establishNow]] <- -1


    return(list(
      candidate = i,
      invadable = value,
      effectAbundance = retAbund[nrow(retAbund), -1] - Abundance,
      effectRichness = Matrix::Matrix(species, nrow = 1, ncol = n),
      effectEstablish = Matrix::Matrix(establish, nrow = 1, ncol = n)
    ))
  })

  # Combine the outputs.
  outputs <- list(
    candidates = unlist(lapply(outputs, function(f) f$candidate)),
    invadable = unlist(lapply(outputs, function(f) f$invadable)),
    effectAbundance = do.call(rbind, lapply(
      outputs, function(f)
        if("effectAbundance" %in% names(f)) {
          f$effectAbundance
        } else {
          Matrix::Matrix(as.numeric(NA), nrow = 1, ncol = n)
        }
    )),
    effectRichness = do.call(rbind, lapply(
      outputs, function(f)
        if("effectRichness" %in% names(f)) {
          f$effectRichness
        } else {
          temp <- rep(as.numeric(NA), n)
          temp[Abundance > EliminationThreshold] <- 0
          Matrix::Matrix(temp, nrow = 1, ncol = n)
        }
    )),
    effectEstablish = do.call(rbind, lapply(
      outputs, function(f)
        if("effectEstablish" %in% names(f)) {
          f$effectEstablish
        } else {
          temp <- rep(as.numeric(NA), n)
          Matrix::Matrix(temp, nrow = 1, ncol = n)
        }
    ))
  )

  return(outputs)
}

CalculateLocalInvadables_KnockOn <- CalculateLocalInvasibles_KnockOn

CalculateLocalInvasibles_KnockOnRunSteady <- function(
    Abundance, # Not Time, All Environments
    PerCapitaDynamics, # function(t, y, parms)
    Environments,
    ArrivalDensity,
    DispersalMatrix,
    Tolerance = sqrt(.Machine$double.eps),
    EliminationThreshold = 1E-4
) {
  stopifnot(
    (n <- length(Abundance)) == ncol(DispersalMatrix),
    n == nrow(DispersalMatrix)
  )

  # Who are you?
  # Why abs(x - 0)? Switching to the more obvious Abund - Thresh for testing.
  candidates <- which(
    #abs
    (Abundance - EliminationThreshold
     #rep(0, n)
    ) < Tolerance
  )

  if (length(candidates) == 0) {return(NULL)}

  # Will the patch accept you? Escape early if it accepts no-one.
  pcd <- PerCapitaDynamics(0, Abundance, list())

  if (!any(pcd[candidates] > 0)) {

    outputs <- list(
      candidates = candidates,
      invadable = rep(FALSE, length(candidates)),
      effectAbundance =
            Matrix::Matrix(as.numeric(NA), nrow = length(candidates), ncol = n),
      effectRichness = {
            temp <- rep(NA, n)
            temp[Abundance > EliminationThreshold] <- 0
            do.call(rbind, replicate(
              length(candidates),
              Matrix::Matrix(temp, nrow = 1, ncol = n)
            ))
          }
    )

    return(outputs)
  }

  # Prepare for mini-assembly process.
  tempDynamics <- compiler::cmpfun(function(t, y, parms) {
    list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
      as.numeric(
        y * PerCapitaDynamics(t, y, parms)
        # Transport: Dispersal means movement of abundance between nodes.
        + DispersalMatrix %*% y
      )
    )
  })

  outputs <- lapply(candidates, function(i) {
    if (pcd[i] < 0) return(list(
      candidate = i,
      invadable = FALSE
    ))

    testAbund <- Abundance
    testAbund[i] <- testAbund[i] + ArrivalDensity

    # We're going to try and run to the steady state.
    # If it fails, we'll test run it, see if anything is removed, and then try
    # again until it works or the test run doesn't change things.
    FLAG <- TRUE

    while(FLAG) {
      capture.output(retAbund <- tryCatch(rootSolve::runsteady(
        y = testAbund,
        func = tempDynamics
      ), error = function(e) {return(FALSE)}))

      if (isTRUE(all.equal(retAbund, FALSE))) {
        temp <- deSolve::lsode(testAbund, 0:1000, tempDynamics)
        testAbund <- temp[nrow(temp), -1]
        testAbund <- ifelse(testAbund > EliminationThreshold, testAbund, 0)
      } else {
        FLAG <- FALSE
      }
    }

    retAbund <- ifelse(retAbund$y <= EliminationThreshold,
                             0, retAbund$y)

    # Must grow at least once over the timescale and not die out by the end.
    value <- (retAbund[i] - ArrivalDensity > 0) # And greater than elim, but...

    names(value) <- NULL

    speciesBase <- which(Abundance > EliminationThreshold)
    speciesNow <- which(retAbund[nrow(retAbund)] > EliminationThreshold)
    speciesGain <- speciesNow[!speciesNow %in% speciesBase]
    speciesPers <- intersect(speciesBase, speciesNow)
    speciesLost <- speciesBase[!speciesBase %in% speciesNow]
    species <- rep(NA, n)
    species[speciesGain] <- 1
    species[speciesPers] <- 0
    species[speciesLost] <- -1

    return(list(
      candidate = i,
      invadable = value,
      effectAbundance = retAbund - Abundance,
      effectRichness = Matrix::Matrix(species, nrow = 1, ncol = n)
    ))
  })

  # Combine the outputs.
  outputs <- list(
    candidates = unlist(lapply(outputs, function(f) f$candidate)),
    invadable = unlist(lapply(outputs, function(f) f$invadable)),
    effectAbundance = do.call(rbind, lapply(
      outputs, function(f)
        if("effectAbundance" %in% names(f)) {
          f$effectAbundance
        } else {
          Matrix::Matrix(as.numeric(NA), nrow = 1, ncol = n)
        }
    )),
    effectRichness = do.call(rbind, lapply(
      outputs, function(f)
        if("effectRichness" %in% names(f)) {
          f$effectRichness
        } else {
          temp <- rep(as.numeric(NA), n)
          temp[Abundance > EliminationThreshold] <- 0
          Matrix::Matrix(temp, nrow = 1, ncol = n)
        }
    ))
  )

  return(outputs)
}



CalculateLocalInvadables_KnockOnRunSteady <- CalculateLocalInvasibles_KnockOnRunSteady
