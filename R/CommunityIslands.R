# Islands and Interactions
# So for our initial interaction, we need to have a graph of patches.
# Each patch contains an ecosystem.
# Patches are effectively homogeneous at this point.
# Heterogeneity is hard to decide how to introduce, given that the obvious
# answer would be to modify the effective species carrying capacities, but
# those are implicit in the model, appearing due to reproduction rates and
# intra-species competition.
# Patches also need their dynamics under which the community is run.
# Then there are dispersal rates which are some function of species and distance
# from the patch perspective.
# From this perspective, islands are closely related sets of patches in that
# their adjacent distances are 0 (which does not mean instantaneous travel!).
#       Patch                        Patch
#   |-----------|                |-----------|
#   | Community |    Dispersal   | Community |
#   | Dynamics  | -------------> | Dynamics  |
#   |           | <------------- |           |
#   |           |                |           |
#   |-----------|                |-----------|
#
#
# Inputs: Species Pool, Species Interaction Matrix, which imply dynamics
# Species present on each patch, Patch distances

CsvRowSplit <- function(csv) {
  return(as.numeric(unlist(strsplit(csv, split = ", "))))
}

FindSteadyStateFromEstimate <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
  Tolerance = 1E-1
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }

  if (is.character(Populations)) {
    init <- CsvRowSplit(Populations)
  } else {
    init <- Populations
  }
  init_min <- min(init)
  init_max <- max(init)

  anyZeroOrNotSame <- TRUE # Set T to make at least one look.
  epsilon <- Tolerance

  # Run and check to see if anyone dies (bad) or changes (not at steady).
  while (anyZeroOrNotSame) {
    init_old <- init

    init <- rootSolve::steady(
      y = init,
      func = Dynamics,
      parms = list(a = InteractionMatrix[com, com],
                   r = Pool$ReproductionRate[com],
                   epsilon = epsilon),
      positive = TRUE
    )$y
    print(init)

    if (any(init < epsilon)) {
      # Someone died, reset to random location.
      print("Died")
      anyZeroOrNotSame <- TRUE
      init[init < epsilon] <- runif(
        n = sum(init < epsilon),
        min = init_min,
        max = init_max
      )
    } else if (any(round(init / init_old, 1) != 1)) {
      # Not done, keep going.
      print("Changed")
      anyZeroOrNotSame <- TRUE
    } else {
      anyZeroOrNotSame <- FALSE
    }
  }

  return(init)
}

Productivity <- function(
  Pool,
  InteractionMatrix,
  Community,
  Populations,
  Dynamics = RMTRCode2::GeneralisedLotkaVolterra,
) {
  if (is.character(Community)) {
    com <- CsvRowSplit(Community)
  } else {
    com <- Community
  }

  if (is.character(Populations)) {
    pop <- CsvRowSplit(Populations)
  } else {
    pop <- Populations
  }

  comMatPos <- InteractionMatrix[com, com]; comMatPos[comMatPos < 0] <- 0
  poolRepPos <- Pool$ReproductionRate[com]; poolRepPos[poolRepPos < 0] <- 0

  parameters <- list(
    a = comMatPos,
    r = poolRepPos
  )

  return(
    sum(Dynamics(0, pop, parameters)[[1]] * pop / sum(pop))
  )
}

IslandLotkaVolterra <- function(t, y, parms) {
  with(as.list(parms), {
    list(as.numeric(y * (r + a %*% y) + d %*% y))
    # as.numeric since the solver doesn't know what Matrix::Matrices are.
  })
}

IslandDynamics <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  DispersalPool, # Species related dispersal rates, which are multiplied by
  DispersalIsland, # Island related dispersal rates. Is a matrix, row = to.
  Dynamics = IslandLotkaVolterra,
  Tolerance = 1E-1,
  Times = seq(from = 0,
              to = 2E4,
              by = 5E2)
  #DispersalRates, # For now, a matrix, to-from notation, for each island.
  #(NOT species on each island, since that requires the user knowing things in advance...)
  #DispersalRates, # List of matrices: column = species, row = (TO other) island, entry = travel rate
) {
  # Sanity checks. #############################################################
  stopifnot(length(Communities) == length(Populations))

  stopifnot(all(
    unlist(lapply(Communities, length)) == unlist(lapply(Populations, length))
  ))

  # Reduce to necessary information. ###########################################
  # Total list of species present.
  # Make sure that formatting is handled.
  CommunitiesNumeric <- lapply(Communities, function(com) {
    if (is.character(com)) {
      CsvRowSplit(com)
    } else if (is.numeric(com)) {
      com
    } else {
      stop("Community is not numeric or string.")
    }
  }
  )

  PopulationsNumeric <- lapply(Populations, function(pop) {
    if (is.character(pop)) {
      CsvRowSplit(pop)
    } else if (is.numeric(pop)) {
      pop
    } else {
      stop("Population is not numeric or string.")
    }
  }
  )

  redCom <- sort(unique(unlist(CommunitiesNumeric))) # Uses recursive unlisting

  redPool <- Pool[redCom,]
  redComMat <- InteractionMatrix[redCom, redCom]

  redComs <- lapply(Communities, function(strCom) {
    which(redCom %in% CsvRowSplit(strCom))
  })

  redPops <- PopulationsNumeric # Unify language, all formatting done.

  # Sanity Check: Lengths should match.
  stopifnot(all(
    unlist(lapply(redComs, length)) == unlist(lapply(redPops, length))
  ))

  # Create the dispersal matrix. ###############################################

  # Alternative: create multiple sparse diagonal matrices, then r and c bind.
  # Loop through, keep track of which island we are on to add in an empty block.

  # Alternative for simple case of just an island, not species, travel matrix.
  # redDisMat <- DispersalRates
  # diag(redDisMat) <- -Matrix::colSums(dispersalMatrix)

  # So the question is how to format the diagonals then.
  dispersalDiags <- NULL
  for (i in length(Communities):-length(Communities)) {
    if (i == 0) next
    # Start in top right band, move to bottom left.
    dispersalDiags <- c(
      dispersalDiags,
      list(c(
        unlist(lapply(
          1:length(Communities),
          function(index,offset,mat,vec) {
            if (offset != 0 &&
                nrow(mat) >= index + offset &&
                index + offset >= 1)
              mat[index, index + offset] * vec
            },
          offset = i, mat = DispersalIsland, vec = DispersalPool
        ))
      ))
      )

  }

  # Amount of travel from j to i is d[i,j]
  # Amount of gain to i from j is d[i,j]
  # Amount of travel from i is d[i,i]
  # This way, we can write the change in y from travel as d %*% y
  # Characterising this as proportions of a population, but that is assuming a
  # normalisation that I do not think is strictly necessary.
  # The matrix is sparse and has colsum = 0, diag < 0, offdiag >= 0.
  dispersalMatrix <- Matrix::bandSparse(
    n = length(redCom) * length(Communities),
    k = length(redCom) * c(1:(length(Communities) - 1),
                           -(1:(length(Communities) - 1))),
    diagonals = dispersalDiags# c(
      # list(c(rep(0.0001, length(redCom)),   # Island 2 -> Island 1
      #        rep(0.0001, length(redCom)))), # Island 3 -> Island 2
      # list(rep(0, length(redCom))),         # Island 3 -> Island 1
      # list(c(rep(0.0001, length(redCom)),   # Island 1 -> Island 2
      #        rep(0.0001, length(redCom)))), # Island 2 -> Island 3
      # list(rep(0, length(redCom)))          # Island 1 -> Island 3
    # )
  )
  stopifnot(all(dispersalMatrix >= 0))
  diag(dispersalMatrix) <- -Matrix::colSums(dispersalMatrix)
  stopifnot(all(Matrix::colSums(dispersalMatrix) == 0))

  parameters <- list(
    r = rep(redPool$ReproductionRate, length(Communities)),
    a = Matrix::bdiag(rep(list(redComMat), length(Communities))),
    d = dispersalMatrix,
    epsilon = Tolerance
  )

  # Technically, a bit of extra work being done here since we already copied the
  # populations per island. The result is somewhat more readable though.
  abundance_init <- unlist(lapply(
    seq_along(redComs),
    function(i, com, pop, numPops) {
      # Interlace 0's with population values
      k <- 1
      retval <- rep(0, numPops)
      for (j in 1:numPops) {
        if (j %in% com[[i]]) {
          retval[j] <- pop[[i]][k]
          k <- k + 1
        }
      }
      return(retval)
    },
    com = redComs,
    pop = redPops,
    numPops = length(redCom)
  ))

  abundance <- deSolve::ode(
    abundance_init,
    times = Times,
    func = dynSys,
    parms = parameters,
    events = list(func = function(t, y, parms) {
      y[y < parms$epsilon] <- 0
      y
    }, time = ts)
  )

  return(abundance)
}
