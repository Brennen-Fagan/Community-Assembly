IslandPreprocess <- function(
  Pool,
  InteractionMatrix,
  Communities, # List containing each Community on each island.
  Populations, # List containing each Population on each island.
  DispersalPool, # Species related dispersal rates
  # Should have length == nrow(Pool). Multiplied by entries of
  DispersalIsland,
  Tolerance,
  Verbose = FALSE
) {
  # Sanity check 1. ############################################################
  stopifnot(length(Communities) == length(Populations))
  
  if (length(DispersalPool) == 1) {
    DispersalPool <- rep(DispersalPool, nrow(Pool))
  }
  
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
  
  # Sanity Check 2. ###########################################################
  if (!((all(
    unlist(lapply(CommunitiesNumeric, length)) ==
    unlist(lapply(PopulationsNumeric, length))
  )))) {
    stop("Entries in Communities and Populations differ.")
  }
  
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
  
  redDisPool <- DispersalPool[redCom]
  
  # Create the dispersal matrix. ###############################################
  
  # Alternative: create multiple sparse diagonal matrices, then r and c bind.
  # Loop through, keep track of which island we are on to add in an empty block.
  
  # Alternative for simple case of just an island, not species, travel matrix.
  # redDisMat <- DispersalRates
  # diag(redDisMat) <- -Matrix::colSums(dispersalMatrix)
  
  # So the question is how to format the diagonals then.
  dispersalDiags <- NULL
  for (i in (length(Communities) - 1):-(length(Communities) - 1)) {
    if (i == 0) next
    # Start in top right band, move to bottom left.
    dispersalDiags <- c(
      dispersalDiags,
      list(c(
        unlist(lapply(
          1:length(Communities),
          function(index, offset, mat, vec) {
            if (offset != 0 &&
                nrow(mat) >= index + offset &&
                index + offset >= 1)
              mat[index, index + offset] * vec
          },
          offset = i, mat = DispersalIsland, vec = redDisPool
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
    k = length(redCom) * c((length(Communities) - 1):1,
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
    epsilon = Tolerance,
    Verbose = Verbose
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
  
  return(list(
    abundance_init = abundance_init,
    redCom = redCom,
    redComs = redComs,
    redPops = redPops,
    redPool = redPool,
    redDisPool = redDisPool,
    redComMat = redComMat,
    parameters = parameters,
    dispersalMatrix = dispersalMatrix
  ))
}