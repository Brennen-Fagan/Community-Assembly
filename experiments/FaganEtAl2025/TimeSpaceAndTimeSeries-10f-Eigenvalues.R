# Introduction: ###############################################################
# On the question of time scale coupling:
# We've taken an approximation about the time scales -- namely that
#       eye_{n \times n} A \approx J
# where eye is the identity matrix (diagonal all 1),
#       A is the interaction matrix, and
#       J is the Jacobian matrix.
# Note it is straightforward to see that the Jacobian looks like
#       J  = diag(x) A + diag(r + A x)
# where x is the abundance vector and r is the intrinsic reproduction/mortality
# rate vector. Note though that if J is evaluated at a fixed point, either
# (r + A x) or x should be 0. If we delete zero population's rows and columns
# we have J \approx diag(x) A. So our approximation is for a unit population.

# We chose this approximation because we know that if the
# immigration/extirpation (external) dynamics are too fast, the local dynamics
# can be overwhelmed resulting in an effectively neutral model, but if the
# external dynamics are too slow we're not only making certain we go to a stable
# fixed point, but we're also spending a lot of time at this stable fixed point
# not doing anything. So we need the two timescales to be coupled somehow.
# Likewise though, we can't survey every Jacobian of every fixed point in the
# system without evaluating them all. Picking x = 0 says only the smallest basal
# should be looked at. This gets more exaggerated with the habitat type
# preferences; matches make the basal species's reproduction rates faster and
# the consumer species's mortality rates slower. Picking x = 1 weights all
# species's interactions evenly, is independent of our experiments, but makes
# the approximation less accurate. Our question here is how much less accurate?

# Strategy:
#   For each Pool-Interaction Matrix combination
#       Calculate the dominant eigenvalue = max(abs(eigen(A)$value))
#       Calculate the summary statistics of r
#   For each simulated abundance time series
#       Calculate the Dynamics
#       Identify the points of "low" change
#       Calculate the dominant eigenvalue = max(abs(eigen(J)$value))
#   Plot
#       Simulation as rows (y), eigenvalues as x,
#       used timescale as a dot,
#       intrinsic r as a square, and
#       a line for the observed dominant (fastest) and slowest eigenvalues.
#   Evaluate
#       How close is our proposed timescale to the observed fixed-point scales?
#       How does this vary based on habitat type and species preferences?

# Parameters: #################################################################
alsoload <- FALSE # if TRUE, try to load all diversity files encountered.
# if FALSE, only try to create new diversity files (and return the outputs).
overwrite <- TRUE

# datfolders <- dir(pattern = "TSTS_Simulations_.+2025-07-30$")
datfolders <- dir(pattern = "TSTS_Simulations_.") # Grab w/&w/o competition.

cargs <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
if (!exists("cores")) {
  if (is.null(cargs) || is.na(cargs)) {
    cores <- 1
  } else {
    cores <- cargs
  }
}

# Libraries: ##################################################################
directory <- '.'
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

library(dplyr)
library(RMTRCode2)
library(parallel)
library(iterators)
library(doParallel)
library(foreach)

source(file.path("R", "checkDominantEigenvalues.R"))

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  
  # R seems to use too much compute despite trying to use 1 core at a time.
  # This should make parallelisation clearer.
  # Suggested environmental variabls from Google Gemini 3.1 Pro.
  parallel::clusterEvalQ(clust, {
    # Environment variables catch standard backends (OpenBLAS, MKL, Apple)
    Sys.setenv(OMP_NUM_THREADS = 1)
    Sys.setenv(OPENBLAS_NUM_THREADS = 1)
    Sys.setenv(MKL_NUM_THREADS = 1)
    Sys.setenv(VECLIB_MAXIMUM_THREADS = 1)
    Sys.setenv(NUMEXPR_NUM_THREADS = 1)
    
    if (requireNamespace("data.table", quietly = TRUE)) {
      # Same idea, but for data.table that tidytable uses.
      data.table::setDTthreads(1)
    }
    
    # RhpcBLASctl as a secondary measure
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
  })
  
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

# Load Data: ##################################################################
datfolders_properties <- strsplit(datfolders, split = "_")
if ( length(datfolders_properties) > 1 &&
     with(list(x = unlist(lapply(datfolders_properties, length))),
          any(x[1] != x)) ){
  stop("Differing folder types implies differing file types.")
}

flag <- datfolders_properties[[1]][1]
if (flag == "TSTS") {
  splitchar <- "_"
} else if (flag == "Data") {
  splitchar <- "-"
} else {
  stop("Folder type not recognized.")
}

poolmats <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Pool"), function(x) {
    names <- load(x)
    return(c(mget(names), "Dir" = dirname(x), "File" = basename(x)))
  })

# Calculations: ###############################################################
allfiles <- dir(datfolders, full.names = TRUE,
                pattern = "(Simulation|Result|Intervention)")

EigenData <- foreach::foreach(
  # id = iterators::iter(1:length(allfiles)),
  # x = iterators::iter(allfiles)
  id = iterators::iter(900:length(allfiles)),
  x = iterators::iter(allfiles[900:length(allfiles)])
  #, .packages = c("dplyr", "RMTRCode2")
) %op% {
  directory <- '.'
  librarypath <- file.path(directory, "Rlibs")
  .libPaths(c(librarypath, .libPaths()))
  library(RMTRCode2); library(dplyr)

  x_properties <- strsplit(basename(x), split = splitchar)
  stopifnot(length(x_properties) == 1)

  filename <- file.path(
    dirname(x),
    if (flag == "TSTS") {
      paste0(c(x_properties[[1]][1],
               "Eigen",
               x_properties[[1]][3:length(x_properties[[1]])]), collapse = "_")
    } else if (flag == "Data") {
      paste0("TSTS_Eigen_",
             paste0(x_properties[[1]][5:length(x_properties[[1]])],
                    collapse = "_"))
    } else {
      paste0("Eigen_", x)
    }
  )

  if(!overwrite && file.exists(filename)) {
    if (alsoload) {
      load(filename)
    }
  } else {
    print(paste(id, filename))
    x_dir <- dirname(x)
    x_poolind <- which(unlist(lapply(poolmats, function(y) y$Dir == x_dir)))
    if(any(x_poolind)) {
      # x_pool <- poolmats[[x_poolind]]$Pool
      DynamicsFunction <- poolmats[[x_poolind]]$DynamicsFunction
      InteractionMatrices <- poolmats[[x_poolind]]$InteractionMatrices
    } else {
      warning(paste("Pool not found, skipping", id, filename)); next()
    }

    # Load result to analyse.
    loaded <- load(x) # names
    stopifnot(length(loaded) == 1)
    loaded <- (get(loaded)) # objects
    # Note: after patch effects, rather than before
    if ("EffectiveReproductionRate" %in% names(loaded$Ellipsis$Affinity)) {
      rprime <- loaded$Ellipsis$Affinity$EffectiveReproductionRate
    } else if ("EffectiveReproductionRateIntervention" %in% names(loaded$Ellipsis$Affinity)) {
      rprime <- loaded$Ellipsis$Affinity$EffectiveReproductionRateIntervention
    } else {
      warning("Reproduction Rate rprime not found."); next()
    }

    # Different handling is needed depending on if arguments are functions.
    # (Could functionalise the calls, but that introduces overhead on each call.)
    if (is.function(rprime)) {
      # Calculate rprime using Parms$Patch (i.e., within patch)
      if (is.function(InteractionMatrices$Mats[[1]])) {
        # Calculate and combine interaction matrices on the fly.
        PerCapitaDynamics <- DynamicsFunction(
          rprime,
          function(t, y, parms) {
            Matrix::bdiag(lapply(
              InteractionMatrices$Mats,
              function(matfunc) {matfunc(t, y, parms)}
            ))
          },
          loaded$NumEnvironments
        )
      }
      else {
        # Just combine the interaction matrices.
        PerCapitaDynamics <- DynamicsFunction(
          rprime,
          Matrix::bdiag(InteractionMatrices$Mats),
          loaded$NumEnvironments
        )
      }
    } else {
      # Treat rprime as constant and explicitly calculated.
      if (is.function(InteractionMatrices$Mats[[1]])) {
        # Calculate and combine interaction matrices on the fly.
        PerCapitaDynamics <- DynamicsFunction(
          rprime,
          function(t, y, parms) {
            Matrix::bdiag(lapply(
              InteractionMatrices$Mats,
              function(matfunc) {matfunc(t, y, parms)}
            ))
          }
        )
      }
      else {
        # Just combine the interaction matrices.
        PerCapitaDynamics <- DynamicsFunction(
          rprime,
          Matrix::bdiag(InteractionMatrices$Mats)
        )
      }
    }

    EigData <- list(
      EigenData = with(loaded, checkAllDominantEigenvalues(
        Abundance = Abundance,                            # Inside loaded
        Events = Events,                                  # Inside loaded
        PerCapitaDynamics = PerCapitaDynamics,            # Assembled from poolmats
        DispersalMatrix = Parameters$DispersalMatrix,     # Inside loaded
        MaxChangeSize = 1e-6 # How much change is "too much" to be @ f.p.
      )),
      EigenInteraction = max(abs(unlist(lapply(
        poolmats[[x_poolind]]$InteractionMatrices$Mats,
        eigen, only.values = TRUE
        )))),
      CharacteristicRate = 1/loaded$ReactionTime, # Should match EigenInteraction.
      ReproductionRate = poolmats[[x_poolind]]$Pool$ReproductionRate, # orig rs
      Ellipsis = loaded$Ellipsis
    )

    if ("ParentRun" %in% names(EigData$Ellipsis))
      EigData$Ellipsis$GrandparentRun <- EigData$Ellipsis$ParentRun
    EigData$Ellipsis$ParentRun <- x

    # So now Diversity contains summary statistics, presence/absence values,
    # and simulation metadata.
    save(EigData, file = filename)
  }
}

if (cores > 1)
  parallel::stopCluster(clust)
