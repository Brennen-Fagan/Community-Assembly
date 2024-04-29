# Double Check that the per capita dynamics function actually yields (roughly?)
# the same as the growth rate that I calculated "by hand" for the nodes.
# Copying and pasting from 6b as necessary.

source(file.path(".", "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(".", "TimeSpaceAndTimeSeries-0-Interventions.R"))

entries <- lapply(strsplit(tag, "_"), strsplit, split = "-")

# appendID <- paste0(
#   # PARAMETERS:
#   interventionPatchDictionaryChoice, "-", # Bundle Inter-Simulation Constants.
#   # Where dynamics would go if necessary.
#   interventionTimeDictionaryChoice, "-",
#   interventionDispersalDictionaryChoice, "-", # Sometimes want to change.
#   interventionDistanceDictionaryChoice
#   , "_",
#   # SEEDS:
#   interventionPatchSeedChoice, "-",
#   # Where dynamics would go if necessary.
#   interventionTimeSeedChoice
# )

as.numeric(entries[[1]][[3]][1]) -> interventionPatchDictionaryChoice
as.numeric(entries[[1]][[3]][2]) -> interventionTimeDictionaryChoice
as.numeric(entries[[1]][[3]][3]) -> interventionDispersalDictionaryChoice
as.numeric(entries[[1]][[3]][4]) -> interventionDistanceDictionaryChoice

# Dictionaries: ###############################################################

# We'll pre-load the pool and patch dynamics. This allows us to infer some
# parameters.
poolMats <- new.env()
poolMats$Pool <- Pool
poolMats$InteractionMatrices <- InteractionMatrices
NumberOfEnvironments <- length(poolMats$InteractionMatrices$Mats)

interventionPatchDictionary <- expand.grid(
  PatchAffinities = c(
    # Detection via if string begins with a numeric or a non-numeric.
    # If numeric, it takes it as a fixed set of affinities.
    # If non-numeric, it attempts to treat the string as a function name.
    # In the latter case, it provides ONLY NumberEnvironments as an argument.
    
    toString(rep(0, NumberOfEnvironments)), # Patches -> {0}
    toString(rep(0.5, NumberOfEnvironments)), # Patches -> {0.5}
    toString(rep(1, NumberOfEnvironments)), # Patches -> {1}
    toString(c(rep(0, NumberOfEnvironments/2),
               rep(1, NumberOfEnvironments/2))), # Patches -> {0, 1} Gradient
    "sample.int.normalized", # Patches -> {0, 1} Unif @ Random
    "patchTypes.0.Half.1", # Patches -> {0, 0.5, 1} Gradient Ring
    "sample.int.3", # Patches -> {0, 0.5, 1} Unif @ Random
    "runifRing", # Patches -> [0, 1] Gradient Ring
    "runif" # Patches -> [0, 1] Unif @ Random
  ),
  InterventionLocation = c(
    # Percentage seems easiest
    NA, # == Random
    1, # Last
    0 # First
  ),
  InterventionPercentage = c(
    0.5,
    1
  ),
  stringsAsFactors = FALSE
)[interventionPatchDictionaryChoice, , drop = FALSE]

interventionTimeDictionary <- data.frame(
  # Time1, Time2; called by eval(str2lang(X)) where X is the string below
  #               and "loaded" is the file that is loaded.
  Time1 = c(
    "median(loaded$Events$Times)",
    "quantile(loaded$Events$Times, p = 0.25)"
  ),
  Time2 = c(
    "1/2 * max(loaded$Events$Times)",
    "quantile(loaded$Events$Times, p = 0.75)"
  ),
  Method = c(# each needs a custom implementation unfortunately!
    "mean",
    "runif"
  ),
  InterventionTimespan = c(
    0 # Instantaneous => Switch
    # Else: Should be numeric > 0, determines timespan for interpolation.
  )
)[interventionTimeDictionaryChoice, ]

# TECH DEBT: Copied from 6a-simulations.R. Should be a common resource in final.
# Furthermore, we may want to add an option to use the old dispersal;
# many of my previous runs were the same pool-patches but different dispersals.
# This implementation forces them to experience a new dispersal, and without a
# proper transition (I don't think I have a dynamic dispersal function?).
interventionDispersalDictionary <- rbind(
  data.frame(Resistance = Inf, Configuration = "None"),
  expand.grid(
    Resistance = 10^c(0:9),
    Configuration = c("Ring", "Line", "Complete")
  ))[ifelse(is.na(interventionDispersalDictionaryChoice),
            1, interventionDispersalDictionaryChoice + 2), ]

interventionDistanceDictionary <- data.frame(
  rhofunction = c( # Take patch
    "rho.2.0.1.euclidean",
    "rho.2.1.2.euclidean",
    "rho.10.1.2.euclidean",
    stringsAsFactors = FALSE
  )
)[interventionDistanceDictionaryChoice, ]

# Instantiate Dispersal Matrix: ###############################################
# Copied from 6a-Simulations.R. Should the configuration "switch" be a function
# in order to make it a common resource as well?
DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
  EnvironmentDistances = with(c(
    interventionDispersalDictionary,
    Environments = NumberOfEnvironments
  ), {
    if (Configuration == "None") {
      DistanceMatrix <- Matrix::sparseMatrix(
        i = Environments, j = Environments, x = 0)
    }
    if (Configuration == "Ring" || Configuration == "Line")
      DistanceMatrix <- Matrix::bandSparse(
        Environments, k = c(-1, 1),
        diagonals = list(rep(Resistance, Environments - 1),
                         rep(Resistance, Environments - 1))
      )
    if (Configuration == "Ring") {
      DistanceMatrix[Environments, 1] <- Resistance
      DistanceMatrix[1, Environments] <- Resistance
    }
    if (Configuration == "Complete") {
      DistanceMatrix <- matrix(Resistance,
                               nrow = Environments,
                               ncol = Environments)
      diag(DistanceMatrix) <- 0
    }
    DistanceMatrix
  }
  ),
  SpeciesSpeeds = poolMats$Pool$Speed
)

##### Post-intervention adjusted intrinsic growth/decay rates: ##############
# TECH DEBT: Copied from 6a-Simulation.R.
rho <- retrieveFunction(interventionDistanceDictionary)

grid <- expand.grid(
  pool = 1:nrow(poolMats$Pool), # Fastest Varying
  patch = 1:nrow(result$Ellipsis$Affinity$PatchAffinitiesIntervention) # Slower Varying.
)
rprime <- result$Ellipsis$Affinity$EffectiveReproductionRateIntervention
rprimeold <- result$Ellipsis$Affinity$EffectiveReproductionRateOld
interventionPatches <- result$Ellipsis$Affinity$PatchInterventions
interventionTime <- result$Ellipsis$Affinity$TimeIntervention

# Need to specify each patch separately with how I've implemented the
# DynamicsFunction. (I recall wanting to have per patch stats.)
# Note, this is inefficient.
if (interventionTimeDictionary$InterventionTimespan == 0) {
  rprimeSwitches <- lapply(1:NumberOfEnvironments, function(i) {
    if (i %in% interventionPatches) {
      switchMatrices(
        rprimeold[
          (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
          ],
        rprime[
          (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
          ],
        switchtime = interventionTime
      )
    } else {
      function(t,...) rprime[
        (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
        ]
    }
  })
  rprimef <- function(t, parms, ...) {
    return(rprimeSwitches[[parms$Patch]](t, parms, ...))
  }
} else {
  rprimeSwitches <- lapply(1:NumberOfEnvironments, function(i) {
    if (i %in% interventionPatches) {
      interpolateMatrices(
        rprimeold[
          (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
          ],
        rprime[
          (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
          ],
        switchtime = interventionTime,
        timespan = interventionTimeDictionary$InterventionTimespan
      )
    } else {
      function(t,...) rprime[
        (i - 1) * nrow(poolMats$Pool) + 1 : nrow(poolMats$Pool)
        ]
    }
  })
  rprimef <- function(t, parms, ...) {
    return(rprimeSwitches[[parms$Patch]](t, parms, ...))
  }
}

interventionPerCapitaDynamics <- with(poolMats, {
  # TECH DEBT: Copied from 6a-simulations.R
  if (is.function(rprimef)) {
    # Calculate rprimef using Parms$Patch
    if (is.function(InteractionMatrices$Mats[[1]])) {
      # Calculate and combine interaction matrices on the fly.
      DynamicsFunction(
        rprimef,
        function(t, y, parms) {
          Matrix::bdiag(lapply(
            InteractionMatrices$Mats,
            function(matfunc) {matfunc(t, y, parms)}
          ))
        },
        NumberOfEnvironments
      )
    }
    else {
      # Just combine the interaction matrices.
      DynamicsFunction(
        rprimef,
        Matrix::bdiag(InteractionMatrices$Mats),
        NumberOfEnvironments
      )
    }
  } else {
    # Treat rprimef as constant and explicitly calculated.
    if (is.function(InteractionMatrices$Mats[[1]])) {
      # Calculate and combine interaction matrices on the fly.
      DynamicsFunction(
        rprimef,
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
      DynamicsFunction(
        rprimef,
        Matrix::bdiag(InteractionMatrices$Mats)
      )
    }
  }
})

Dynamics <- function(t, y, parms) {
  list( # Reaction: PerCapitaDynamics includes interactions and reproduction.
    as.numeric(
      y * interventionPerCapitaDynamics(t, y, parms)
      # Transport: Dispersal means movement of abundance between nodes.
      + DispersalMatrix %*% y
    )
  )
}

timeSimulation <- times[timestep] * result$ReactionTime
timestepResult <- which(result$Abundance[, 1] == timeSimulation)


expected <- 
  Dynamics(timeSimulation, result$Abundance[timestep,-1], parms = list())[[1]] / 
  result$Abundance[timestep,-1]
