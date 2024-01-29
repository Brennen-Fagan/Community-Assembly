library(RMTRCode2)
library(dplyr)
library(Matrix)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

# Parameters: ##################################################################
Species <- c(Basal = 34, Consumer = 66) * 2
Environments <- 10
EventsEach <- Environments * ceiling(sum(Species) * (log(sum(Species) + 0))) *
  # 1 # DEFAULT
  2 # Data_2024-01-29; Binary Patch Niche. Make sure we're extra burnt in.

EventRateModifiers <- c(1, 1) # Immigration, Extirpation

PoolPatchNiche <- TRUE
PoolPatchNicheSplit <- c(0, 0.5, 0.5) # Core, Niche1, Niche2, ...
#                        Core is always present (Generalist).
stopifnot(!PoolPatchNiche ||
            sum(PoolPatchNicheSplit) == 1)
AdjustImmigration <- TRUE #
stopifnot(!AdjustImmigration ||
            all(PoolPatchNicheSplit[-1] == PoolPatchNicheSplit[2]))
PoolPatchNicheIntervention <- TRUE

LMParameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
LMLogBodySize <- c(-2, -1, -1, 0)

PerIslandDistance <- 10^c(Inf) # 10^5 # Inf # 10^0
SpeciesSpeeds <- 1
Space <- match.arg("Ring", c("None", "Ring", "Line", "Full"))

EliminationThreshold <- 10^-4 # Below which species are removed from internals
ArrivalDensity <- EliminationThreshold * 4 * 10 ^ 3 # Traill et al. 2007
ExtinctionProportion <- 1

MaximumTimeStep <- 1 # Maximum time solver can proceed without elimination.
BetweenEventSteps <- 10 # Number of steps to reach next event to smooth.

CalculatePoolAndMatrices <- TRUE
dir <- paste0("Data_", "2024-01-29")#Sys.Date())
# getSrcDirectory(function(){})

if (!dir.exists(dir)) {
  dir.create(dir, showWarnings = FALSE)
}

# > runif(3) * 1e8
# [1] 11365664 91994571 20423344 # Data_2023-07-06
# [1] 65566924 64305636 14447307 # Data_2023-09-22
# [1] 71113291 29907014 76606233 # Data_2023-09-23
# [1] 53606086 70944574 40035408 # Data_2023-09-24
# [1]  5005152 70117044 42048254 # Data_2023-09-25
# [1] 77776934  9954265 47259175 # Data_2023-09-26
# [1] 38427042 12032489 28665115 # Data_2024-01-17; Low PValues (0, 0.16) Matrix
# [1] 75027622 64713671 21957601 # Data_2024-01-18; High P (0.84, 1) Matrix
# [1] 54497638 90525137 12496702 # Data_2024-01-29; Binary Patch Niche
seeds <- c(
  # 11365664, 91994571, 20423344 # Data_2023-07-06
  # 65566924, 64305636, 14447307 # Data_2023-09-22
  # 71113291, 29907014, 76606233 # Data_2023-09-23
  # 53606086, 70944574, 40035408 # Data_2023-09-24
  #  5005152, 70117044, 42048254 # Data_2023-09-25
  # 77776934,  9954265, 47259175 # Data_2023-09-26
  # 38427042, 12032489, 28665115 # Data_2024-01-17
  # 75027622, 64713671, 21957601 # Data_2024-01-18
  54497638, 90525137, 12496702, 34126575, 87083934 # Data_2024-01-29
)
PoolSeed <- seeds[1]
EnvironmentSeed <- seeds[2]
HistorySeed <- seeds[3]
if (PoolPatchNiche) {
  NicheSeed <- seeds[4]
  if (PoolPatchNicheIntervention) {
    InterventionSeed <- seeds[5]
  }
}

ConstrainTruncNormPs <- # Need to be NULL or two ordered probabilities.
  NULL # DEFAULT
# c(0, pnorm(-1)) # Data_2024-01-17; -1 Std. Dev. or more Extreme.
# c(pnorm(+1), 1) # Data_2024-01-18; +1 Std. Dev. or more Extreme.

# Setup: #######################################################################

# https://stats.stackexchange.com/a/313138
# Note the edit to the y.dual. Possible change from diag's 2017 behaviour?
complement_whuber <- function(y, rho, x, threshold=1e-12) {
  #
  # Process the arguments.
  #
  if(!is.matrix(y)) y <- matrix(y, ncol=1)
  d <- ncol(y)
  n <- nrow(y)
  y <- scale(y, center=FALSE) # Makes computations simpler
  if (missing(x)) x <- rnorm(n)
  #
  # Remove the effects of `y` on `x`.
  #
  e <- residuals(lm(x ~ y))
  #
  # Calculate the coefficient `sigma` of `e` so that the correlation of
  # `y` with the linear combination y.dual %*% rho + sigma*e is the desired
  # vector.
  #
  y.dual <- with(
    svd(y),
    (n-1)*u %*%
      diag(ifelse(d > threshold, 1/d, 0), # Failed in the ncol(y) == 1 case.
           nrow = ncol(u), ncol = ncol(v)) %*%
      t(v))
  sigma2 <- c((1 - rho %*% cov(y.dual) %*% rho) / var(e))
  #
  # Return this linear combination.
  #
  if (sigma2 >= 0) {
    sigma <- sqrt(sigma2)
    z <- y.dual %*% rho + sigma*e
  } else {
    warning("Correlations are impossible.")
    z <- rep(0, n)
  }
  return(z)
}

# TODO: Could be better by not providing the entire pool.
# Works better the larger the pool is. More entries entails more flexibility.
calculateNumericNiche <- function( # Wrapper for complement_whuber.
  # NOTE: FUNCTION USES LINEARITY; PROBABLY DOESN'T MAKE SENSE FOR CATEGORICAL
  pool, # Pool to be mutated.
  nicheValues, # Numeric values to modify. No required distribution.
  correlations, # intensity of *linear* relationship with pool.
  correlationColumns # Columns to use for the correlations. Characters.
) {# Returns the pool with a new column

  stopifnot(length(correlationColumns) == length(correlations))

  niches <- complement_whuber(
    pool[, correlationColumns], correlations, nicheValues
  )

  return(niches)
}

# calculate rather than mutate because we're not directly adding to the Pool.
calculateCategoricalNiche <- function(
  pool, # Pool to be mutated.
  functions, # Output probabilities for each category of the trait. List.
  correlationColumns, # Columns input for each function.
  seed
) {# Returns the pool with a new column

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    set.seed(seed)
  }

  probs <- matrix(unlist(lapply(
    seq_along(functions),
    function(i, pool) functions[[i]](pool[, correlationColumns]),
    pool = pool
  )), ncol = length(functions)) # populates down columns.
  # Each row will be a species.

  niches <- apply(probs, MARGIN = 1, FUN = function(ps) {
    which(rmultinom(1, 1, ps) == 1)
  })

  if (!is.null(names(functions))) {
    niches <- names(functions)[niches]
  }

  if (!is.null(seed)) {
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }

  return(niches)
}

## Pools and Interaction Matrices: #############################################
if (CalculatePoolAndMatrices) {
  Pool <- RMTRCode2::LawMorton1996_species(
    Basal = Species[1],
    Consumer = Species[2],
    Parameters = LMParameters,
    LogBodySize = LMLogBodySize,
    seed = PoolSeed
  )

  if (PoolPatchNiche) {
    PoolPatchNicheFunctions <- lapply(
      PoolPatchNicheSplit, function(p) function(x) {rep(p, nrow(x))}
    )
    names(PoolPatchNicheFunctions) <-
      c("Core", paste0("Niche", 1:(length(PoolPatchNiche) - 1)))

    Pool <- Pool %>% dplyr::mutate(
      Niche_Cat = calculateCategoricalNiche(
        Pool,
        functions = PoolPatchNicheFunctions,
        correlationColumns = c("ID", "Size", "Type"), # for later usage.
        seed = NicheSeed
      )
    )
  }

  InteractionMatrices <- RMTRCode2::CreateEnvironmentInteractions(
    Pool = Pool, NumEnvironments = Environments,
    ComputeInteractionMatrix = RMTRCode2::LawMorton1996_CommunityMat,
    Parameters = LMParameters,
    EnvironmentSeeds = EnvironmentSeed,
    ConstrainP = ConstrainTruncNormPs
  )
  save(Pool, InteractionMatrices,
       file = file.path(dir, paste0(
         "MNA-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
} else {
  load(file = file.path(dir, paste0(
    "MNA-ExampleOutcome-PoolMats-Env", Environments, ".RData")))
}

## Events: #####################################################################

# Note: eigenvalues of block matrices are the eigenvalues of the blocks.
CharacteristicRate <- max(unlist(lapply(
  InteractionMatrices$Mats, function(m) {abs(eigen(m)$values)}
)))

ArrivalRateCorrection <- 1
if (PoolPatchNiche && AdjustImmigration) {
  # Multiplier to increase rate to compensate for the decreased eligibles.
  # If we have (e.g.) 1/4 ineligible, then rate 3/4th as high.
  # So we correct by dividing by the effective eligible pool.
  ArrivalRateCorrection <- 1 / (sum(PoolPatchNicheSplit[1:2]))
}

Events <- RMTRCode2::CreateAssemblySequence(
  Species = sum(Species),
  NumEnvironments = Environments,
  ArrivalEvents = EventsEach * EventRateModifiers[1] * ArrivalRateCorrection,
  ArrivalRate = CharacteristicRate * EventRateModifiers[1] * ArrivalRateCorrection,
  ArrivalFUN = RMTRCode2::ArrivalFUN_Example2,
  ExtinctEvents = EventsEach * EventRateModifiers[2],
  ExtinctRate = CharacteristicRate * EventRateModifiers[2],
  ExtinctFUN = RMTRCode2::ExtinctFUN_Example2,
  HistorySeed = HistorySeed
)

print(combinations <-
        table(Events$Events$Species,
              Events$Events$Environment,
              Events$Events$Type))
if(any(combinations == 0)) {warning(
  "Exists a species which doesn't immigrate to an environment."
)}

## Dynamics: ###################################################################

IntMat <- Matrix::bdiag(InteractionMatrices$Mats)
PerCapitaDynamics <- RMTRCode2::PerCapitaDynamics_Type1(
  Pool$ReproductionRate, IntMat,
  NumEnvironments = Environments
)

records <- foreach::foreach(
  dist = iterators::iter(PerIslandDistance),
  # .export = c(
  #   "Pool", "Environments",
  #   "InteractionMatrices",
  #   "Events",
  #   "PerCapitaDynamics",
  #   "EliminationThreshold",
  #   "ArrivalDensity",
  #   "ExtinctionProportion",
  #   "MaximumTimeStep",
  #   "BetweenEventSteps"),
  .export = "reprate",
  .packages = "RMTRCode2"
) %do% {
  print(paste(dist, "in"))
  # table(Pool)
  pool <- data.frame(n = 1:sum(Species))#Hack
  # print(paste(dist, "hack"))
  ### Spatial/Dispersal: #########################################################
  if (Space == "None") {
    DistanceMatrix <- Matrix::sparseMatrix(
      i = Environments, j = Environments, x = 0)
  }
  if (Space == "Ring" || Space == "Line")
    DistanceMatrix <- Matrix::bandSparse(
      Environments, k = c(-1, 1),
      diagonals = list(rep(dist, Environments - 1),
                       rep(dist, Environments - 1))
    )
  if (Space == "Ring") {
    DistanceMatrix[Environments, 1] <- dist
    DistanceMatrix[1, Environments] <- dist
  }
  if (Space == "Grid") {
    # Given matrix(1:4, nrow = 2), trying 1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4.
    # I.e. matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0), nrow = 4)
    # Use divisor closest to but <= square root for number of rows.
  }
  if (Space == "Full") {
    DistanceMatrix <- matrix(1, nrow = Environments, ncol = Environments)
    diag(DistanceMatrix) <- 0
  }

  DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
    EnvironmentDistances = DistanceMatrix,
    SpeciesSpeeds = rep(SpeciesSpeeds, nrow(pool))
  )

  # print(paste(dist, "disp mat"))
  ## Run: #######################################################################
  print(Sys.time())
  if (exists("MultipleNumericalAssembly_Dispersal")) {
    theFun <- MultipleNumericalAssembly_Dispersal
  } else {
    theFun <- RMTRCode2::MultipleNumericalAssembly_Dispersal
  }
  # print(paste(dist, "fun"))
  print(record <- Sys.time())

  # !!!Technical Debt!!!
  # TODO: Implement correctly by creating a patch dataframe structure
  # that is used inside (a derivative of) MNA_Dispersal and handled by
  # EliminationAndNeutralEvents. The new versions should be able to
  # check patch specific identity and refuse it, or to have the PerCapDyn
  # create a penalty or other functional form to interact with continuous
  # and categorical niches and how well species niche matches patch niche.
  #
  # Instead, what we're doing here is to create a copy of the events,
  # provide one to the simulation and remove the success column from the other,
  # join the success column after the simulation to the preserved one,
  # and then put the preserved one in place of the actually used.
  # That way, the one provided by the simulation can have entries removed
  # (due to mismatch between categorical pool and patch niches).
  # In the process, we'll need half full time on the characteristic scale
  # in order to set when the intervention takes place.

  EventsUnfiltered <- Events
  timeMax <- max(Events$Events$Time) * CharacteristicRate
  timeSwitch <- (timeMax - 1000) / 2 + 1000 # 1000 as a techdebt burn-in

  if (PoolPatchNicheIntervention) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    # Pick a random patch as control, rest as experiment.
    # is adding 1:5 (or w/e) okay? Yes, we're assuming contiguous patches.
    control <- ((sample.int(Environments, 1) + 1:(Environments / 2) ) %% Environments) + 1
    # print(paste(bootstrapID, ":", toString(control)))
    experiment <- c(1:Environments)[!c(1:Environments) %in% control]
    # print(paste(bootstrapID, ":", toString(experiment)))
    if (exists("oldSeed")) {
      .Random.seed <- oldSeed
    }
  } else {
    control <- 1:Environments
  }

  patchesIdentities <- do.call(rbind, lapply(1:Environments, function(i) {
    data.frame(
      Patch = i,
      Type = if(i %in% control) 1 else c(1, 2),
      TimeMin = if(i %in% control) -Inf else c(-Inf, timeSwitch),
      TimeMax = if(i %in% control) Inf else c(timeSwitch, Inf)
    )
  }))

  Events$Events <- Events$Events %>% dplyr::left_join(
    patchesIdentities, by = c("Environment" = "Patch"), suffix = c("", "_Patch")
  ) %>% dplyr::filter(# Reduce initial over-joining.
    Times >= TimeMin,
    Times <= TimeMax
  ) %>% dplyr::left_join(
    Pool, by = c("Species" = "ID"), suffix = c("", "_Pool")
  ) %>% dplyr::filter(
    Type_Patch == Niche_Cat
  ) %>% dplyr::select(
    Times, Species, Environment, Type, Success
  )

  result <- theFun(
    pool, NumEnvironments = Environments,
    InteractionMatrices = InteractionMatrices,
    Events = Events,
    PerCapitaDynamics = PerCapitaDynamics,
    DispersalMatrix = DispersalMatrix,
    EliminationThreshold = EliminationThreshold,
    ArrivalDensity = ArrivalDensity,
    ExtinctionProportion = ExtinctionProportion,
    MaximumTimeStep = MaximumTimeStep,
    BetweenEventSteps = BetweenEventSteps,
    Verbose = TRUE
  )

  EventsUnfiltered$Events <- dplyr::left_join(
    EventsUnfiltered$Events %>% dplyr::select(-Success),
    result$Events,
    by = c("Times", "Species", "Environment", "Type")
  )

  # Added Missed Events with NA Success.
  result$Events <- EventsUnfiltered$Events

  save(result,
       file = file.path(dir, paste0(
         "MNA-ExampleExtProp-Result-Env", Environments,
         "-", Space,
         "-", gsub(round(log10(dist)),
                   pattern = "-", replacement = "_", fixed = TRUE),
         "-", EventRateModifiers[1], "-", EventRateModifiers[2],
         "-ExtProp", ExtinctionProportion,
         if(PoolPatchNicheIntervention) "-Intervention",
         ".RData")
       )
  )
  print(Sys.time())

  return(Sys.time() - record)
}

#parallel::stopCluster(clust)
