# Abstract: ####################################################################
# This file lays out some of the "what" simulations we are looking at doing,
# as well as creating the csv we will reference to initialise each simulation.

library(RMTRCode2)

# https://stackoverflow.com/a/15373917
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

# Initialisation: ##############################################################
thisDirectory <- dirname(thisFile())
thisFilePrefix <- "MNA-Master"
thisFileSuffix <- "-Cases.csv"
thisSeparator <- "|"

# > runif(1) * 1E8
# [1] 38427042
if (exists(".Random.seed")) {
  old.seed <- .Random.seed
}
set.seed(38427042)

# Parameters: ##################################################################
systemBase <- list()
systemMods <- list()

### Universal Parameters: ######################################################
# So we are currently thinking 10 histories for each pool-environment combo and
# 10 pool-environment combinations for each parameter set.
systemBase$historiesPerSystem <- 10
systemBase$systemsPerParamSet <- 10
systemBase$environsPerSystem <- 10

# Threshold below which species are removed during calculations.
systemBase$eliminationThreshold <- 10^-4
# Distance above a threshold that a new immigration community starts at.
# Traill et al. 2007 used for inspiration.
systemBase$arrivalDensity <- systemBase$eliminationThreshold * 4 * 10 ^ 3

# Maximum time solver can proceed without checking for elimination.
systemBase$maximumTimeStep <- 1
# Minimum number of steps to reach next event to smooth.
systemBase$betweenEventSteps <- 30

systemBase$species <- c(Basal = 34, Consumer = 66)
systemBase$speciesSpeeds <- 1

# See the coupon collector's problem.
systemBase$eventNumberFunc <- function(Envs, Spec, Const = 5) {
  Envs * ceiling(Spec * (log(Spec + Const)))
}

### Framework Parameters: ######################################################
# (Effectively, the ecological dynamics parameters.)
# We are interested in varying the trophic structure with one of the following.
#   1.   No Change.
#   2.   Another Basal Layer (Below).
#   3.   Another Consumer Layer (Above).
#   4.   Both 2. and 3.
#   5.   No noise.

# In order:
# Interaction Strength, Preferred Prey Size, Willingness to Deviate,
# Efficiency, Basal Equilibrium Biomass Density, and Impact of Noise.
systemBase$LM1996ParamSet <- c(0.01, 10, 0.5, 0.2, 100, 0.1)
# Log body sizes: Basal Min, Basal Max, Consumer Min, Consumer Max
systemBase$LM1996BodySize <- c(-2, -1, -1, 0)

systemMods <- list(
  PoolBodySizes = matrix(c(
     0, 0, 0, 0,
    -1, 0, 0, 0,
     0, 0, 0, 1,
    -1, 0, 0, 1,
     0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  InteractionsNoiseMultiplier = c(1, 1, 1, 1, 0)
)

### Neutral Dynamics Parameters: ###############################################
# We are interested in varying the history structure with one of the following.
#   I.    No Change.
#   II.   1/10 Extinction Rate.
#   III.  10 * Extinction Rate.
#   IV.   No Extinctions.
#   V.    1/10 Arrival Rate.
#   VI.   10 * Arrival Rate.
#   VII.  1/10 Extinction Rate. 1/10 Arrival Rate.
#   VIII. 1/10 Extinction Rate. 10 * Arrival Rate.
#   IX.   10 * Extinction Rate. 1/10 Arrival Rate.
#   X.    10 * Extinction Rate. 10 * Arrival Rate.
# (These are modified from the characteristic rate of the largest perturbation.)

systemMods$NeutralRateMultipliers <- matrix(c(
  # Note reversed from above.
   1,    1,
   1,    0.1,
   1,   10,
   1,    0,
   0.1,  1,
  10,    1,
   0.1,  0.1,
  10,    0.1,
   0.1, 10,
  10,   10
), ncol = 2, byrow = TRUE)

### Spatial Dynamics Parameters: ###############################################
# For the spatial component, we are interested in the following structures.
#   x.   Ring Spatial coupling, distance Inf. (Equivalent to None)
#   a.   Ring Spatial coupling, distance 1E9.
#   b.   Ring Spatial coupling, distance 1E8.
#   c.   Ring Spatial coupling, distance 1E7.
#   d.   Ring Spatial coupling, distance 1E6.
#   e.   Ring Spatial coupling, distance 1E5.
#   f.   Ring Spatial coupling, distance 1E4.
#   g.   Ring Spatial coupling, distance 1E3.
#   h.   Ring Spatial coupling, distance 1E0.

systemMods$SpaceDistanceMultiplier <- 10^c(Inf, 9:3, 0)

# Cases: #######################################################################
# Some of these are special cases and are reduced.
#   A = (1, c(IV, VII, VIII, IX, X), c(x, a:h)) * 10 * 10          =>  4500
#   B = (5, c(I, IV), c(x, a:h)) * 10 * 10                         =>  1800
# The remainder are grouped together.
#   C = (c(1, 2, 3, 4), c(I, II, III, V, VI), c(x, a:h)) * 10 * 10 => 18000
# For a grand total of 21,600 runs.
# As an initial estimate, we saw about 30 minutes per run previously
# and it seems reasonable to get about 30 cores from Viking.
#   24,300 * 0.5 / 30 = 360 hours = 30 days.

# We use indices to refer to rows in the above parameters.
cases <- list()

### Case A: ####################################################################
#   Effectively a supplemental figure, but what happens when we combine the
#   changes in rates or we see an extreme "no elimination" case.
#   We don't anticipate that there would be anything *interesting* here, mind,
#   but we do think they are necessary to have touched base with and verified
#   that nothing interesting actually happens.
#   (1, c(IV, VII, VIII, IX, X), c(a, b, c, d))

cases$A <- expand.grid(
  Framework = 1,
  Neutral = c(4, 7:10),
  Space = 1:9
)

### Case B: ####################################################################
#   What happens when we have precise control of the system. No noise, no
#   randomness, just means of distributions doing precisely what we asked.
#   For reference, this also seems like a good "no elimination" case to be
#   aware of so that we know how much the eliminations and noise together
#   drive the ecological dynamics versus how much the framework does.
#   (5, c(I, IV), c(a, b, c, d))

cases$B <- expand.grid(
  Framework = 5,
  Neutral = c(1, 4),
  Space = 1:9
)

### Case C: ####################################################################
#   This case has all the things that we think are going to be useful and/or
#   interesting from the perspective of the paper.
#   Effectively, this is what happens if we allow our system to change through
#   time (I hesitate to say evolve since I am not wanting to use that
#   terminology) using different parameter values that explore the relationship
#   between the different types of "space" in the system and differing variants
#   of the pool.
#   (To be clear, I do not expect the pool's shape to really change the output,
#   but it would be very important to have that under control if something were
#   to go wrong with it. I would much rather that we look at a different
#   framework, but I agree we should clean this one out to some extent first.)
#   (c(1, 2, 3, 4), c(I, II, III, V, VI), c(a, b, c, d))

cases$C <- expand.grid(
  Framework = 1:4,
  Neutral = c(1:3, 5:6),
  Space = 1:9
)

# Make sure you know the cases you are working with...
stopifnot(sum(sapply(cases, nrow)) == 243)

# Preparations: ################################################################
### Assign Random Seeds: #######################################################
# Assign to each case historiesPerSystem * systemsPerParamSet history seeds,
# and systemsPerParamSet sets of enviromentsPerPool environment seeds.
cases <- lapply(
  cases,
  function(case, base) {
    case$PoolSeeds <- replicate(nrow(case), paste0(runif(
      base$systemsPerParamSet) * 1E8,
      collapse = thisSeparator))
    case$HistorySeeds <- replicate(nrow(case), paste0(runif(
      base$historiesPerSystem * base$systemsPerParamSet) * 1E8,
      collapse = thisSeparator))
    case$EnvironmentSeeds <- replicate(nrow(case), paste0(runif(
      base$environsPerSystem * base$systemsPerParamSet) * 1E8,
      collapse = thisSeparator))
    case
  }, base = systemBase
  )

# Creation: ####################################################################
### Cases: #####################################################################
lapply(
  seq_along(cases),
  function(i, case, name) {
    write.table(
      cases[[i]], file = file.path(
        thisDirectory,
        paste0(thisFilePrefix, name[i], thisFileSuffix)
      ), sep = ","
    )
  }, case = cases,  name = names(cases)
)

### Parameters: ################################################################
save(systemBase, systemMods, thisSeparator, file = file.path(
  thisDirectory,
  paste0(thisFilePrefix, "Parameters", ".RData")
))

# Cleanup: #####################################################################
if (exists("old.seed")) {
  set.seed(old.seed)
}
