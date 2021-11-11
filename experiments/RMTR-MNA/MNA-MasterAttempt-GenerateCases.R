# Abstract: ####################################################################
# This file lays out some of the "what" simulations we are looking at doing,
# as well as creating the csv we will reference to initialise each simulation.

library(RMTRCode2)

# Parameters: ##################################################################
systemBase <- list()
systemMods <- list()

### Universal Parameters: ######################################################
# So we are currently thinking 10 histories for each pool-environment combo and
# 10 pool-environment combinations for each parameter set.
systemBase$historiesPerSystem <- 10
systemBase$systemsPerParamSet <- 10

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
  ), ncol = 4),
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
), ncol = 2)

### Spatial Dynamics Parameters: ###############################################
# For the spatial component, we are interested in the following structures.
#   a.   Ring Spatial coupling, distance 1E9.
#   b.   Ring Spatial coupling, distance 1E6.
#   c.   Ring Spatial coupling, distance 1E3.
#   d.   Ring Spatial coupling, distance 1E0.

systemMods$SpaceDistanceMultiplier <- c(1E9, 1E6, 1E3, 1E0)

# Cases: #######################################################################
# Some of these are special cases and are reduced.
#   A = (1, c(IV, VII, VIII, IX, X), c(a, b, c, d)) * 10 * 10          => 2000
#   B = (5, c(I, IV), c(a, b, c, d)) * 10 * 10                         =>  800
# The remainder are grouped together.
#   C = (c(1, 2, 3, 4), c(I, II, III, V, VI), c(a, b, c, d)) * 10 * 10 => 8000
# For a grand total of 10,800 runs.
# As an initial estimate, we saw about 30 minutes per run previously
# and it seems reasonable to get about 30 cores from Viking.
#   10,800 * 0.5 / 30 = 180 hours = 7.5 days.

# We use indices to refer to rows in the above parameters.
cases <- list()

# Case A: ######################################################################
#   Effectively a supplemental figure, but what happens when we combine the
#   changes in rates or we see an extreme "no elimination" case.
#   We don't anticipate that there would be anything *interesting* here, mind,
#   but we do think they are necessary to have touched base with and verified
#   that nothing interesting actually happens.
#   (1, c(IV, VII, VIII, IX, X), c(a, b, c, d))

cases$A <- expand.grid(
  Framework = 1,
  Neutral = c(4, 7:10),
  Space = 1:4
)

# Case B: ######################################################################
#   What happens when we have precise control of the system. No noise, no
#   randomness, just means of distributions doing precisely what we asked.
#   For reference, this also seems like a good "no elimination" case to be
#   aware of so that we know how much the eliminations and noise together
#   drive the ecological dynamics versus how much the framework does.
#   (5, c(I, IV), c(a, b, c, d))

cases$B <- expand.grid(
  Framework = 5,
  Neutral = c(1, 4),
  Space = 1:4
)

# Case C: ######################################################################
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
  Neutral = 1:5,
  Space = 1:4
)

# Make sure you know the cases you are working with...
stopifnot(sum(sapply(cases, nrow)) == 108)

# Assign to each case historiesPerSystem * systemsPerParamSet history seeds,
# and systemsPerParamSet sets of
