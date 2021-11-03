# This file lays out some of the "what" simulations we are looking at doing,
# as well as creating the csv we will reference to initialise each simulation.

# So we are currently thinking 10 histories for each pool-environment combo and
# 10 pool-environment combinations for each parameter set.

# We are interested in varying the trophic structure with one of the following.
#   1.   No Change.
#   2.   Another Basal Layer (Below).
#   3.   Another Consumer Layer (Above).
#   4.   Both 2. and 3.
#   5.   Strictly Basal species.

# We are interested in varying the history structure with one of the following.
#   I.   No Change.
#   II.  1/10 Extinction Rate.
#   III. 10 * Extinction Rate.
#   IV.  No Extinctions.
#   V.   1/10 Arrival Rate.
#   VI.  10 * Arrival Rate.

# For the spatial component, we are interested in the following structures.
#   a.   No Spatial coupling, distance Inf. (Eq. to Ring/2-Ring, dist. Inf.)
#   b.   Ring Spatial coupling, distance 1E9.
#   c.   Ring Spatial coupling, distance 1E6.
#   d.   Ring Spatial coupling, distance 1E3.
#   e.   Ring Spatial coupling, distance 1E0.
#   f.   2-Ring Spatial coupling, distance 1E9.
#   g.   2-Ring Spatial coupling, distance 1E6.
#   h.   2-Ring Spatial coupling, distance 1E3.
#   i.   2-Ring Spatial coupling, distance 1E0.

# Finally, we might also consider the Anthropocene impact.
# For this, we will be varying parameters midway through the experiment.
# We do not change the pool or environment; 1 - 5 stay the same.
# The options that we consider then are:
#   Any of I - VI to any of I - VI.
#   Any of a - e to any of a - e.
#   Any of a, f - i to any of a, f - i as well as any of a, f - i with the rings
#   changed from 1 - 5, 6 - 10 to 3 - 7, 8 - 2.

# As can be seen, this is a lot of combinations, roughly:
#   10 * 10 * 5 * (6 * 6 * 9 + 6 * 5 * 5 + 6 * 5 * 10) = 387,000 runs.

# We can reduce this somewhat if we acknowledge that runs with no treatment
# are equivalent between the three alternation cases above.
#   10 * 10 * 5 * (6 * 5 * 9 + 6 * 5 * 4 + 6 * 5 * 9 + 6 * 9) = 357,000 runs.

# There is the potential to reduce this further by considering the treatment
# runs as stemming from the same base system. In this case, we would perform
#   10 * 10 * 5 * 6 * 9 = 27,000 "half-runs".
# Then each run would be finished with its original conditions as
#   -> 27,000 "half-runs"
# as well as being used with treatments:
#   27,000 * 5 "half-runs" for rate variation
#   10 * 10 * 5 * 6 * (1 * 8 + 4 * 4 + 4 * (4 + 4)) = 168,000 "half-runs"
#     for space variations.
# This results again in 357,000 as expected, but half-runs instead of full runs.
# This is more coding intensive, but could save us a significant amount of time.
# The question is merely of trade-offs.
