# Follows TimeSpaceAndTimeSeries-1a-Simulation.R
# Takes in the Simulations, which have documented:
#  abundances (matrix)
#  intervention times and locations (data.frame)
#  intervention types and magnitudes (list?)
#  event time series (data.frame)
#  sampling random seeds (to maintain a high degree of reproducibility)
# and performs sampling immediately prior to, during, and after intervention.
# We then present the usual set of summary statistics in order to verify
# functionality.
# Together, this file and TimeSpaceAndTimeSeries-1a-Simulation.R should produce
# output similar to TimeSpaceAndTimeSeries-1-Bootstraps and *Intervention.R