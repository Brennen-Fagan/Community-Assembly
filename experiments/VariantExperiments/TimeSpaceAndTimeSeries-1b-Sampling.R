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

# Note that this file should be agnostic to the intervention ideally, outside
# of the information suggested above.

samplingQuantity <- 100 # Not guaranteed!
samplingTimeScaleLogarithmic <- TRUE
# calculationsPlotLong <- FALSE

### Sampling Regime: ##########################################################
# The times (with t = 0 == the intervention time) at which we should sample.
if(samplingTimeScaleLogarithmic) {
  # This version is symmetric on the log scale, centred on 1 time unit,
  # and ends at the time gap. Number of sampling times not guaranteed.
  # The centre is chosen for its relevance to the characteristic time scale.
  samplingTimes <- c(0, unique(exp(c(
    seq(from = log(1),
        to = -log(samplingMaxTime),
        length.out = floor(samplingQuantity/2)),
    seq(from = log(1),
        to = log(samplingMaxTime),
        length.out = ceiling(samplingQuantity/2))
  ))))
} else {
  samplingTimes <- seq(from = 0,
                       by = samplingMaxTime/samplingQuantity,
                       to = samplingMaxTime)
}
