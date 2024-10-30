library(dplyr)
library(lhs)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9-Dictionaries.R"))

# Latin Hypercube Sampling Experiments over the parameters
# Jon is particularly keen on k3, but Ines has no strong priors.
# The unfortunate part of this is that parameters were bundled together
# to try to keep the size of the dictionaries smaller.
# Now I need to expand them, and possibly adjust some code.
# On the bright side, we're not looking into interventions right now, because
# the interest is in getting the base number of species up.
# > runif(1)*1e8
# [1] 99139559

# Columns of interest, across three dictionaries:
paramsOfInterest <- c(
  # poolpatchDictionaryOrigin
  "PoolK1InteractionEffectiveness",
  "PoolK2ConsumerSizeAdvantage",
  "PoolK3ConsumerPredationRange",
  "PoolK4ConsumerEfficiency",
  "PoolK5BasalBiomass",
  "PoolK6CoefOfVariation",
  "PoolBasalLogBodySize",
  "PoolConsumerLogBodySize",
  # dynamicsDictionaryOrigin
  "InteractionK1InteractionEffectiveness",
  "InteractionK2ConsumerSizeAdvantage",
  "InteractionK3ConsumerPredationRange",
  "InteractionK4ConsumerEfficiency",
  "InteractionK5BasalBiomass",
  "InteractionK6CoefOfVariation",
  "InteractionEliminationThreshold",
  # eventsDictionaryOrigin
  "ColonizationPropaguleSize"
)

LatinHypercube <- withRandom({
  lhs::optimumLHS(n = 100, # Number of experiments
                  k = length(paramsOfInterest),  # Number of Parameters
                  maxSweeps = 3, eps = 0.1, verbose = FALSE)
}, seed = 99139559)
LatinHypercube <- as.data.frame(LatinHypercube)
colnames(LatinHypercube) <- paramsOfInterest

LatinHypercube <- LatinHypercube %>% dplyr::mutate(
  dplyr::across(.cols = dplyr::starts_with("Pool"), 
                .fns = function(colval) {
                  colname <- dplyr::cur_column()
                  vals <- sort(unique(poolpatchDictionaryOrigin[[colname]]))
                  colval <- ceiling(colval*length(vals)) # unif -> integer.
                  colval <- vals[colval]
                  return(colval)
                })
)
LatinHypercube <- LatinHypercube %>% dplyr::mutate(
  dplyr::across(.cols = dplyr::starts_with("Interaction"), 
                .fns = function(colval) {
                  colname <- dplyr::cur_column()
                  vals <- sort(unique(dynamicsDictionaryOrigin[[colname]]))
                  colval <- ceiling(colval*length(vals)) # unif -> integer.
                  colval <- vals[colval]
                  return(colval)
                })
)
LatinHypercube <- LatinHypercube %>% dplyr::mutate(
  dplyr::across(.cols = dplyr::starts_with("Colonization"), 
                .fns = function(colval) {
                  colname <- dplyr::cur_column()
                  vals <- sort(unique(eventsDictionaryOrigin[[colname]]))
                  colval <- ceiling(colval*length(vals)) # unif -> integer.
                  colval <- vals[colval]
                  return(colval)
                })
)

experiments <- apply(LatinHypercube, MARGIN = 1, function(params) {
  
})
