set.seed(55650078)

basal <- c(3, 10, 30, 100, 300, 1000)
consumer <- c(3, 10, 30, 100, 300, 1000) * 2
events <- (max(basal) + max(consumer)) * 2
runs <- 100

logBodySize <- c(-3, 2, -2, 3) # Overlapping body sizes.
parameters <- c(0.01, 10, 0.5, 0.2, 100, 0.1)

# Need to rerun seedsPrep to get the random number generation right for seedsRun
seedsPrep <- runif(2 * length(basal) * length(consumer)) * 1E8
seedsRun <- runif(runs * length(basal) * length(consumer)) * 1E8

Competition = 0
Mutualism = 0 # Suggestion from Coyte 2015, in [0, 1]
CompetitionBasal = 0
Connectance = 1 # in [0, 1]
DiagParam = 0
