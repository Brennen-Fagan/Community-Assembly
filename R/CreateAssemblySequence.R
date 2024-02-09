CreateAssemblySequence <- function(
  Species, NumEnvironments,
  # See https://en.wikipedia.org/wiki/Coupon_collector%27s_problem#Extensions_and_generalizations
  ArrivalEvents = NULL, # Recommend, n = nrow(Pool), n (ln n + 5)
  ArrivalRate = NULL, # If NULL, characteristic time of average interaction matrix.
  ArrivalFUN = NULL, # Takes Events and Rate, Returns Event Times.
  ExtinctEvents = NULL, # Note it is very possible to have misses.
  ExtinctRate = NULL, # If NULL, set to arrival rate.
  ExtinctFUN = NULL, # Takes Events and Rate, Returns Event Times.
  HistorySeed = NULL # Use only one seed, this controls all "external" dynamics.
) {
  ### Set up the History, i.e. arrival and extinction events. ##################
  if (is.null(HistorySeed)) {
    HistorySeed <- 1E8 * runif(1)
  }
  
  if (!is.null(HistorySeed)) {
    if (exists(".Random.seed")) {
      oldSeed <- .Random.seed
    }
    
    set.seed(HistorySeed)
    
    if (is.null(ArrivalEvents)) {
      ArrivalEvents <- 10
    }
    if (is.null(ExtinctEvents)) {
      ExtinctEvents <- ArrivalEvents
    }
    
    if (is.null(ArrivalRate)) { # Roughly
      ArrivalRate <- 1 # / min(abs(Re(eigen(mean(InteractionMatrices)))))
    }
    if (is.null(ExtinctRate)) {
      ExtinctRate <- ArrivalRate
    }
    
    # Retrieve Times.
    if (is.null(ArrivalFUN)) {
      ArrivalFUN <- ArrivalFUN_Example
    }
    if (is.null(ExtinctFUN)) {
      ExtinctFUN <- ExtinctFUN_Example
    }
    Arrivals <- ArrivalFUN(ArrivalEvents, ArrivalRate)
    Extincts <- ExtinctFUN(ExtinctEvents, ExtinctRate)
    
    # Retrieve Species. Retrieve Locations.
    # Note that length(Arrivals) == ArrivalEvents is not guaranteed.
    Events <- data.frame(
      Times = c(Arrivals, Extincts),
      Species = sample.int(Species,
                           size = length(Arrivals) + length(Extincts),
                           replace = TRUE),
      Environment = sample.int(NumEnvironments,
                               size = length(Arrivals) + length(Extincts),
                               replace = TRUE),
      Type = c(rep("Arrival", length(Arrivals)),
               rep("Extinct", length(Extincts))),
      Success = NA
      # Arrival AND Invadable <- TRUE
      # Arrival AND Uninvadable <- FALSE
      # Extinct AND Present     (immediately before extinction) <- TRUE
      # Extinct AND Not Present (immediately before extinction) <- FALSE
    )
    
    # Place in temporal order.
    Events <- with(Events, Events[order(Times), ])
    
    if (exists("oldSeed")) {
      .Random.seed <<- oldSeed
    }
  }
  
  return(list(Events = Events, Seed = HistorySeed))
}