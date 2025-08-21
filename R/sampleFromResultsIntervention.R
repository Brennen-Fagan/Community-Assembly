sampleFromResultsIntervention <- function(
  baseAbundance,
  interventionAbundance,
  sampling, # Time and Patch to take sample from (with Type for IDing).
  control, # Control Patches
  interventionTime, # Intervention Time Period
  nSpecies, # Number of species
  samplingPerAbundance, # Convert from Abundance to sample-able individuals.
  samplingFailureRate, # Pr(Researcher Doesn't See Sample)
  Pool
) {
  # ASSUME:
  #  The intervention abundance and the result abundance are the same format
  #  with respect to the number of species (i.e. within environment columns.).

  # sampling <-
  sampling %>% dplyr::mutate(
    # Descriptions
    Control = ifelse(
      Patch %in% control | TimeBase < interventionTime,
      "Control", "Experiment"
    )
  ) %>% dplyr::rowwise(
  ) %>% dplyr::mutate(
    # Retrieve values to begin sampling:

    # Location in resultAbundance (note: which.max finds first time after.)
    TimeActualRow = which.max(baseAbundance[, 1] > TimeBase) - 1,
    TimeActual = baseAbundance[TimeActualRow, 1], # First res.Time > Time
    TimeAlternateRow = which.max(baseAbundance[, 1] > TimeIntervention) - 1,
    TimeAlternate = interventionAbundance[TimeAlternateRow, 1]
  ) %>% dplyr::group_by(
    # dplyr::rowwise doesn't work with group_modify (issue 6870 on github)
    TimeBase, Patch
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      if (x$Control == "Control") {
        temp <- baseAbundance[
          x$TimeActualRow, # First Col. = Time
          1 + nSpecies * (y$Patch - 1) + 1:nSpecies
          ]

        return(x %>% dplyr::mutate(
          # Abundances to know number of events:
          SamplingAbundance = sum(temp),

          # Identities and Species Weights:
          SamplingNonZeroSpecies = toString(which(temp > 0)),
          SamplingNonZeroAbundances = toString(nonzero(temp))
        ))
      } else {
        temp <- interventionAbundance[
          x$TimeAlternateRow, # First Col. = Time
          1 + nSpecies * (x$PatchIntervention - 1) + 1:nSpecies
          ]

        return(x %>% dplyr::mutate(
          # Abundances to know number of events:
          SamplingAbundance = sum(temp),

          # Identities and Species Weights:
          SamplingNonZeroSpecies = toString(which(temp > 0)),
          SamplingNonZeroAbundances = toString(nonzero(temp))
        ))
      }
    }
  ) %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    SamplingEvents = rpois(n = nrow(sampling),
                           lambda = SamplingAbundance * samplingPerAbundance),
    SamplingObserved = rbinom(n = nrow(sampling),
                              size = SamplingEvents,
                              prob = 1 - samplingFailureRate)
  ) %>% dplyr::group_by(
    TimeBase, Patch
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      if (x$SamplingObserved) {
        draws <- rmultinom(
          1, size = x$SamplingObserved,
          prob = strsplit(x$SamplingNonZeroAbundances, ", ", fixed = TRUE)[[1]]
        )[, 1]

        IDs <- strsplit(x$SamplingNonZeroSpecies, ", ", fixed = TRUE)[[1]]
        drawTypes <- Pool$Type[as.numeric(IDs)] # (ab)Using ID = Row Number.
        IDs <- rep(IDs, times = draws)
        Types <- Pool$Type[as.numeric(IDs)]

        x %>% dplyr::mutate(
          SamplingIDs = toString(IDs),
          SamplingTypes = toString(Types),
          SamplingAlpha = sum(draws > 0),
          SamplingAlphaType1 = table(drawTypes[draws > 0])[1],
          SamplingAlphaType2 = table(drawTypes[draws > 0])[2]
        )
      } else {
        x %>% dplyr::mutate(
          SamplingIDs = "",
          SamplingTypes = "",
          SamplingAlpha = 0,
          SamplingAlphaType1 = 0,
          SamplingAlphaType2 = 0
        )
      }
    }
  )
}
