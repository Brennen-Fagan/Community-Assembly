sampleFromResults2 <- function(
  resultAbundance,
  sampling, # Time and Patch to take sample from (with Type for IDing).
  control, # Control Patches
  intervention, # Intervention Time Period
  nSpecies, # Number of species
  samplingPerAbundance, # Convert from Abundance to sample-able individuals.
  samplingFailureRate, # Pr(Researcher Doesn't See Sample)
  PoolTypes # Vector of Species types, assume abundance colnumber is a proxy.
) {
  if(!is.null(PoolTypes)) stopifnot(sum(PoolTypes) == nSpecies ||
                                      length(PoolTypes) == nSpecies)

  # sampling <-
  sampling %>% dplyr::mutate(
    # Descriptions
    Control = ifelse(
      Patch %in% control | Time < intervention,
      "Control", "Experiment"
    )
  ) %>% dplyr::rowwise(
  ) %>% dplyr::mutate(
    # Retrieve values to begin sampling:

    # Location in resultAbundance (note: which.max finds first time after.)
    TimeActualRow = which.max(resultAbundance[, 1] > Time) - 1,
    TimeActual = resultAbundance[TimeActualRow, 1] # First res.Time > Time
  ) %>% dplyr::group_by(
    # dplyr::rowwise doesn't work with group_modify (issue 6870 on github)
    Time, Patch
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      temp <- resultAbundance[
        x$TimeActualRow, # First Col. = Time
        1 + nSpecies * (y$Patch - 1) + 1:nSpecies
        ]

      x %>% dplyr::mutate(
        # Abundances to know number of events:
        SamplingAbundance = sum(temp),

        # Identities and Species Weights:
        SamplingNonZeroSpecies = toString(which(temp > 0)),
        SamplingNonZeroAbundances = toString(nonzero(temp))
      )
    }
  ) %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    SamplingEvents = rpois(n = nrow(sampling),
                           lambda = SamplingAbundance * samplingPerAbundance),
    SamplingObserved = rbinom(n = nrow(sampling),
                              size = SamplingEvents,
                              prob = 1 - samplingFailureRate)
  ) %>% dplyr::group_by(
    Time, Patch
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      if (x$SamplingObserved) {
        draws <- rmultinom(
          1, size = x$SamplingObserved,
          prob = strsplit(x$SamplingNonZeroAbundances, ", ", fixed = TRUE)[[1]]
        )[, 1]

        IDs <- strsplit(x$SamplingNonZeroSpecies, ", ", fixed = TRUE)[[1]]

        if(!is.null(PoolTypes)) {
          if(sum(PoolTypes) == nSpecies)
            drawTypes <- cut(
              as.numeric(IDs), # (ab)Using ID = Row Number.
              c(0, cumsum(PoolTypes)),
              labels = if(!is.null(names(PoolTypes))) names(PoolTypes)
            )
          # Pool$Type[as.numeric(IDs)]
          else if (length(PoolTypes) == nSpecies)
            drawTypes <- PoolTypes[as.numeric(IDs)]
        }

        IDs <- rep(IDs, times = draws)

        if(!is.null(PoolTypes)) {
          if(sum(PoolTypes) == nSpecies)
            Types <- cut(
              as.numeric(IDs), # (ab)Using ID = Row Number.
              c(0, cumsum(PoolTypes)),
              labels = if(!is.null(names(PoolTypes))) names(PoolTypes)
            ) # different since change in IDs from unique to freq. dependent.
          else if (length(PoolTypes) == nSpecies)
            Types <- PoolTypes[as.numeric(IDs)]
        }

        x <- x %>% dplyr::mutate(
          SamplingIDs = toString(IDs),
          SamplingAlpha = sum(draws > 0)
        )
        if (!is.null(PoolTypes)) {
          # Need to make sure each type is presented, even if not in the pop.
          if (!is.table(PoolTypes)) {
            # All Types:
            TypesAsFrame0 <- data.frame(table(PoolTypes))
          } else {
            TypesAsFrame0 <- data.frame(PoolTypes)
            colnames(TypesAsFrame0)[1] <- "PoolTypes"
          }
          TypesAsFrame0 <- TypesAsFrame0 %>% dplyr::mutate(
            PoolTypes = as.character(PoolTypes)
          )
          # Types in Sample:
          TypesAsFrame1 <-
            data.frame(table(drawTypes[draws > 0])) %>% dplyr::mutate(
              Var1 = as.character(Var1)
            )
          TypesAsFrame0 <- dplyr::left_join(
            TypesAsFrame0, TypesAsFrame1,
            by = c("PoolTypes" = "Var1"), suffix = c("Orig", "")
          )
          TypesAsFrame2 <- t(data.frame(TypesAsFrame0$Freq))
          colnames(TypesAsFrame2) <- paste0("SamplingAlphaType",
                                            TypesAsFrame0$PoolTypes)


          x <- dplyr::bind_cols(x %>% dplyr::mutate(
            SamplingTypes = toString(Types)
          ), TypesAsFrame2)
        }
        return(x)
      } else {
        x <- x %>% dplyr::mutate(
          SamplingIDs = "",
          SamplingAlpha = 0
        )
        if (!is.null(PoolTypes)) {
          if (!is.table(PoolTypes)) {
            # All Types:
            TypesAsFrame0 <- data.frame(table(PoolTypes))
          } else {
            TypesAsFrame0 <- data.frame(PoolTypes)
            colnames(TypesAsFrame0)[1] <- "PoolTypes"
          }
          TypesAsFrame0 <- TypesAsFrame0 %>% dplyr::mutate(
            PoolTypes = as.character(PoolTypes),
            Freq = 0
          )
          TypesAsFrame2 <- t(data.frame(TypesAsFrame0$Freq))
          colnames(TypesAsFrame2) <- paste0("SamplingAlphaType",
                                            TypesAsFrame0$PoolTypes)

          x <- dplyr::bind_cols(x %>% dplyr::mutate(
            SamplingTypes = ""
          ), TypesAsFrame2)
        }
        return(x)
      }
    }
  )
}
