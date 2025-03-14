
# Functions: ##################################################################

### General: ##################################################################
nonzero <- function(vec) {
  vec[vec!=0]
}

# Need to use stackoverflow.com/a/47012149 to convert our arguments to a list.
callMeMaybe2 <- function(listofcharargs) {
  Charargs = unlist(listofcharargs)
  if(is.null(Charargs)) return(alist())
  eval(parse(
    text = paste0("alist(",
                  paste(parse(text = Charargs),
                        collapse = ","),")")
  ))
}

# We need a function to generate control random values.
# https://stackoverflow.com/a/59875367 Gwang-Jin Kim and
# https://stackoverflow.com/a/14324316 Romain Francois (same question).
withRandom <- function(expr, seed) {
  if (exists(".Random.seed")) {
    oldSeed <- .Random.seed
    on.exit({.Random.seed <<- oldSeed})
  }
  set.seed(seed)
  expr
}

# Function that automatically increments after use. Note that we need to store
# the result if we want to use the same index twice.
indexFactory <- function() {
  index <- 1
  function() {
    on.exit(index <<- index + 1)
    return(index)
  }
}

### Randomness: ###############################################################
# Run runif and organise in a smooth-ish ring.
runifRing <- function(n, ...) {
  indices <- if (n %% 2) {
    # Odd (why?)
    c(1, seq(from = 2, by = 2, to = n), seq(from = n, by = -2, to = 2))
  } else {
    # Even.
    c(1, seq(from = 2, by = 2, to = n), seq(from = n - 1, by = -2, to = 2))
  }
  sort(runif(n, ...))[indices]
}

# Discrete niche samplers.
sample.int.normalized <- function(n, slots = 2) {
  (sample.int(slots, size = n, replace = TRUE) - 1) / (slots - 1)
}
sample.int.3 <- purrr::partial(sample.int.normalized, slots = 3)

### Parameters: ###############################################################
# reproductive rate r' = r * rho ^ (sign(r)), where rho is a function of
# distance between patch and species in some sort of trait way.
rhofunction <- function(
  base = 2, offset = 0, multiplier = 1, metric = "euclidean"
) {
  force(base);force(offset);force(multiplier)
  function(m, n) {
    base ^ (offset - multiplier * dist(
      matrix(c(m, n), byrow = TRUE, nrow = 2), method = metric)
    )
  }
}

rho.noop <- function(m, n) {1}
rho.2.0.1.euclidean <- rhofunction()
rho.2.1.2.euclidean <- rhofunction(2, 1, 2)
rho.5.0.1.euclidean <- rhofunction(5, 0, 1)
rho.5.1.2.euclidean <- rhofunction(5, 1, 2)
rho.10.0.1.euclidean <- rhofunction(10, 0, 1)
rho.10.1.2.euclidean <- rhofunction(10, 1, 2)

# Easy ring gradients.
patchTypes.0.Half.1 <- function(n) {
  toString(
    c(0,
      rep(0, floor((n - 2)/3)),
      rep(0.5, ceiling((n - 2)/6)),
      1,
      rep(1, floor((n - 2)/3)),
      rep(0.5, ceiling((n - 2)/6))
    )
  )
}

# Why? so we can have single argument functions with partials.
repFixed <- function(value = 0.5) {
  force(value)
  function(n) {rep(value, n)}
}
rep_0 <- repFixed(0)
rep_0.25 <- repFixed(0.25)
rep_0.5 <- repFixed()
rep_0.75 <- repFixed(0.75)
rep_1 <- repFixed(1)

evensplit <- function(values = c(0, 1)) {
  force(values)
  function(n) {
    c(rep(values, times = floor(n / length(values))),
      if (n %% length(values) != 0) {
        values[1:(n %% length(values))]
      })
  }
}
evensplit_01 <- evensplit()
evensplit_0.51 <- evensplit(c(0.5, 1))

gradientline_01 <- function(n) {
  c(rep(0, ceiling(n/2)), rep(1, floor(n/2)))
}
gradientline_0half1 <- function(n) {
  left <- rep(0, floor(n / 3)); right <- rep(1, floor(n/3))
  return(c(left, rep(0.5, n - length(left) - length(right)), right))
}

# Coupon Collector's Problem
# I think this is probably higher accuracy than the previous version.
defaultEvents <- function(
  NumberOfEnvironments, NumberOfSpecies, constant = 3
) {
  ceiling(
    NumberOfEnvironments * NumberOfSpecies * (
      log(NumberOfEnvironments * NumberOfSpecies) + constant
    )
  )
}

### Loading: ##################################################################
loadSimulation <- function(filepath1, filepath2 = NULL) {
  file1 <- load(filepath1)
  if (length(file1) > 1) {
    warning("filepath1 has more than one object. Defaulting to first.")
  }
  file1 <- get(file1[1])

  EntriesRequired <- c("Events", "Abundance")
  if (any(! EntriesRequired %in% names(file1))) {
    error("Events or Abundance not found in filepath1 object 1.")
  }

  if (!is.null(filepath2)) {
    file2 <- load(filepath2)
    if (length(file2) > 1) {
      warning("filepath2 has more than one object. Defaulting to first.")
    }
    file2 <- get(file2[1])
    if (any(! EntriesRequired %in% names(file2))) {
      error("Events or Abundance not found in filepath2 object 1.")
    }

    EntriesToCheck <- !names(file2) %in% c("Events", "Abundance")

    stopifnot(isTRUE(all.equal(file1[EntriesToCheck],
                               file2[EntriesToCheck])))

    file1$Abundance <- file1$Abundance[
      file1$Abundance[, 1] < min(file2$Abundance[, 1]),
      ]

    file1$Events <- rbind(file1$Events, file2$Events)
    file1$Events <- file1$Events %>% dplyr::distinct()
    file1$Abundance <- rbind(file1$Abundance, file2$Abundance)
  }

  return(file1)
}

retrieveFunction <- function(funcstring) {
  funcstring <- strsplit(funcstring, split = "::")
  if (length(funcstring) > 1) {
    stop(paste0("Too many functions provided in string: ", length(funcstring)))
  } else {
    funcstring <- funcstring[[1]]
  }
  if (length(funcstring) > 2) {
    stop(paste0("Too many parts to function provided: ",
                length(funcstring)))
  } else if (length(funcstring) == 2) {
    funcstring <- get(funcstring[2], envir = loadNamespace(funcstring[1]))
  } else if (length(funcstring) == 1) {
    funcstring <- get(funcstring[1])
  } else {
    stop("No parts found for function.")
  }
  return(funcstring)
}

### Conversion: ###############################################################
convertDispersalDictToDistMatrix <- function(dispersalDictionary, nEnv) {
  with(c(
    dispersalDictionary,
    Environments = nEnv
  ), {
    if (Configuration == "None") {
      DistanceMatrix <- Matrix::sparseMatrix(
        i = Environments, j = Environments, x = 0)
    }
    if (Configuration == "Ring" || Configuration == "Line")
      DistanceMatrix <- Matrix::bandSparse(
        Environments, k = c(-1, 1),
        diagonals = list(rep(Resistance, Environments - 1),
                         rep(Resistance, Environments - 1))
      )
    if (Configuration == "Ring") {
      DistanceMatrix[Environments, 1] <- Resistance
      DistanceMatrix[1, Environments] <- Resistance
    }
    if (Configuration == "Complete") {
      DistanceMatrix <- matrix(Resistance,
                               nrow = Environments,
                               ncol = Environments)
      diag(DistanceMatrix) <- 0
    }
    return(DistanceMatrix)
  }
  )
}

### Sampling: #################################################################
sampleFromResults <- function(
  resultAbundance,
  sampling, # Time and Patch to take sample from (with Type for IDing).
  control, # Control Patches
  intervention, # Intervention Time Period
  nSpecies, # Number of species
  samplingPerAbundance, # Convert from Abundance to sample-able individuals.
  samplingFailureRate # Pr(Researcher Doesn't See Sample)
) {
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
        1 + nSpecies * (y$Patch - 1) + 1:nSpecies]

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
        drawTypes <- Pool$Type[as.numeric(IDs)] # (ab)Using ID = Row Number.
        IDs <- rep(IDs, times = draws)
        Types <- Pool$Type[as.numeric(IDs)] # NOTE: Bug here! Pool not defined.

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

### Computations: #############################################################
# Okay, next, one common topic is essentially native vs invasive.
# In truth all species are more or less the same, originating from the true
# regional pool. Our researchers wouldn't know this.
# Instead, they would both use their control to estimate the natives.
computeSpeciesInControl <- function(sampling,
                                    Time = "Time",
                                    IDColumn = "SamplingIDs",
                                    OutPrefix = "Sampling") {
  IDColNum <- which(colnames(sampling) == IDColumn)

  controlSpecies <- sampling %>% dplyr::filter(
    Control == "Control"
  ) %>% dplyr::pull(
    IDColNum
  ) %>% strsplit(
    ", ", fixed = TRUE
  ) %>% unlist(
  ) %>% unique()

  sampling %>% dplyr::group_by(
    dplyr::across(dplyr::all_of(Time)), Patch
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      splitIDs <- strsplit(x[[IDColumn]], ", ", fixed = T)[[1]]
      iDsInControl <- splitIDs %in% controlSpecies

      # Note: group_by (Patch) enforces Alpha / Local scale
      x <- x %>% dplyr::mutate(
        IDsNative = toString(splitIDs[iDsInControl]),
        AbundanceNative = sum(iDsInControl),
        AlphaNative = length(unique(splitIDs[iDsInControl])),
        IDsInvasive = toString(splitIDs[!iDsInControl]),
        AbundanceInvasive = sum(!iDsInControl),
        AlphaInvasive = length(unique(splitIDs[!iDsInControl]))
      )

      colnames(x)[(ncol(x)-5):ncol(x)] <-
        paste0(OutPrefix, colnames(x)[(ncol(x)-5):ncol(x)])

      return(x)
    }
  )
}

addDistanceColumns <- function(bootstrapSamples, mindig = 1, Time = "Time") {
  bootstrapSamples %>% dplyr::mutate(
    ##### Temporal Distance: ##################################################
    TimeSinceStart = round(# Rare 1 != 1 issue.
      TimeSinceStart,
      digits = mindig
    )
  ) %>% dplyr::group_by(
    Type, Control, Bootstrap, Patch
  ) %>% dplyr::arrange(
    TimeSinceStart
  ) %>% dplyr::mutate(
    TimeGapNumber = seq_along(TimeSinceStart),
    TimeGapNumber = ifelse(
      Type == "Time series" & Control == "Control",
      rev(TimeGapNumber), TimeGapNumber
    ) # i.e.  5, 4, 3, 2, 1, ___, 1, 2, 3, 4, 5
  ) %>% dplyr::ungroup(
    ##### Spatial Distance: #################################################
  ) %>% dplyr::group_by(
    Type, Control, Bootstrap, TimeSinceStart
  ) %>% dplyr::mutate(
    DistanceFromCenter = dplyr::case_when(
      min(Patch) == 1 & max(Patch) == 10 ~ {
        ifelse(Patch < 5, Patch + 10, Patch) - median(
          ifelse(Patch < 5, Patch + 10, Patch)
        )
      },
      TRUE ~ Patch - median(Patch)
    ),
    DistanceFromCenterExpRev = ifelse(Control == "Experiment" &
                                        Type == "Space for time",
                                      -DistanceFromCenter,
                                      DistanceFromCenter)
  ) %>% dplyr::ungroup(
    ##### Define control species: ###########################################
  ) %>% dplyr::group_by(
    DistanceFromCenterExpRev, Bootstrap, TimeGapNumber, Type
  ) %>% dplyr::group_modify(
    .f = ~ computeSpeciesInControl(.x, Time = Time)
  ) %>% dplyr::ungroup(
  )
}

# NOTE: Breaks idiom; should be renamed to be a verb.
SpeciesStringsToBeta <- function(
  .x, .y,
  SpeciesColumn = "SamplingNonZeroSpecies",
  AbundanceColumn = "SamplingNonZeroAbundances",
  Method = "Jaccard",
  PresenceAbsence = TRUE
) {
  stopifnot(length(.x$Patch) == length(unique(.x$Patch)))

  .x$Patch <- as.character(.x$Patch)

  flagMeaningless <- .x$Patch[.x[[AbundanceColumn]] == ""]

  uniqueSpecies <- sort(unique(unlist(strsplit(
    .x[[SpeciesColumn]],
    split = ", "))))

  comdatmat <- matrix(0,
                      nrow = length(.x$Patch),
                      ncol = length(uniqueSpecies))
  colnames(comdatmat) <- uniqueSpecies
  rownames(comdatmat) <- .x$Patch

  for(i in seq_along(.x$Patch)) {
    # print(strsplit(.x$SamplingNonZeroSpecies[i], split = ", ")[[1]])
    # print(strsplit(.x$SamplingNonZeroAbundances[i], split = ", ")[[1]])
    if (SpeciesColumn == AbundanceColumn) {
      # Abundance = Species => List of IDs with repetition for abundance.
      vals <- table(strsplit(.x[[AbundanceColumn]][i], split = ", ")[[1]])
      indices <- names(vals)
    } else {
      # Else, abundance in same order as species detected.
      vals <- as.numeric(strsplit(.x[[AbundanceColumn]][i], split = ", ")[[1]])
      indices <- strsplit(.x[[SpeciesColumn]][i], split = ", ")[[1]]
    }

    comdatmat[.x$Patch[i], indices] <- vals

  }

  if (PresenceAbsence) comdatmat <- comdatmat > 0

  # Short for Jaccard, our default.
  Jacs <- vegan::vegdist(method = tolower(Method), x = comdatmat)

  data.frame(
    Beta = Jacs[1:length(Jacs)],
    Patch1 = rep(attr(Jacs, "Labels"), (length(attr(Jacs, "Labels"))-1):0),
    Patch2 = attr(Jacs, "Labels")[
      sequence(from = seq_along(attr(Jacs, "Labels"))[-1],
               nvec = (length(attr(Jacs, "Labels")) - 1):1)
      ]
  ) %>% dplyr::mutate(
    Meaningless = Patch1 %in% flagMeaningless | Patch2 %in% flagMeaningless
  )

}


ConvertPreparedToBeta <- function(
  prepared, columns, method, presenceabsence, postfixes, indicator
) {
  # Lazily didn't implement. Requires editing combine step.
  stopifnot(length(columns) == 3)
  # Lazily, also didn't remove the dependence on Bootstrap and Type.
  # This should probably be a complement.

  # # Prepared looks like, with indicator as the first grouping variable.
  # bootstrapSamplesPairedBeta_source %>% dplyr::filter(
  #   TimeGapNumber == preferredTimeGap
  # ) %>% dplyr::group_by(
  #   DistanceFromCenterExpRev, Bootstrap, Type
  # )
  # # or like
  # bootstrapSamplesPairedBeta_source %>% dplyr::filter(
  #   DistanceFromCenterExpRev == 0
  # ) %>% dplyr::group_by(
  #   TimeGapNumber, Bootstrap, Type
  # )

  # Generate the Beta Diversities Across the Columns
  retval <- lapply(
    columns,
    function(Column, dat, method) {
      temp <- dat %>% dplyr::group_modify(
        .f = SpeciesStringsToBeta,
        SpeciesColumn = Column,
        AbundanceColumn = Column,
        PresenceAbsence = presenceabsence,
        Method = method
      )
      colnames(temp)[colnames(temp) == "Beta"] <- method
      return(temp)
    },
    dat = prepared,
    method = method
  )

  # Rename in preparation for combining.
  retval <- lapply(
    seq_along(retval),
    function(i, dat, postfixes) {
      indices <- which(colnames(dat[[i]]) %in% c(method, "Meaningless"))
      colnames(dat[[i]])[indices] <-
        paste0(colnames(dat[[i]])[indices], postfixes[i])
      dat[[i]]
    },
    dat = retval, postfixes = postfixes
  )

  # Combine
  retval <- dplyr::full_join(
    retval[[1]],
    retval[[2]],
    by = c(dplyr::all_of(indicator), "Bootstrap", "Type", "Patch1", "Patch2")
  ) %>% dplyr::full_join(
    retval[[3]],
    by = c(dplyr::all_of(indicator), "Bootstrap", "Type", "Patch1", "Patch2")
  ) %>% tidyr::pivot_longer(
    cols = -c(dplyr::all_of(indicator), "Bootstrap", "Type", "Patch1", "Patch2"),
    names_to = c("Measure", "Subset"),
    names_sep = ", "
  ) %>% tidyr::pivot_wider(
    names_from = "Measure", values_from = "value"
  )

  return(retval)
}


### Iterate ThinAndCalculate for non-binary: ##################################
thinAbundance <- function(abundance, events, threshold,
                          preferred_rows_per_event) {
  bythin <- floor((nrow(abundance)
                   / nrow(events))
                  / preferred_rows_per_event)

  abundance <- abundance[seq(from = 1,
                             to = nrow(abundance),
                             by = bythin), ]

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- abundance[, -1] < threshold # & abundance[, -1] > 0
  abundance[, -1][toEliminate] <- 0

  return(abundance)
}

thinAbundanceEqualTimeSteps <- function(
  abundance, threshold, preferredTimeStep, preferredStart = NULL
) {
  time <- abundance[, 1]
  consistentDistance <- max(diff(time), preferredTimeStep)
  if (consistentDistance != preferredTimeStep) {
    warning(paste0("preferredTimeStep is too small, minimum:",
                   consistentDistance))
  }
  if (is.null(preferredStart)) preferredStart <- min(time)
  targets <- seq(from = preferredStart,
                 to = max(time), by = consistentDistance)
  rows <- unique( # Just In Case?
    sapply(targets, function(x, y) {which.max(y >= x)}, y = time)
  )

  abundance <- abundance[rows, ]
  abundance[, 1] <- targets
  # Technically an approximation, but we should be high resolution
  # enough for it not to be a problem.

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- abundance[, -1] < threshold # & abundance[, -1] > 0
  abundance[, -1][toEliminate] <- 0

  return(abundance)
}

hillWrapper <- function(env, i) {
  if (length(dim(env)) >= 1 && ncol(env) == 1) {
    metrics <- do.call(rbind, lapply(env, function(i)
      if(i > 0) {
        c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
      } else {
        c(NA, NaN, NaN, NA, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
      }
    ))
    metrics <- data.frame(metrics)
    colnames(metrics) <-  c("0", "0.25", "0.5", "1", "2", "4",
                            "8", "16", "32", "64", "Inf")
  } else {
    metrics <- vegan::renyi(env, hill = TRUE) # returns data.frame, but only if
    if (length(dim(env)) <= 1 || nrow(env) == 1) {# nrow > 1.
      tempnames <- names(metrics)
      metrics <- data.frame(matrix(metrics, nrow = 1))
      colnames(metrics) <- tempnames
    }
  }
  names(metrics) <- paste0("Hill:", names(metrics))
  cbind(Environment1 = i,
        Environment2 = NA,
        metrics,
        stringsAsFactors = FALSE) %>% tidyr::pivot_longer(
          cols = tidyr::contains("Hill"),
          names_to = "Metric", values_to = "Value"
        )
}

# Based on RMTRCode2::Calculate_Diversity and calculateAbundanceMetrics as well
# as the betapart procedure.
# The main difference is to try to provide a uniform layout from the beginning
calculateDiversityMetrics <- function(abundance, nspecies, nenvironments) {
  stopifnot(nenvironments >= 1)

  envs <- lapply(
    1:nenvironments,
    function(i, abund, numSpecies) {
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      env <- matrix(env, nrow = nrow(abund))
      return(env)
    },
    abund = abundance,
    numSpecies = sum(nspecies)
  )

  time <- abundance[, 1]

  # Format: data.frame: Time, Environment 1, Environment 2, Metric, Value
  # If an environment is not appropriate, we instead use NA.

  # Alpha: ####################################################################
  diversityAlpha <- dplyr::bind_rows(lapply(
    1:nenvironments,
    function(i) {
      vals <- hillWrapper(env = envs[[i]], i = i)
      cbind(Time = rep(time, each = nrow(vals)/length(time)), vals)
    }
  )) %>% dplyr::mutate(
    Metric = paste0("Alpha ", Metric)
  )

  # Gamma: ####################################################################
  if (nenvironments > 1) {
    envgamma <- envs[[1]]
    if (nenvironments > 1) {
      for (i in 2:nenvironments) {
        envgamma <- envgamma + envs[[i]]
      }
    }

    diversityGamma <-
      hillWrapper(env = envgamma, i = NA) %>% dplyr::mutate(
        Metric = paste0("Gamma ", Metric)
      )
    diversityGamma <- cbind(
      Time = rep(time, each = nrow(diversityGamma)/length(time)),
      diversityGamma
    )
  }

  # Gamma Temporal: ###########################################################
  diversityGammaTemporal <- dplyr::bind_rows(lapply(
    1:nenvironments,
    function(i) {
      # Gamma (spatial) adds across environments; this adds across times.
      vals <- hillWrapper(env = colSums(envs[[i]]), i = i)
      cbind(Time = NA, vals)
    }
  )) %>% dplyr::mutate(
    Metric = paste0("GammaTemporal ", Metric)
  )

  # Beta Spatial: #############################################################
  if (nenvironments > 1) {
    diversityBetaSpace <- apply(
      abundance,
      MARGIN = 1, # Rows
      function(row, nenv) {
        thistime <- row[1]

        # List with three components:
        #   balance/turnover, gradient/nestedness, and total.
        distsBC <- betapart::beta.pair.abund(
          x = matrix(row[-1], nrow = nenv, byrow = TRUE),
          index.family = "bray"
        )
        distsJ <- betapart::beta.pair(
          x = matrix((row[-1] > 0) + 0, nrow = nenv, byrow = TRUE),
          index.family = "jaccard"
        )

        dataf <- expand.grid(
          Environment1 = 1:nenv,
          Environment2 = 1:nenv
        ) %>% dplyr::filter(
          Environment1 < Environment2
        ) %>% dplyr::arrange(
          Environment1, Environment2
        ) %>% dplyr::mutate(
          Time = thistime,
          SpaceBrayCurtisBalance = as.vector(distsBC$beta.bray.bal),
          SpaceBrayCurtisGradient = as.vector(distsBC$beta.bray.gra),
          SpaceBrayCurtis = as.vector(distsBC$beta.bray),
          SpaceJaccardTurnover = as.vector(distsJ$beta.jtu),
          SpaceJaccardNestedness = as.vector(distsJ$beta.jne),
          SpaceJaccard = as.vector(distsJ$beta.jac)
        )

        return(dataf)
      },
      nenv = nenvironments
    ) %>% dplyr::bind_rows() %>% tidyr::pivot_longer(
      cols = SpaceBrayCurtisBalance:SpaceJaccard,
      names_to = "Metric", values_to = "Value"
    )
  }

  # Beta Temporal: ############################################################
  # ASSUMING ALREADY EQUAL TIME DIFFERENCES.

  diversityBetaTime <- dplyr::bind_rows(lapply(
    1:nenvironments, function(i) {
      lapply(2:nrow(envs[[i]]), function(r) {
        target <- rbind(envs[[i]][r-1, ], envs[[i]][r, ])
        distsBC <- betapart::beta.pair.abund(
          x = target, index.family = "bray")
        distsJ <- betapart::beta.pair(
          x = (target > 0) + 0, index.family = "jaccard")
        expand.grid(
          Environment1 = i,
          Environment2 = NA
        ) %>% dplyr::mutate(
          Time = time[r],
          TimeBrayCurtisBalance = as.vector(distsBC$beta.bray.bal),
          TimeBrayCurtisGradient = as.vector(distsBC$beta.bray.gra),
          TimeBrayCurtis = as.vector(distsBC$beta.bray),
          TimeJaccardTurnover = as.vector(distsJ$beta.jtu),
          TimeJaccardNestedness = as.vector(distsJ$beta.jne),
          TimeJaccard = as.vector(distsJ$beta.jac)
        )
      }) %>% dplyr::bind_rows()
    }
  )) %>% tidyr::pivot_longer(
    cols = TimeBrayCurtisBalance:TimeJaccard,
    names_to = "Metric", values_to = "Value"
  )

  return(dplyr::bind_rows(
    diversityAlpha,
    if (nenvironments > 1) diversityGamma,
    if (nenvironments > 1) diversityBetaSpace,
    diversityBetaTime
  ))
}

# Like Calculating Diversities, but not for binary data only.
calculateAbundanceMetrics <- function(abundance, nspecies, nenvironments) {
  stopifnot(nenvironments >= 1)

  envs <- lapply(
    1:nenvironments,
    function(i, abund, numSpecies) {
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
    },
    abund = abundance,
    numSpecies = sum(nspecies)
  )

  time <- abundance[, 1]

  ### Alpha Diversity: ##################################################
  diversity_alpha <- lapply(
    1:nenvironments,
    function(i, abund, numSpecies) {
      env <- envs[[i]]
      env_basal <- env[, 1:nspecies[1], drop = FALSE]
      env_consumer <- env[, nspecies[1] + 1:nspecies[2], drop = FALSE]
      # fix a bug in vegan for single species ecosystems...
      if(ncol(env_consumer) == 1) env_consumer <- cbind(env_consumer, 0)
      metrics <- vegan::renyi(env, hill = TRUE)
      metrics_basal <- vegan::renyi(env_basal, hill = TRUE)
      metrics_consumer <- vegan::renyi(env_consumer, hill = TRUE)
      names(metrics) <- paste0(names(metrics), ", All")
      names(metrics_basal) <- paste0(names(metrics_basal), ", Basal")
      names(metrics_consumer) <- paste0(names(metrics_consumer), ", Consumer")
      cbind(Time = time,
            metrics, metrics_basal, metrics_consumer,
            Environment = i,
            stringsAsFactors = FALSE)
    },
    abund = abundance,
    numSpecies = sum(nspecies)
  )

  diversity_alpha <- dplyr::bind_rows(diversity_alpha)

  print("alpha")
  ### Gamma Diversity: ##################################################
  # In contrast, need to combine across environments.
  envgamma <- envs[[1]]
  if (nenvironments > 1) {
    for (i in 2:nenvironments) {
      envgamma <- envgamma + envs[[i]]
    }
  }

  env_basal <- envgamma[, 1:nspecies[1], drop = FALSE]
  env_consumer <- envgamma[, nspecies[1] + 1:nspecies[2], drop = FALSE]
  # fix a bug in vegan for single species ecosystems...
  if(ncol(env_consumer) == 1) env_consumer <- cbind(env_consumer, 0)

  metrics <- vegan::renyi(envgamma, hill = TRUE)
  metrics_basal <- vegan::renyi(env_basal, hill = TRUE)
  metrics_consumer <- vegan::renyi(env_consumer, hill = TRUE)
  names(metrics) <- paste0(names(metrics), ", All")
  names(metrics_basal) <- paste0(names(metrics_basal), ", Basal")
  names(metrics_consumer) <- paste0(names(metrics_consumer), ", Consumer")
  diversity_gamma <- cbind(Time = time,
                           metrics, metrics_basal, metrics_consumer,
                           Environment = NA,
                           stringsAsFactors = FALSE)

  print("gamma")
  ### Beta Diversity (Jaccard, Space): ##################################
  diversity_beta <- apply(
    abundance,
    MARGIN = 1, # Rows
    function(row, envs) {
      thistime <- row[1]

      # Vegan complains about rows with all 0's.
      # The warning is generic, so we cannot silence it specifically.
      dists <- suppressWarnings(vegan::vegdist(
        method = "bray",
        x = matrix(row[-1], nrow = envs, byrow = TRUE)
      ))

      dataf <- expand.grid(
        Env1 = 1:envs,
        Env2 = 1:envs
      ) %>% dplyr::filter(
        Env1 < Env2
      ) %>% dplyr::arrange(
        Env1, Env2
      ) %>% dplyr::mutate(
        Time = thistime,
        BrayCurtis = as.vector(dists)
      )

      return(dataf)
    },
    envs = nenvironments
  ) %>% dplyr::bind_rows()

  print("beta")
  ### Return Diversities: ###############################################
  return(list(
    alpha = diversity_alpha,
    beta = diversity_beta,
    gamma = diversity_gamma
  ))
}

calculateColExtMetrics <- function(sim) {
  # Might be easier to take the binary matrix and to say
  # 0->1 and in the successful events ~ colonization
  # 0->1 and not in the successful events ~ dispersal
  # 1->0 and in the successful events ~ neutral extirpation
  # 1->0 and not in the successful events ~ dynamic extirpation
  # But note that more than one of these can co-occur.
  binaryMatrix <-
    sim$Abundance[, -1] > sim$Parameters$EliminationThreshold
  ColExtMatrix <- # Next - Current so timing of event lines up with the change
    binaryMatrix[2:nrow(sim$Abundance), ] -
    binaryMatrix[1:(nrow(sim$Abundance) - 1), ]

  nspecies <- ncol(ColExtMatrix)

  # Reattach the times
  ColExtMatrix <- cbind(sim$Abundance[1:(nrow(sim$Abundance) - 1), 1],
                        ColExtMatrix)

  allEvents <- lapply(1:nrow(ColExtMatrix), function(i, mat, binmat, events) {
    time <- mat[i, 1]
    changes <- which(mat[i, -1] != 0)
    event <- events %>% dplyr::filter(Times == time, Success) # N.B. nrow can >1

    if (nrow(event) == 0 && length(changes) == 0) {
      # Done, nothing to report.
      return(NULL)
    }

    environments <- ((changes - 1) %/% nspecies) + 1 # 1:20 w/10 -> 10 1s, 10 2s
    species <- ((changes - 1) %% nspecies) + 1 # 1:20 w/10 -> c(1:10, 1:10)
    stopifnot(environments <= sim$NumEnvironments)

    if (nrow(event) == 1 && length(changes) == 1 &&
        event$Species == species && event$Environment == environments &&
        ((mat[i, changes + 1] == 1 && event$Type == "Arrival") ||
         (mat[i, changes + 1] == -1 && event$Type == "Extinct")
        )) {
      # Trivial, everything is already recorded
      return(event)
    }

    # Bundle information together for easier processing.
    changesdf <- data.frame(
      change = changes,
      type = mat[i, changes+1], # +1 column for time
      species = species,
      environments = environments
    ) %>% dplyr::mutate(index = dplyr::row_number())

    # Else: complicated, some combination of events or non-neutral events.
    # Proceed along events, removing any changes that are already accounted for.
    # (Essentially, we're doing a complicated implementation of a filtering join
    # where we're trying to filter out things that have already been addressed.)
    # Simultaneously, if an event occurs and the result is not registered as a
    # change, then something must have countered it. If neutral events cannot
    # occur simultaneously, then it must be a non-neutral event.
    if (nrow(event) > 0) {
      for (j in 1:nrow(event)) {
        thisEvent <- event[j, ]
        thisChange <- changesdf %>% dplyr::filter(
          species == thisEvent$Species,
          environments == thisEvent$Environment
        )

        if (nrow(thisChange) == 1 &&
            ((thisChange$type == 1 && thisEvent$Type == "Arrival") ||
             (thisChange$type == -1 && thisEvent$Type == "Extinct")
            )) {
          # This change is already in event(s) as this event.
          # Remove it from the data.frame of event(s) to add.
          changesdf <- changesdf %>% dplyr::filter(index != thisChange$index)
        } else if (
          nrow(thisChange) == 0 # then an event happened, but no change recorded.
          # result depends on whether a population is there or not.
        ) {
          speciesPresent <- with(
            thisEvent,
            binmat[i, Species + (Environment - 1) * nspecies] # No Time Col.
          )
          newEvent <- data.frame(
            Times = time,
            Species = thisEvent$Species,
            Environment = thisEvent$Environment,
            Type = if(thisEvent$Type == "Arrival" && speciesPresent) {
              "Present"
            } else if (thisEvent$Type == "Arrival" && !speciesPresent) {
              "Dynamic Loss"
            } else if (thisEvent$Type == "Extinct" && speciesPresent) {
              "Dispersal"
            } else if (thisEvent$Type == "Extinct" && !speciesPresent) {
              "Dispersal2" # This shouldn't occur, Success should be FALSE.
            } else {
              "Oops"
            },
            Success = TRUE
          )
          if (newEvent$Type == "Present") {
            # Instead, just replace the current event to present.
            # (Excess work, but clearer code!)
            event[j, ]$Type <- "Present"
          } else {
            event <- rbind(event, newEvent)
          }
        }
      }
    }


    # Proceed along changes, adding them to events.
    if (nrow(changesdf) > 0) {
      for (j in 1:nrow(changesdf)) {
        thisChange <- changesdf[j, ]
        thisEvent <- event %>% dplyr::filter(
          Species == thisChange$species,
          Environment == thisChange$environments
        )
        if (nrow(thisEvent) == 0) {
          # Event definitely not recorded, create a new one and add it.
          event <- rbind(event, data.frame(
            Times = time,
            Species = thisChange$species,
            Environment = thisChange$environments,
            Type = switch(thisChange$type + 2, # -1:1 + 2 -> 1:3
                          "Dynamic Loss",
                          "Oops",
                          "Dispersal"
            ),
            Success = TRUE
          ))
        } else if (nrow(thisEvent) == 1 &&
                   thisEvent$Species == thisChange$species &&
                   thisEvent$Environment == thisChange$environments &&
                   ((thisChange$type == 1 && thisEvent$Type == "Arrival") ||
                    (thisChange$type == -1 && thisEvent$Type == "Extinct")
                   )) {
          # Trivial, everything is already recorded.
          # Note this also should already have been caught!
          next()
        } else {
          # Event recorded, but opposite of what is observed -- undone twice.
          # Assuming that events take place at distinct times, which should be
          # guaranteed from the code we are using, this shouldn't happen.
          # Discretisation means that it can happen though, and we can have
          # something considered a false positive: arrival events are tracked as
          # positive (because the species is added to the local population) even
          # if the population is declining and about to be removed.
          if (nrow(thisEvent) == 1 &&
              thisEvent$Species == thisChange$species &&
              thisEvent$Environment == thisChange$environments &&
              (thisChange$type == -1 && thisEvent$Type == "Arrival")
              ) {
            # "Poorly Timed False Positive Arrival"
            event <- rbind(event, data.frame(
              Times = time,
              Species = thisChange$species,
              Environment = thisChange$environments,
              Type = "Dynamic Loss",
              Success = TRUE
            ))
            # We also should account for that the Arrival is actually a Present.
            event[event$Times == thisEvent$Times &
                    event$Species == thisEvent$Species &
                    event$Environment == thisEvent$Environment &
                    event$Type == thisEvent$Type ,]$Type <- "Present"
          } else {
            stop("Double-check assumptions: a 'triple event' occurred.")
          }
        }
      }
    }

    return(event)
  },
  mat = ColExtMatrix,
  binmat = binaryMatrix,
  events = sim$Events)

  allEvents <- dplyr::bind_rows(allEvents) %>% dplyr::arrange(Times)
  # Adding in typings is then dependent on how we are breaking up affinity
  # and what information from the pool is desired. I think it might make
  # more sense to do that outside of this function.
  return(allEvents)
}

### Biodiversity Intactness: ##################################################

# if (calculationsBII) {
# Adjusting
# https://adrianadepalma.github.io/BII_tutorial/bii_example.html#Compositional_Similarity
# to match our formatting.
# Note we get both the old and new versions.

getJacAbSym <- function(s1, s2, s1Control, s2Control, data){

  # "get the list of species that are present in site 1 (i.e., their abundance was greater than 0)"
  s1species <- data %>%

    # "filter out[sic] the SSBS that matches s1" (pristine)
    # (Note extra filter to make sure pristine/control because timeseries.)
    dplyr::filter(Patch == s1,
                  Control == s1Control ) %>%

    # "filter out[sic] the species where the Measurement (abundance) is greater than 0"
    dplyr::filter(Abundance > 0) %>%

    # "get the unique species from this dataset"
    dplyr::distinct(ID) %>%

    # "pull out the column into a vector"
    dplyr::pull(ID)

  # "for site 2, get the total abundance of species that are also present in site 1"

  s2abundance_s1species <- data %>%

    # "filter out[sic] the SSBS that matches s2"
    # (If timeseries, s1 might be s2, in which case grab not control data.)
    dplyr::filter(Patch == s2, Control == s2Control ) %>%

    # "filter out[sic] the species that are also present in site 1"
    dplyr::filter(ID %in% s1species) %>%

    # "pull out the Measurement into a vector"
    dplyr::pull(Abundance) %>%

    # "calculate the sum"
    sum()

  # "calculate the total abundance of all species in site 2"
  s2_sum <- data %>%

    # "filter out[sic] the SSBS that matches s2"
    dplyr::filter(Patch == s2) %>%

    # "pull out the measurement column (the abundance)"
    dplyr::pull(Abundance) %>%

    # "calculate the sum"
    sum()


  # "Now calculate the compositional similarity
  # this is the number of individuals of species also found in s1, divided by the total abundance in s2
  # so that equates to the proportion of individuals in s2 that are of species also found in s1"

  sor <- s2abundance_s1species / s2_sum


  # "if there are no taxa in common, then sor = 0
  # if abundances of all taxa are zero, then similarity becomes NaN."
  return(sor)

}

# NOTE: SIMILARITIES
get_bray <- function(s1, s2, s1Control, s2Control, data){

  require(betapart)

  sp_data <- data %>%

    # filter patches to the pair we're interested in
    dplyr::filter(paste(Patch, Control) %in%
                    c(paste(s1, s1Control),
                      paste(s2, s2Control))
    ) %>%

    dplyr::mutate(Patch = paste(Patch, Control)) %>%

    # pull out the site name, species name and abundance information
    dplyr::select(Patch, ID, Abundance) %>%

    # pivot so that each column is a species and each row is a Patch
    tidyr::pivot_wider(names_from = ID, values_from = Abundance,
                       values_fill = 0 # For some reason missing in orig. func?
    ) %>%

    # set the rownames to the Patch and then remove that column
    tibble::column_to_rownames("Patch")

  # if one of the sites doesn't have any individuals in it
  # i.e. the row sum is 0
  if(sum(rowSums(sp_data) == 0, na.rm = TRUE) == 1){
    # then the similarity between sites should be 0
    bray <- 0
    # if both sites have no individuals
  }else if(sum(rowSums(sp_data) == 0, na.rm = TRUE) == 2){
    # then class the similarity as NA
    bray <- NA
    # otherwise if both sites have individuals, calculate the balanced bray-curtis
    # as similarity (1-bray)
  }else{
    bray <- 1 -
      betapart::bray.part(sp_data) %>%
      purrr::pluck("bray.bal") %>%
      purrr::pluck(1)
  }
  return(bray)
}

# NOTE: DISSIMILARITIES
get_bray_all <- function(s1, s2, s1Control, s2Control, data){
  require(betapart)
  sp_data <- data %>%
    # filter patches to the pair we're interested in
    dplyr::filter(paste(Patch, Control) %in%
                    c(paste(s1, s1Control),
                      paste(s2, s2Control))
    ) %>%
    dplyr::mutate(Patch = paste(Patch, Control)) %>%
    # pull out the site name, species name and abundance information
    dplyr::select(Patch, ID, Abundance) %>%
    # pivot so that each column is a species and each row is a Patch
    tidyr::pivot_wider(names_from = ID, values_from = Abundance,
                       values_fill = 0 # For some reason missing in orig. func?
    ) %>%
    # set the rownames to the Patch and then remove that column
    tibble::column_to_rownames("Patch")

  bray <- betapart::bray.part(sp_data) %>% unlist() %>% t() %>% as.data.frame()
  return(bray)
}

inv_logit <- function(f, a) {
  a <- (1 - 2*a)
  (a * (1 + exp(f) ) + (exp(f) - 1)) / (2 * a * (1 + exp(f) ))
}
# }
