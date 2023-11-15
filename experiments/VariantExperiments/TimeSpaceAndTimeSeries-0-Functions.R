
# Functions: ##################################################################
nonzero <- function(vec) {
  vec[vec!=0]
}

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

# Okay, next, one common topic is essentially native vs invasive.
# In truth all species are more or less the same, originating from the true
# regional pool. Our researchers wouldn't know this.
# Instead, they would both use their control to estimate the natives.
computeSpeciesInControl <- function(sampling, Time = "Time") {
  controlSpecies <- sampling %>% dplyr::filter(
    Control == "Control"
  ) %>% dplyr::pull(
    SamplingIDs
  ) %>% strsplit(
    ", ", fixed = TRUE
  ) %>% unlist(
  ) %>% unique()

  sampling %>% dplyr::group_by(
    dplyr::across(dplyr::all_of(Time)), Patch
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      splitSamplingIDs <- strsplit(x$SamplingIDs, ", ", fixed = T)[[1]]
      iDsInControl <- splitSamplingIDs %in% controlSpecies

      x %>% dplyr::mutate(
        SamplingIDsNative = toString(splitSamplingIDs[iDsInControl]),
        SamplingAbundanceNative = sum(iDsInControl),
        SamplingAlphaNative = length(unique(splitSamplingIDs[iDsInControl])),
        SamplingIDsInvasive = toString(splitSamplingIDs[!iDsInControl]),
        SamplingAbundanceInvasive = sum(!iDsInControl),
        SamplingAlphaInvasive = length(unique(splitSamplingIDs[!iDsInControl]))
      )
    }
  )
}

addDistanceColumns <- function(bootstrapSamples, mindig = 1, Time = "Time") {
  bootstrapSamples %>% dplyr::mutate(
    ##### Temporal Distance: ####################################################
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
      ##### Spatial Distance: #####################################################
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
      ##### Define control species: ###############################################
    ) %>% dplyr::group_by(
      DistanceFromCenterExpRev, Bootstrap, TimeGapNumber, Type
    ) %>% dplyr::group_modify(
      .f = ~ computeSpeciesInControl(.x, Time = Time)
    ) %>% dplyr::ungroup(
    )
}

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
