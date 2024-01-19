# Introduction: ###############################################################
# As a sequel to TimeSpaceAndTimeSeries-1-Bootstrap.R and -1-Intervention.R,
# we are now introducing a simulation based intervention.
# Please see the previous files for some design choices, although we aim to
# improve design at each stage.
# For file management, we'll split up the sampling into a second file.
# We'll pass through as records relevant simulation/intervention details.
# TODO: Consider if the introduction should be split off into a file 0.

# Create a Master seed, which we'll use to generate simulation seeds.
# runif(1) * 1e8
simulationsNumber <- 1
simulationsDictionaryChoice <-
# 1 # Example set with infinite disp. res.
 2 # Single Run with Interaction Matrix with PConstrain set to c(0, pnorm(-1))
# 3 # Single Run with Interaction Matrix with PConstrain set to c(pnorm(+1), 1)

interventionChoice <-
# 1. # The change does not occur, but the scientists don't know that!
 2. # The change occurs, is discontinuous, and local, e.g. forest fire.
# 3. # The change is slow and local, e.g. gradual deforestation.
# 4. # The change occurs on all patches, e.g. climate change.
# 5. # Two changes, one following 2, the other 4.
# 6. # Two changes, one following 3, the other 4.
# 2, 3 are interaction matrix changes.
# 4 is inspired from Amarasekare's 2019(?) paper on temperature -> pred-prey.

interventionParameterChoice <- # Valid for some interventions
# 1. # 2,3: Intervention Matrix with Base Parameters
# 2. # 2,3: Intervention Matrix with High Variance:Mean Ratio (0.1 -> 0.2)
# 3. # 2,3: Intervention Matrix with Low Variance:Mean Ratio (0.1 -> 0)
# 4. # NI, 4: Base parameters for changing LM1996.
# 5. # NI, 5: Base parameters for changing Amarasekare2019.
 6. # 2,3: Intervention Matrix with Constrained Sampling (pnorm(1), 1)
# 7. # 2,3: Intervention Matrix with Constrained Sampling (0, pnorm(-1))

interventionTimeSpan <- 200 # for interventionChoices 3 and 4 only.
# 20 seems like a good number to afford strong sampling on the transition
# with still consistent sampling after the transition is effectively finished
# for samplingMaxTime = 250 and Quantity = 100.

#TODO MOVE TO SEQUEL.
# Note: these are on the characteristic scale, rather than the normal scale.
samplingMaxTime <- 250 # This is roughly twice the observed burn-in time from 0.
samplingQuantity <- 100 # Not guaranteed!
samplingTimeScaleLogarithmic <- TRUE
# calculationsPlotLong <- FALSE

cores <- 1

# NOTE: Might be able to improve by splitting and collapsing Sources 1 and 2.
# Goal => create a list of simulations, each represented as a list of files.
simulationsDictionary <- data.frame(
  SourceSimulations = c(
    # Sets of files, as a string, that can be split with ", ".
    # NOTE: Files should be generated with same parameters for comparability.
    # NOTE: Use file.path if you need files in directories.
    paste0(c(
      file.path("Data_2023-07-06",
                "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData"),
      file.path("Data_2023-09-23",
                "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData"),
      file.path("Data_2023-09-24",
                "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData"),
      file.path("Data_2023-09-25",
                "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData"),
      file.path("Data_2023-09-26",
                "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData")
    ), collapse = ", "),
    paste0(c(
      file.path("Data_2024-01-17",
                "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData")
    ), collapse = ", "),
    paste0(c(
      file.path("Data_2024-01-18",
                "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData")
    ), collapse = ", ")
  ),
  SourceSimulations2 = c(
    # Sets of files, as a string, that can be split with ", ".
    # NOTE: Files here should be second halves of files in SourceSimulations.
    # NOTE: Use "NA" if a corresponding second half does not exist.
    # NOTE: Use file.path if you need files in directories.
    paste0(c("NA", "NA", "NA", "NA", "NA"), collapse = ", "),
    paste0(c("NA"), collapse = ", "),
    paste0(c("NA"), collapse = ", ")
  ),
  SourcePools = c(
    # Sets of files, as a string, that can be split with ", ".
    # NOTE: Should correspond one-to-one with the SourceSimulations.
    paste0(c(
      file.path("Data_2023-07-06",
                "MNA-ExampleOutcome-PoolMats-Env10.RData"),
      file.path("Data_2023-09-23",
                "MNA-ExampleOutcome-PoolMats-Env10.RData"),
      file.path("Data_2023-09-24",
                "MNA-ExampleOutcome-PoolMats-Env10.RData"),
      file.path("Data_2023-09-25",
                "MNA-ExampleOutcome-PoolMats-Env10.RData"),
      file.path("Data_2023-09-26",
                "MNA-ExampleOutcome-PoolMats-Env10.RData")
    ), collapse = ", "),
    paste0(c(
      file.path("Data_2024-01-17",
                "MNA-ExampleOutcome-PoolMats-Env10.RData")
    ), collapse = ", "),
    paste0(c(
      file.path("Data_2024-01-18",
                "MNA-ExampleOutcome-PoolMats-Env10.RData")
    ), collapse = ", ")
  ),
  Seeds = c(# runif(1) * 1e8 # Unsure if this should be in the interventionDic.
    69148611,
    35181045,
    79180721
  )
)[simulationsDictionaryChoice, ]

interventionDictionary <- data.frame(
  DynamicsFunction = c(
    "RMTRCode2::PerCapitaDynamics_Type1",
    "RMTRCode2::PerCapitaDynamics_Type1",
    "RMTRCode2::PerCapitaDynamics_Type1",
    NA,
    NA,
    "RMTRCode2::PerCapitaDynamics_Type1",
    "RMTRCode2::PerCapitaDynamics_Type1"
  ),
  InteractionMatrixFunction = c(
    "RMTRCode2::LawMorton1996_CommunityMat",
    "RMTRCode2::LawMorton1996_CommunityMat",
    "RMTRCode2::LawMorton1996_CommunityMat",
    NA,
    NA,
    "RMTRCode2::LawMorton1996_CommunityMat",
    "RMTRCode2::LawMorton1996_CommunityMat"
  ),
  InteractionMatrixArguments = c(
    # I can't find a better way to create partials here then as strings.
    # Suggestions welcome!
    "Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.1)",
    "Parameters = c(0.01, 10, 0.5, 0.2, 100, 0.2)",
    "Parameters = c(0.01, 10, 0.5, 0.2, 100, 0)",
    NA,
    NA,
    "ConstrainP = c(pnorm(1), 1)",
    "ConstrainP = c(0, pnorm(-1))"
  ),
  DispersalResistance = c(
    Inf, Inf, Inf, NA, NA, Inf, Inf
  ),
  DispersalFormat = c(
    "Ring", "Ring", "Ring", NA, NA, "Ring", "Ring"
  )
)[interventionParameterChoice, ]

# Libraries: ##################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(patchwork)

library(RMTRCode2)

library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

directory <- '.' # VariousExperiments folder.
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))

# Parallelization: ############################################################
if (cores > 1) {
  clust <- parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  `%op%` <- foreach::`%dopar%`
} else {
  `%op%` <- foreach::`%do%`
}

# Files: ######################################################################
directory <- "." # Should be "VariantExperiments"
datfolder <- file.path(
  directory,
  paste0(
    "TSTS_Simulations_", # Separate the Name (TSTS) and Type (Simulations)
    simulationsDictionaryChoice, "-", # Bundle the Triple
    interventionChoice, "-",
    interventionParameterChoice, "_", # Separate the Date
    Sys.Date())
)
if (!dir.exists(datfolder)) {
  dir.create(datfolder)
}

# Convert the simulationsDictionary into proper file pairs.
.sSeed <- simulationsDictionary$Seeds
.sFile1 <- strsplit(simulationsDictionary$SourceSimulations, split = ", ")[[1]]
.sFile2 <- strsplit(simulationsDictionary$SourceSimulations2, split = ", ")[[1]]
.sPool <- strsplit(simulationsDictionary$SourcePools, split = ", ")[[1]]

stopifnot(length(.sFile1) == length(.sFile2),
          length(.sFile1) == length(.sPool))

.s <- foreach::foreach(
  id = seq_along(.sFile1),
  file1 = iterators::iter(.sFile1),
  file2 = iterators::iter(.sFile2),
  fpool = iterators::iter(.sPool)
) %op% {
  simulation <- loadSimulation(file1, if(file2 != "NA") file2)

  # Make sure we clear out any extinct abundances:
  simulation$Abundance[, -1] <- ifelse(
    simulation$Abundance[, -1] <
      simulation$Parameters$EliminationThreshold,
    0, simulation$Abundance[, -1]
  )

  # Interesting problem: we need everything to be on the same time scales.
  # We already use the reaction times as our best proxy/control of the time
  # scale, so we'll reformat all of the times using those.
  simulation$Abundance[, 1] <-
    simulation$Abundance[, 1] / simulation$ReactionTime
  simulation$Events$Times <-
    simulation$Events$Times / simulation$ReactionTime

  load(fpool)

  # Some files we might be interested in have multiple pools bundled.
  # This code uses the file names to extract the correct one.
  if (exists("pools") && !exists("Pool")) {
    IDNumbers <- sub(basename(file1), pattern = ".RData", replacement = "")
    IDNumbers <- strsplit(IDNumbers, split = "-", fixed = TRUE)[[1]]
    IDNumbers <- IDNumbers[
      (which(IDNumbers == "Prepared") + 1):length(IDNumbers)
      ]
    IDNumbers <- as.numeric(IDNumbers)
    Pool <- pools[[cases$Parameters[IDNumbers[2]]]][[
      cases$System[IDNumbers[2]]
      ]]
    # Why 2? 1 == File / Main Case, 2 == row of cases (derived from CSV row)
    # 3 == The seeds used, 4 == which part of the simulation was saved.
    # Note 3 & 4 are optional if all simulations of a row were saved together.
  }

  lastTimeSampleable <- simulation$Events$Times[
    length(simulation$Events$Times)
    ]
  lastTimeSampleable <-
    lastTimeSampleable - samplingMaxTime * 2
  # The * 2 makes sure we're also outside the burnout as well as that we have
  # enough time for sampling.

  return(list(
    ID = id,
    File1 = file1,
    File2 = file2,
    FPool = fpool,
    NEnvs = simulation$NumEnvironments,
    NSpec = nrow(Pool),
    Simulation = simulation,
    Pool = Pool,
    InteractionMatrices = InteractionMatrices,
    IntervenableTimes = c(
      125, # Identified by inspection.
      lastTimeSampleable
    )
  ))
}

# Setup and Framing: ##########################################################
# (Supposing omniscience == know all information about the present, !future:)
# Two researchers, an omniscient recorder, and a dimension hopping omniscient
# recorder come upon a set of sites after having been told that the
# local government has designated them for some form of land use change
# (e.g. fertilizer, w/e). The first researcher arrives with some time
# before the change is implemented and will compare the sites before and after
# the land use change. The second researcher is too late and needs to compare
# instead with surrounding sites that did not undergo the land use change.
# Meanwhile the omniscient recorder knows the actual changes (no sampling), and
# the dimension hopper *also* knows what would have happened if the land use
# change did not, in fact, occur.
#
# Furthermore, suppose that there are 6 universes, each with a separate change:
#   1. The change does not occur, but the scientists don't know that!
#   2. The change occurs, is discontinuous, and local, e.g. forest fire.
#   3. The change is slow and local, e.g. gradual deforestation.
#   4. The change occurs on all patches, e.g. climate change.
#   5. Two changes, one following 2, the other 4.
#   6. Two changes, one following 3, the other 4.
#
# How do the four observers' results differ between each other in each universe?
# Consider the impact of sites monitored (incl. quantity),
#                        monitoring effectiveness (e.g. sampling intensity),
#                        and length of monitoring.
#
# For comparability, interventions are all considered to start at the same time,
# but the samples we compare may be midway through or after interventions, esp.
# when comparing slow to fast or immediate interventions.

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

### Interventions: ############################################################
actTrivially <- function(elem) {force(elem); function(...){elem}}

switchMatrices <- function(matrix1, matrix2, switchtime) {
  stopifnot(dim(matrix1) == dim(matrix2))
  force(matrix1);force(matrix2);force(switchtime)
  function(t, ...) {
    if (t < switchtime) {
      matrix1
    } else if (t > switchtime) {
      matrix2
    } else {
      (matrix1 + matrix2)/2
    }
  }
}

switchMatrixLists <- function(matrixlist1, matrixlist2, switchtimes) {
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(length(matrixlist1) == length(switchtimes) ||
              length(switchtimes) == 1)

  Matrix::bdiag(lapply(seq_along(matrixlist1), function(i, m1, m2, st) {
    switchMatrices(m1[[i]], m2[[i]], st[i])
  },
  m1 = matrixlist1, m2 = matrixlist2,
  st = if (length(switchtimes) == 1) {
    rep(switchtimes, length(matrixlist1))
  } else {switchtimes}
  ))
}

switchMatrixLists2 <- function(matrixlist1, matrixlist2, switchtime) {
  # Note: only accepts 1 switchtime.
  #       multiple would require steps.
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(length(matrixlist1) == length(switchtime) ||
              length(switchtime) == 1)
  stopifnot(isTRUE(all.equal(lapply(matrixlist1, dim),
                             lapply(matrixlist2, dim))))

  matrix1 <- Matrix::bdiag(matrixlist1)
  matrix2 <- Matrix::bdiag(matrixlist2)

  switchMatrices(matrix1, matrix2, switchtime = switchtime)
}

# Choose a sigmoidal interpolation between the matrices.
# Problem: we want no deviation at 0 or t = timespan, but not too much in mid.
# Crit. Value = first numb. in tanh. 4 => 2% dev at edge, ~75% done in mid. 50%.
# Require: Output is a function.
# Require: Output: t = timespan / 2 => out = (mat1 + mat2) / 2.
interpolateMatrices <- function(matrix1, matrix2, timespan, switchtime = 0) {
  stopifnot(dim(matrix1) == dim(matrix2))
  force(matrix1);force(matrix2);force(timespan);force(switchtime)
  function(t, ...) {
    # if (t < switchtime) {
    # matrix1
    # } else {
    matrix1 + (matrix2 - matrix1) * (
      tanh(4 * ( (t - switchtime) / timespan - 0.5)) + 1
    ) / 2
    # }
  }
}

interpolateMatrixLists <- function(
  matrixlist1, matrixlist2, timespans, switchtimes
) {
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(length(matrixlist1) == length(timespans) ||
              length(timespans) == 1)
  stopifnot(length(matrixlist1) == length(switchtimes) ||
              length(switchtimes) == 1)

  Matrix::bdiag(lapply(
    seq_along(matrixlist1), function(i, m1, m2, ts, st) {
      interpolateMatrices(m1[[i]], m2[[i]], ts[i], st[i])
    },
    m1 = matrixlist1, m2 = matrixlist2,
    ts = if (length(timespans) == 1) {
      rep(timespans, length(matrixlist1))
    } else {timespans},
    st = if (length(switchtimes) == 1) {
      rep(switchtimes, length(matrixlist1))
    } else {switchtimes}
  ))
}

interpolateMatrixLists2 <- function(matrixlist1, matrixlist2,
                                    timespan, switchtime) {
  stopifnot(length(matrixlist1) == length(matrixlist2))
  stopifnot(#length(matrixlist1) == length(switchtime) ||
              length(switchtime) == 1)
  stopifnot(#length(matrixlist1) == length(timespan) ||
              length(timespan) == 1)
  stopifnot(isTRUE(all.equal(lapply(matrixlist1, dim),
                             lapply(matrixlist2, dim))))

  matrix1 <- Matrix::bdiag(matrixlist1)
  matrix2 <- Matrix::bdiag(matrixlist2)

  interpolateMatrices(matrix1, matrix2,
                      timespan = timespan, switchtime = switchtime)
}

### Retreive the stored dynamics function: ####################################
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

dynFunc <- retrieveFunction(interventionDictionary$DynamicsFunction)

# Need to use stackoverflow.com/a/47012149 to convert our arguments to a list.
callMeMaybe2 = function(listofcharargs) {
  Charargs = unlist(listofcharargs)
  if(is.null(Charargs)) return(alist())
  eval(parse(
    text = paste0("alist(",
                  paste(parse(text = Charargs),
                        collapse = ","),")")
  ))
}

intMatFunc <- purrr::partial(
  retrieveFunction(interventionDictionary$InteractionMatrixFunction),
  ... =, # otherwise, partialised arguments occur first.
  !!!callMeMaybe2(interventionDictionary$InteractionMatrixArguments)
  # !!! suggested by Bing Chat. Not obvious how it evaluates the list.
  # Prompt 'R purrr, using a list of arguments to partialize a function'.
  # 2024/01/19
)


# Intervention Simulation: ####################################################
# If we've coded things correctly, we should be able to re-run our original
# RMTRCode2::MultipleNumericalAssembly_Dispersal(...) simulation.
# Recall that this function
#   1) can start from arbitrary initial populations and times,
#   2) does not need the interaction matrices themselves,
#   3) only needs the pool size to not change, the rest of the pool can,
#   4) can take the fixed events list, and
#   5) can take any PerCapitaDynamics that accepts (t, y, parms) arguments.
# Downsides:
#   1) we'll have to fix a characteristic rate. Pick the most frequent sampling.
#   2) The DispersalMatrix will need to remain fixed.

# In order to charge ahead then, we need to
#   1) pick a base file and history,
#   2) reformat it to meet the MNAD criteria,
#   3) create a PerCapitaDynamics function to pass to MNAD.

# WARNING: doRNG generates random seeds in a way that makes the *order* of the
# iterations (but not their order of evaluation) important.
# Changing the order of iterations (e.g. id = simulationsNumber:1) WILL change
# the random seed assigned to a given iteration!!!
results <- foreach::foreach(
  id = 1:simulationsNumber,
  s = iterators::iter(.s, recycle = TRUE), # File and History (not random!).
  .options.RNG = simulationsDictionary$Seeds,
  .combine = "rbind",
  .packages = c("dplyr")
) %dorng% {
  fullID <- paste0(simulationsDictionaryChoice, "-",
                   interventionChoice, "-",
                   interventionParameterChoice, "-",
                   id, "-",
                   s$ID)
  filename <- file.path(datfolder, paste0(
    "TSTS_Simulation_", fullID, ".RData"
  ))
  if (file.exists(filename)) {
    load(filename)
    if(exists("interventionSimulation")) {
      return(interventionSimulation)
    } else {
      stop(paste("Bad file:", filename))
    }
  }

  # Pick a random patch as control, rest as experiment.
  # is adding 1:5 (or w/e) okay? Yes, we're assuming contiguous patches.
  control <- ((sample.int(s$NEnvs, 1) + 1:(s$NEnvs / 2) ) %% s$NEnvs) + 1
  # print(paste(bootstrapID, ":", toString(control)))
  experiment <- c(1:s$NEnvs)[!c(1:s$NEnvs) %in% control]
  # print(paste(bootstrapID, ":", toString(experiment)))

  # Pick a random time to start the intervention and sampling from.
  timeIntervention <- runif(
    n = 1, min = s$IntervenableTimes[1], max = s$IntervenableTimes[2]
  )

  timeInterventionRow <- which.max(
    s$Simulation$Abundance[, 1] > timeIntervention
  ) - 1

  stopifnot(timeInterventionRow > 0)

  # Create the intervention per capita dynamics.
  # DynFunc follows archetype with formal arguments:
  #  1) ReproductionRate, 2) InteractionMatrices, 3) NumEnvironments.
  # TODO: Add support (how? partial completion?) for other arguments,
  #       e.g. Mutualistic1's SpeciesTypes and Saturations.
  # Below is roughly the idea, but we need to act differently between
  # the control and the experimental cases.
  # TODO: Turn this into a function?
  # NOTE: Times inside this function are on the original scale,
  #       NOT the characteristic.
  PerCapDyn <- switch(
    interventionChoice,
    # 1. # The change does not occur, but the scientists don't know that!
    dynFunc(actTrivially(s$Pool$ReproductionRate),
            actTrivially(Matrix::bdiag(s$InteractionMatrices$Mats)),
            s$NEnvs),
    # 2. # The change occurs, is discontinuous, and local, e.g. forest fire.
    {
      experimentMats <- RMTRCode2::CreateEnvironmentInteractions(
        Pool = s$Pool, NumEnvironments = length(experiment),
        ComputeInteractionMatrix = intMatFunc
      )
      mats <- list()
      eMindex <- 1
      for(i in 1:s$NEnvs) {
        if (i %in% control) {
          mats[[i]] <- s$InteractionMatrices$Mats[[i]]
        } else if (i %in% experiment) {
          mats[[i]] <- experimentMats$Mats[[eMindex]]
          eMindex <- eMindex + 1
        } else {
          stop("An environment is not allocated to control or experiment.")
        }
      }
      dynFunc(actTrivially(s$Pool$ReproductionRate),
              switchMatrixLists2(
                s$InteractionMatrices$Mats,
                mats,
                switchtime = timeIntervention * s$Simulation$ReactionTime
              ),
              s$NEnvs)
    },
    # 3. # The change is slow and local, e.g. gradual deforestation.
    {
      experimentMats <- RMTRCode2::CreateEnvironmentInteractions(
        Pool = s$Pool, NumEnvironments = length(experiment),
        ComputeInteractionMatrix = intMatFunc
      )
      mats <- list()
      eMindex <- 1
      for(i in 1:s$NEnvs) {
        if (i %in% control) {
          mats[[i]] <- s$InteractionMatrices$Mats[[i]]
        } else if (i %in% experiment) {
          mats[[i]] <- experimentMats$Mats[[eMindex]]
          eMindex <- eMindex + 1
        } else {
          stop("An environment is not allocated to control or experiment.")
        }
      }
      dynFunc(actTrivially(s$Pool$ReproductionRate),
              interpolateMatrixLists2(
                s$InteractionMatrices$Mats,
                mats,
                timespan = interventionTimeSpan * s$Simulation$ReactionTime,
                switchtime = timeIntervention * s$Simulation$ReactionTime
              ),
              s$NEnvs)
    },
    # 4. # The change occurs on all patches, e.g. climate change.
    {stop("Not Implemented"); DynFunc()},
    # 5. # Two changes, one following 2, the other 4.
    {stop("Not Implemented"); DynFunc()},
    # 6. # Two changes, one following 3, the other 4.
    {stop("Not Implemented"); DynFunc()}
  )

  # Perform the intervention simulation.

  ## Events
  # Note formatting is a list containing a data.frame named Events.
  # Note that this is in normalized time units.
  eventsPostIntervention <- list(
    Events = s$Simulation$Events %>% dplyr::filter(
      Times > s$Simulation$Abundance[timeInterventionRow, 1]
      ))
  # Why not timeIntervention? To make sure that we don't miss out on an event.
  # Possibly unnecessary.
  eventsPostIntervention$Events$Success <- NA

  # TODO: Move this out of the loop to make it more accessible as an
  #       intervention choice. We'll then want to convert to a function
  #       in order to enable varying temporal changes.
  ## Distance Matrix
  if (interventionDictionary$DispersalFormat == "None") {
    DistanceMatrix <- Matrix::sparseMatrix(
      i = s$NEnvs, j = s$NEnvs, x = 0)
  }
  if (interventionDictionary$DispersalFormat == "Ring" ||
      interventionDictionary$DispersalFormat == "Line")
    DistanceMatrix <- Matrix::bandSparse(
      s$NEnvs, k = c(-1, 1),
      diagonals = list(
        rep(interventionDictionary$DispersalResistance, s$NEnvs - 1),
        rep(interventionDictionary$DispersalResistance, s$NEnvs - 1))
    )
  if (interventionDictionary$DispersalFormat == "Ring") {
    DistanceMatrix[s$NEnvs, 1] <-
      interventionDictionary$DispersalResistance
    DistanceMatrix[1, s$NEnvs] <-
      interventionDictionary$DispersalResistance
  }
  if (interventionDictionary$DispersalFormat == "Full") {
    DistanceMatrix <- matrix(1, nrow = s$NEnvs, ncol = s$NEnvs)
    diag(DistanceMatrix) <- 0
  }

  DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(
    EnvironmentDistances = DistanceMatrix,
    SpeciesSpeeds = if(is.null(s$Pool$Speed)) {
      rep(1, nrow(s$Pool))
    } else {
      rep(SpeciesSpeeds, nrow(s$Pool))
    }
  )

  ## Pick fastest sampling characteristic rate.
  if (exists("experimentMats")) {
    charRate <- max(1/s$Simulation$ReactionTime,
                    unlist(lapply(experimentMats$Mats, function(mat) {
                      abs(eigen(mat, only.values = TRUE)$values)
                    })))
  } else {
    charRate <- 1/s$Simulation$ReactionTime
  }

  # Undo the time scale normalisation on the event and simulation times.
  # This is because transforms weren't performed on the simulation parameters
  # (i.e. the species interactions, dispersal, and species properties).
  s$Simulation$Abundance[, 1] <-
    s$Simulation$Abundance[, 1] * s$Simulation$ReactionTime
  eventsPostIntervention$Events$Times <-
    eventsPostIntervention$Events$Times * s$Simulation$ReactionTime

  ## Simulation
  interventionSimulation <- RMTRCode2::MultipleNumericalAssembly_Dispersal(
    PopulationInitial = s$Simulation$Abundance[timeInterventionRow, -1],
    TimeInitial = s$Simulation$Abundance[timeInterventionRow, 1],
    # This induces a slight delay compared to starting AT the intervention.
    Pool = s$Pool,
    NumEnvironments = s$NEnvs,
    Events = eventsPostIntervention,
    PerCapitaDynamics = PerCapDyn,
    DispersalMatrix = DispersalMatrix,
    EliminationThreshold = s$Simulation$Parameters$EliminationThreshold,
    ArrivalDensity = s$Simulation$Parameters$ArrivalDensity,
    MaximumTimeStep = s$Simulation$Parameters$MaximumTimeStep,
    BetweenEventSteps = s$Simulation$Parameters$BetweenEventSteps,
    CharacteristicRate = charRate,
    Verbose = TRUE,
    # Using the ellipsis pass through feature:
    ReactionTimePrev = s$Simulation$ReactionTime,
    FullID = fullID,
    Intervention = list(# NOTE: characteristic scale!
      Patches = experiment,
      Time = timeIntervention,
      OriginalRow = timeInterventionRow,
      TimeSpan = if (interventionChoice == 3) interventionTimeSpan
    )
  )

  # NOTE: For time consistency, we need to use the characteristic timescale,
  # as opposed to the evaluation timescale that is newly resolved.
  # TODO: Ask Jon, Susan about this choice to verify.
  interventionSimulation$Abundance[, 1] <-
    interventionSimulation$Abundance[, 1] / s$Simulation$ReactionTime
  interventionSimulation$Events$Times <-
    interventionSimulation$Events$Times / s$Simulation$ReactionTime


  # Save results so we can sample in the next file.
  save(interventionSimulation, file = filename)

  return(interventionSimulation)
}

# Useful double check:
# par(mfrow = c(5, 5))
# if (!exists("timeInterventionRow")) {
#   timeInterventionRow <- which.max(
#     s$Simulation$Abundance[, 1] >
#       interventionSimulation$Ellipsis$Intervention$Time
#   ) - 1
# }
# for(j in which(
#   apply(
#     s$Simulation$Abundance[timeInterventionRow:nrow(s$Simulation$Abundance),],
#     MARGIN = 2,
#     FUN = function(x) any(x > s$Simulation$Parameters$EliminationThreshold)) |
#   apply(
#     interventionSimulation$Abundance, MARGIN = 2,
#     FUN = function(x) any(x > s$Simulation$Parameters$EliminationThreshold)))
# ) {
#
#   plot(
#     s$Simulation$Abundance[timeInterventionRow:nrow(s$Simulation$Abundance), 1],
#     s$Simulation$Abundance[timeInterventionRow:nrow(s$Simulation$Abundance), j],
#     type = "l", col = "red", main = j,
#     ylim = c(0, max(
#       s$Simulation$Abundance[timeInterventionRow:nrow(s$Simulation$Abundance), j],
#       interventionSimulation$Abundance[, j]
#       )))
#   lines(
#     interventionSimulation$Abundance[, 1],
#     interventionSimulation$Abundance[, j],
#     col = "blue", lty = 2)
# }
