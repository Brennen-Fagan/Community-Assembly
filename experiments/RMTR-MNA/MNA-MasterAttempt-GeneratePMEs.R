# Abstract: ####################################################################
# This file takes the output of "MNA-MasterAttempt-GenerateCases.R" and converts
# it to a set of pools, matrices and events to prevent bottlenecking during a
# submission to Viking.

library(RMTRCode2)

# https://stackoverflow.com/a/15373917
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

# Initialisation: ##############################################################
# Identify the files to take.
thisDirectory <- dirname(thisFile())
theCases <- dir(thisDirectory,
                pattern = "MNA-Master[A-Z]-Cases[.]csv",
                full.names = TRUE)
theParams <- dir(thisDirectory,
                 pattern = "MNA-MasterParameters[.]RData",
                 full.names = TRUE)

# Identify objects to load.
params <- c("systemBase", "systemMods", "thisSeparator")
for (p in params) {assign(p, NULL)}

# Load all parameters.
loaded <- unlist(lapply(theParams, load, .GlobalEnv))

# Make sure they are all loaded correctly.
stopifnot(
  params %in% loaded,
  !is.null(params)
)

# Repeat for the cases.
# We assign each case an environment within which we
# will do our calculations (create pools, matrices, events).
# This means most calculations will look like zips.
# All calculations will also be ignorant of the number
# of cases that we are actually considering.
theEnvironments <- lapply(theCases, function(case) {
  tempenv <- new.env()
  assign("cases", read.csv(case, as.is = TRUE), envir = tempenv)
  return(tempenv)
})

# Generate Pools: #############################################################
theEnvironments <- lapply(
  theEnvironments, function(env) {
    with(env, {
      # pools will be a row-associated list of seed-associated lists of data frames.
      pools <- lapply(
        1:nrow(cases), function(row, b, m, cs) {
          crow <- cases[row, ]

          # Using inheritance from global for this sep., sysBase, sysMods.
          seeds <- strsplit(crow$PoolSeeds,
                            split = thisSeparator,
                            fixed = TRUE)[[1]]

          lapply(seeds, function(s) {
            RMTRCode2::LawMorton1996_species(
              Basal = b$species[1],
              Consumer = b$species[2],
              Parameters = b$LM1996ParamSet,
              LogBodySize = b$LM1996BodySize +
                m$PoolBodySizes[crow$Framework, ],
              seed = as.numeric(s)
            )
          })
        }, cs = cases, b = systemBase, m = systemMods
      )
    })

    return(env)
  }
)

# Generate Matrices: ##########################################################
theEnvironments <- lapply(
  theEnvironments, function(env) {
    with(env, {
      # matrs will be a row-associated list of pool-associated lists of
      # seed-associated matrices.
      matrs <- lapply(
        1:nrow(cases), function(row, b, m, cs, p) {
          crow <- cases[row, ]
          poolGrp <- pools[[row]]

          # Using inheritance from global for this sep., sysBase, sysMods.
          seeds <- strsplit(crow$EnvironmentSeeds,
                            split = thisSeparator,
                            fixed = TRUE)[[1]]

          # Group using number of environments per system.
          seedGrps <- lapply(
            1:b$systemsPerParamSet, function(i, s) {
            s[1:b$environsPerSystem
                  + (i - 1) * b$environsPerSystem]
          }, s = seeds)

          lapply(
            seq_along(seedGrps), function(i, s, pl) {
              RMTRCode2::CreateEnvironmentInteractions(
              Pool = pl[[i]],
              NumEnvironments = b$environsPerSystem,
              ComputeInteractionMatrix = RMTRCode2::LawMorton1996_CommunityMat,
              Parameters = b$LM1996ParamSet * c(1, 1, 1, 1, 1,
                                                m$InteractionsNoiseMultiplier[
                                                  crow$Framework
                                                ]),
              EnvironmentSeeds = as.numeric(s[[i]])
            )
          }, s = seedGrps, pl = poolGrp)
        }, cs = cases, b = systemBase, m = systemMods, p = pools
      )
    })

    return(env)
  }
)

# Generate Events: ############################################################
theEnvironments <- lapply(
  theEnvironments, function(env) {
    with(env, {
      # evnts will be a row-associated list of pool-associated lists of
      # seed-associated event data frames.
      evnts <- lapply(
        1:nrow(cases), function(row, b, m, cs, p) {
          crow <- cases[row, ]
          poolGrp <- pools[[row]]
          matrGrp <- matrs[[row]]

          # Using inheritance from global for this sep., sysBase, sysMods.
          seeds <- strsplit(crow$HistorySeeds,
                            split = thisSeparator,
                            fixed = TRUE)[[1]]

          # Group using number of environments per system.
          seedGrps <- lapply(
            1:b$systemsPerParamSet, function(i, s) {
              s[1:b$historiesPerSystem
                    + (i - 1) * b$historiesPerSystem]
            }, s = seeds)

          # We need to identify the appropriate rate to couple the events to.
          rateGrp <- lapply(
            matrGrp,
            function(Mat) {
              # Maximum magnitude eigenvalue should correspond to the (timescale
              # of the) fluctuations caused by the dominant eigenvector.
              max(unlist(lapply(
                Mat$Mats, function(m) abs(eigen(m, only.values = TRUE)$values)
                )))
            }
          )


          lapply(
            seq_along(seedGrps), function(i, s, pl, r) {
              nrowpl <- nrow(pl[[i]])
              CreateAssemblySequence(
                Species = nrowpl,
                NumEnvironments = b$environsPerSystem,
                ArrivalEvents = b$eventNumberFunc(
                  b$environsPerSystem, spec = nrowpl, const = 7
                  # More than 5 since we are not seeing enough mixing.
                ),
                ArrivalRate = r[[i]],
                ArrivalFUN = ArrivalFUN_Example,
                ExtinctEvents = b$eventNumberFunc(
                  b$environsPerSystem, spec = nrowpl, const = 7
                  # More than 5 since we are not seeing enough mixing.
                ),
                ExtinctRate = r[[i]],
                ExtinctFUN = ExtinctFUN_Example,
                HistorySeed = s[[i]]
              )
            }, s = seedGrps, pl = poolGrp, r = rateGrp)
        }, cs = cases, b = systemBase, m = systemMods, p = pools
      )
    })

    return(env)
  }
)
