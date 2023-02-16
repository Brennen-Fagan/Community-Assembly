# Abstract: ####################################################################
# This file takes the output of "MNA-MasterAttempt-GenerateCases.R" and converts
# it to a set of pools, matrices and events to prevent bottlenecking during a
# submission to Viking.

library(RMTRCode2)
library(parallel)
library(foreach)
library(doParallel)
library(iterators)
library(dplyr)
library(tidyr)

cores <- 14 #parallel::detectCores() - 1
clust <- parallel::makeCluster(cores)
doParallel::registerDoParallel(clust)

# https://stackoverflow.com/a/15373917
thisFile <- function() {
  # cmdArgs <- commandArgs(trailingOnly = FALSE)
  # needle <- "--file="
  # match <- grep(needle, cmdArgs)
  # if (length(match) > 0) {
  #   # Rscript
  #   return(normalizePath(sub(needle, "", cmdArgs[match])))
  # } else {
  #   # 'source'd via R console
  #   return(normalizePath(sys.frames()[[1]]$ofile))
  # }
  return(".")
}

# Initialisation: ##############################################################
print("Initialisation.")
# Identify the files to take.
thisDirectory <- "."#dirname(thisFile())
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
print("Environments.")
theEnvironments <- lapply(theCases, function(case) {
  tempenv <- new.env()
  assign("name", basename(case), envir = tempenv)
  assign("cases", read.csv(case, as.is = TRUE), envir = tempenv)
  return(tempenv)
})

# Generate Pools: #############################################################
### Rearrange: ################################################################
# We rearrange the cases as we go within the environments so that
# each row links a pool, set of matrices, and history,
# as opposed to each row collecting all pools, matrices,
# and histories for a given set of parameters.
print("Rearrange for pools.")
theEnvironments <- lapply(theEnvironments, function(env) {

  assign("pools", vector("list", nrow(env$cases)), envir = env)
  assign("matrs", vector("list", nrow(env$cases)), envir = env)
  assign("evnts", vector("list", nrow(env$cases)), envir = env)

  with(env, {
    # systemBase$systemsPerParamSet systems
    # with 1 pool per system,
    # systemBase$environsPerSystem environments per system
    # with 1 matrix per environment
    # and systemBase$historiesPerSystem histories per system.

    # The environment seeds.
    envrseeds <- matrix(unlist(strsplit(
      cases$EnvironmentSeeds, split = thisSeparator, fixed = TRUE
    )),
    nrow = systemBase$environsPerSystem,
    # ncol = systemBase$systemsPerParamSet * length(params)
    # fills by column = by system x by parameter set
    )

    # The environment seeds.
    histseeds <- matrix(unlist(strsplit(
      cases$HistorySeeds, split = thisSeparator, fixed = TRUE
    )),
    nrow = systemBase$historiesPerSystem,
    # ncol = systemBase$systemsPerParamSet * length(params)
    # fills by column = by system x by parameter set
    )

    # The parameters for each system and pool seeds
    cases <- cases %>% dplyr::mutate(
      Parameters = 1:nrow(cases)
    ) %>% tidyr::separate(
      # Break System-Pool seeds into columns
      PoolSeeds,
      into = paste0("P", as.character(1:systemBase$systemsPerParamSet)),
      sep = paste0("[", thisSeparator, "]"),
    ) %>% tidyr::pivot_longer(
      # Columns to Rows
      cols = paste0("P", as.character(1:systemBase$systemsPerParamSet)),
      names_to = "System", values_to = "PoolSeed"
    ) %>% dplyr::mutate(
      # Colname to index
      System = as.numeric(substring(System, 2)),
      PoolSeed = as.numeric(PoolSeed)
    )
  })

  return(env)
})

### Calculate: ################################################################
print("Generate Pools.")
print(Sys.time())
print(system.time(
  foreach::foreach(
    env = iterators::iter(theEnvironments)
  ) %:% foreach::foreach(
    crow = iterators::iter(with(env, cases), by = "row")
  ) %do% { # Faster without the overhead.
    with(env,
         pools[[crow$Parameters]][[crow$System]] <-
           RMTRCode2::LawMorton1996_species(
             Basal = systemBase$species[1],
             Consumer = systemBase$species[2],
             Parameters = systemBase$LM1996ParamSet,
             LogBodySize = systemBase$LM1996BodySize +
               systemMods$PoolBodySizes[crow$Framework, ],
             seed = crow$PoolSeed
           )
    )
  }))

# Generate Matrices: ##########################################################
### Rearrange: ################################################################
print("Rearrange for matrices. ")
theEnvironments <- lapply(theEnvironments, function(env) {

  with(env, {
    # systemBase$systemsPerParamSet systems
    # systemBase$environsPerSystem environments per system
    # with 1 matrix per environment
    cases$EnvironmentSeeds <- apply(envrseeds, MARGIN = 2, toString)
  })

  return(env)
})

### Calculate: ################################################################
print("Generate Matrices.")
print(Sys.time())
print(system.time(
  matrices <- foreach::foreach(
    envIndex = iterators::iter(
      rep(1:length(theEnvironments),
          times = unlist(lapply(theEnvironments,
                                function(env) {nrow(env$cases)}))
      )
    ),
    crow = iterators::iter(
      dplyr::bind_rows(lapply(theEnvironments,
                              function(env) env$cases)),
      by = "row"
    )
  ) %dopar% {
    # print(envIndex)
    env <- theEnvironments[[envIndex]]

    templist <- list()

    # env$matrs[[crow$Parameters]][[crow$System]] <-
    templist[[paste(envIndex, crow$Parameters, crow$System)]] <-
      RMTRCode2::CreateEnvironmentInteractions(
        Pool = env$pools[[crow$Parameters]][[crow$System]],
        NumEnvironments = systemBase$environsPerSystem,
        ComputeInteractionMatrix = RMTRCode2::LawMorton1996_CommunityMat,
        Parameters = systemBase$LM1996ParamSet * c(
          1, 1, 1, 1, 1,
          systemMods$InteractionsNoiseMultiplier[
            crow$Framework
          ]),
        EnvironmentSeeds = as.numeric(
          strsplit(crow$EnvironmentSeeds, split = ", ", fixed = TRUE)[[1]]
        )
      )

    return(templist)
  }))
matrices <- matrices[order(unlist(lapply(matrices, names)))]

for (mat in matrices) {
  indices <- as.numeric(strsplit(names(mat), split = " ", fixed = TRUE)[[1]])
  theEnvironments[[indices[1]]]$matrs[[indices[2]]][[indices[3]]] <- mat[[1]]
}

# Generate Events: ############################################################
### Rearrange: ################################################################
print("Rearrange for events. ")
theEnvironments <- lapply(theEnvironments, function(env) {

  with(env, {
    # systemBase$systemsPerParamSet systems
    # systemBase$environsPerSystem environments per system
    # with 1 matrix per environment
    cases$HistorySeeds <- apply(histseeds, MARGIN = 2, toString)
  })

  return(env)
})

### Calculate: ################################################################
print("Generate Events.")
print(Sys.time())
print(system.time(
  events <- foreach::foreach(
    envIndex = iterators::iter(
      rep(1:length(theEnvironments),
          times = unlist(lapply(theEnvironments,
                                function(env) {nrow(env$cases)}))
      )
    ),
    crow = iterators::iter(
      dplyr::bind_rows(lapply(theEnvironments,
                              function(env) env$cases)),
      by = "row"
    )
  ) %dopar% {
    # print(envIndex)
    env <- theEnvironments[[envIndex]]

    # We need to identify the appropriate rate to couple the events to.
    rate <-
        # Maximum magnitude eigenvalue should correspond to the (timescale
        # of the) fluctuations caused by the dominant eigenvector.
        max(unlist(lapply(
          env$matrs[[crow$Parameters]][[crow$System]]$Mats,
          function(m) abs(eigen(m, only.values = TRUE)$values)
        )))

    nrowpl <- nrow(env$pools[[crow$Parameters]][[crow$System]])

    templist <- list()

    # env$matrs[[crow$Parameters]][[crow$System]] <-
    templist[[paste(envIndex, crow$Parameters, crow$System)]] <-
      lapply(strsplit(crow$HistorySeeds, split = ", ", fixed = TRUE)[[1]],
             function(s) {
               RMTRCode2::CreateAssemblySequence(
                 Species = nrowpl,
                 NumEnvironments = systemBase$environsPerSystem,
                 ArrivalEvents = systemBase$eventNumberFunc(
                   systemBase$environsPerSystem, Spec = nrowpl, Const = 7
                   # More than 5 since we are not seeing enough mixing.
                 ) * systemMods$NeutralRateMultipliers[
                   crow$Neutral, 1
                 ],
                 ArrivalRate = rate * systemMods$NeutralRateMultipliers[
                   crow$Neutral, 1
                 ],
                 ArrivalFUN = RMTRCode2::ArrivalFUN_Example2,
                 ExtinctEvents =  systemBase$eventNumberFunc(
                   systemBase$environsPerSystem, Spec = nrowpl, Const = 7
                   # More than 5 since we are not seeing enough mixing.
                 ) * systemMods$NeutralRateMultipliers[
                   crow$Neutral, 2
                 ],
                 ExtinctRate = rate * systemMods$NeutralRateMultipliers[
                   crow$Neutral, 2
                 ],
                 ExtinctFUN = RMTRCode2::ExtinctFUN_Example2,
                 HistorySeed = as.numeric(s)
               )
             }
      )

    return(templist)
  }))
events <- events[order(unlist(lapply(events, names)))]

for (evt in events) {
  indices <- as.numeric(strsplit(names(evt), split = " ", fixed = TRUE)[[1]])
  theEnvironments[[indices[1]]]$evnts[[indices[2]]][[indices[3]]] <- evt[[1]]
}

# Cleanup: ####################################################################
parallel::stopCluster(clust)

lapply(theEnvironments,
       function(env) {
         newname <- strsplit(env$name, split = ".", fixed = TRUE)[[1]]
         newname <- paste0(newname[-length(newname)], "-Prepared.RData")
         save(file = file.path(thisDirectory, newname),
              envir = env, cases, evnts, matrs, pools)
       }
)
