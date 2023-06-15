# Borrowing from Viking_HandleOutput_Invadability2.R

library(dplyr)        # Data Manipulation
library(Matrix)       # Sparse Matrices
library(foreach)
library(iterators)
library(tidyr)        # Data Pivotting
library(ggplot2)      # 2-D Plot
library(gridExtra)

source("Mutualism.R")

# Parameters: #################################################################

# Problems with X11
options(bitmapType = "cairo")

by_for_thinning <- 10 # time steps
divide_time_by <- 1E1 # time units
burn_in <- 1E1 # time units

preferred_rows_per_event <- 1.5

burnout <- 60 # In original units.
numEnvs <- 10

# Load: #######################################################################
# All .RData
files_dat <- dir(
  path = "Data_2023-05-19", # "Data_2022-09-16",
  pattern = "Mutualism-Example.+[.]RData$", # "MNA[-]ExampleOutcome[-].+[.]RData$",
  full.names = TRUE
)

# Separate out PoolMats
files_dat_PM <- files_dat[
  grepl(x = files_dat,
        pattern = "PoolMats",
        fixed = TRUE)
  ]
files_dat <- files_dat[
  !grepl(x = files_dat,
         pattern = "PoolMats",
         fixed = TRUE)
  ]

PoolsMats <- sapply(
  files_dat_PM,
  load_safe,
  simplify = FALSE, USE.NAMES = TRUE
)

numSpecies <- table(PoolsMats[[1]]$Pool$Type)
# Put Producers first.
numSpecies <- numSpecies[c(length(numSpecies), 2:(length(numSpecies)) - 1)]

IntMat <- Matrix::bdiag(
  PoolsMats[[1]]$InteractionMatrices$Mats
)
PerCapitaDynamics <- PerCapitaDynamics_Mutualistic1(
  PoolsMats[[1]]$Pool$ReproductionRate, IntMat,
  NumEnvironments = numEnvs,
  SpeciesTypes = numSpecies, Saturations = PoolsMats[[1]]$Saturation
)

# For each file, we have a different dispersal matrix.
Properties <- do.call(
  rbind, strsplit(basename(files_dat), split = "-", fixed = TRUE)
)[, -c(1:4)]

colnames(Properties) <- c("Arrangement", "Space", "Imm", "Ext", "Proportion")
Properties <- data.frame(Properties)

results_Invadability <- foreach::foreach(
  file = iterators::iter(files_dat),
  theParameters = iterators::iter(Properties, by = "row"),
  .packages = c("dplyr")
) %do%{
  print(file)

  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    fileContents <- get(fileContents)

    print("loaded")

    PerIslandDistance <- theParameters$Space
    PerIslandDistance <- 10^(as.numeric(sub(pattern = "_", replacement = "-",
                                            x = PerIslandDistance, fixed = TRUE)))

    # Ring Dynamics
    DistanceMatrix <- Matrix::bandSparse(
      numEnvs, k = c(-1, 1),
      diagonals = list(rep(PerIslandDistance, numEnvs - 1),
                       rep(PerIslandDistance, numEnvs - 1))
    )
    DistanceMatrix[numEnvs, 1] <- PerIslandDistance
    DistanceMatrix[1, numEnvs] <- PerIslandDistance

    DispersalMatrices <- RMTRCode2::CreateDispersalMatrix(
      EnvironmentDistances = DistanceMatrix,
      SpeciesSpeeds = rep(1, nrow(PoolsMats[[1]]$Pool))
    )

    # fileContents is a result
    Invadability <- list(
      Invadabilities = list(thinAndCalculateInvadabilities(
        fileContents,
        dyn = PerCapitaDynamics,
        dis = DispersalMatrices,
        preferred_rows_per_event = preferred_rows_per_event
      ))
    )

    return(Invadability)
  })
}

results_TimeBeta <- foreach::foreach(
  file = iterators::iter(files_dat),
  theParameters = iterators::iter(Properties, by = "row"),
  .packages = c("dplyr")
) %do%{
  print(file)

  tryCatch({
    # 1 Load the Data: ########################################################
    # Files contain a single large list object
    fileContents <- load(file)

    # Move "pointer" to the large list object
    fileContents <- get(fileContents)

    print("loaded")

    TimeJaccard <- list(
      Betas = list(Calculate_TimeJaccard(
        fileContents,
        nspecies = numSpecies,
        minTime = max(fileContents$Events$Times)/101
      ))
    )

    TimeJaccard$Betas <- lapply(TimeJaccard$Betas, function(b) {
      b %>% dplyr::mutate(
        Time  = Time / divide_time_by
      )
    })

    return(TimeJaccard)
  })
}

names(results_Invadability) <- basename(files_dat)
names(results_TimeBeta) <- basename(files_dat)

# Saved to Invasibility_TemporalJac_2023-05-19.RData.

# Plotting: ###################################################################
### Temporal Beta: ############################################################
##### Prep Work: ##############################################################
results_TimeBeta <- lapply(results_TimeBeta, function(x) x[[1]][[1]])
results_TimeBeta <- dplyr::bind_rows(results_TimeBeta, .id = "File")

results_TimeBeta <- results_TimeBeta %>% tidyr::separate(
  sep = "[-]", col = File, into = c(
    "Mutualism", "ExampleExtProp", "Result", "EnvNumber", "Gridtype",
    "Space", "Immigration", "Extirpation", "ExtirpationProportion"
  )
) %>% dplyr::select(
  -Mutualism, -ExampleExtProp, -Result
) %>% dplyr::mutate(
  EnvNumber = as.numeric(
    gsub(pattern = "[a-zA-Z]", replacement = "", EnvNumber)
  ),
  Space = 10^as.numeric(
    gsub(pattern = "[_]", replacement = "-", Space)
  ),
  ExtirpationProportion = as.numeric(
    gsub(pattern = "(ExtProp)|(.RData)", replacement = "",
         ExtirpationProportion)
  )
)

spaceToDispersal <- data.frame(
  Space = unique(results_TimeBeta$Space)
) %>% dplyr::mutate(
  Dispersal = 1 - exp( -2 / as.numeric(Space) ),
  Dispersal = paste0(
    formatC(Dispersal))
) %>% dplyr::arrange(
  as.numeric(Dispersal)
) %>% dplyr::mutate(
  # Edit to change the order of the factor levels
  Dispersal = factor(Dispersal, ordered = TRUE,
                     levels = Dispersal)#,
  #Space = factor(Space, ordered = TRUE, levels = rev(Space))
)

results_TimeBeta <- results_TimeBeta %>% dplyr::left_join(
  spaceToDispersal, by = "Space"
)

##### Plot: ###################################################################
ggplot2::ggplot(
  results_TimeBeta %>% dplyr::filter(
    Time > burn_in / divide_time_by, Time < burnout / divide_time_by
  ),
  ggplot2::aes(
    x = Dispersal,
    y = Value,
    fill = Measurement
  )
) + ggplot2::geom_boxplot(
)

### Invasibility: #############################################################
##### Prep Work: ##############################################################
# Echoing Viking_HandleInvadability_Combine.R
results_Invadability <- lapply(
  results_Invadability,
  function(x) {
    target <- x$Invadabilities[[1]]$invadability # only 1 run.
    numSpec <- ncol(target) / numEnvs

    data.frame(
      Environment = 1:numEnvs,
      Local = unlist(lapply(1:numEnvs, function(j) {
        sum(target[1, 1:numSpec + numSpec * (j - 1)])
      }))
    )
  }
)
results_Invadability <- dplyr::bind_rows(results_Invadability, .id = "File")

results_Invadability <- results_Invadability %>% tidyr::separate(
  sep = "[-]", col = File, into = c(
    "Mutualism", "ExampleExtProp", "Result", "EnvNumber", "Gridtype",
    "Space", "Immigration", "Extirpation", "ExtirpationProportion"
  )
) %>% dplyr::select(
  -Mutualism, -ExampleExtProp, -Result
) %>% dplyr::mutate(
  EnvNumber = as.numeric(
    gsub(pattern = "[a-zA-Z]", replacement = "", EnvNumber)
  ),
  Space = 10^as.numeric(
    gsub(pattern = "[_]", replacement = "-", Space)
  ),
  ExtirpationProportion = as.numeric(
    gsub(pattern = "(ExtProp)|(.RData)", replacement = "",
         ExtirpationProportion)
  )
)

results_Invadability <- results_Invadability %>% dplyr::left_join(
  spaceToDispersal, by = "Space" # From TimeBeta workings above
)

##### Plot: ###################################################################

ggplot2::ggplot(
  results_Invadability,
  ggplot2::aes(
    x = Dispersal,
    y = Local,
    color = factor(Environment, ordered = TRUE, levels = 1:10),
    group = Environment
  )
) + ggplot2::geom_point(
  alpha = 0.5
) + ggplot2::scale_color_viridis_d(
  name = "Environment"
) + ggplot2::labs(
  y = "Local Species Invasibility"
) + ggplot2::geom_line(
  alpha = 0.3
) + ggplot2::geom_point(
  data = results_Invadability %>% dplyr::group_by(
    Dispersal
  ) %>% dplyr::summarise(
    Local = median(Local)
  ),
  color = "red",
  shape = 4, size = 4,
  group = "Median"
) + ggplot2::geom_point(
  data = results_Invadability %>% dplyr::group_by(
    Dispersal
  ) %>% dplyr::summarise(
    Local = mean(Local)
  ),
  color = "red",
  shape = 5, size = 4,
  group = "Mean"
) + ggplot2::geom_text(
  data = data.frame(
    Local = c(6.25, 18.75),
    Dispersal = 0.75,
    Text = c("Median", "Mean")
  ),
  mapping = ggplot2::aes(
    label = Text
  ),
  group = "Text", color = "black"
)
