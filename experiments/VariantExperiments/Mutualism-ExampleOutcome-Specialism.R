# Try to detect specialism vs generalism.
# Easiest way to do this might be to consider the evenness:
#   of intra-guild losses of a node,
#   of inter-guild gains for a node.
# We'll take the relative Shannon's Entropy (Pielou's index)
#   -Sum_i x_i log (x_i) / log(N)
#   where i is species index, x_i is relative effect and N is number of effects.
# Counterargument: this doesn't work for comparing between a species that feeds
# evenly on 2 species but no others and one that feeds mostly on 1, but rarely
# on 100 others.
# "The Fix": our networks aren't sparse; indeed, they are complete within patch.
# As such, the counterargument isn't a problem. Each species has a minimal
# contribution from every other species, so we're comparing apples:apples.
# Common set-up: ##############################################################
library(rlang)
library(reshape2)
library(dplyr)        # Data Manipulation
library(Matrix)       # Sparse Matrices
library(foreach)
library(iterators)
library(tidyr)        # Data Pivotting
library(ggplot2)      # 2-D Plot
library(gridExtra)
library(RMTRCode2)

source("Mutualism.R")

### Parameters: ###############################################################

# Problems with X11
options(bitmapType = "cairo")

by_for_thinning <- 10 # time steps
divide_time_by <- 1E1 # time units
burn_in <- 1E1 # time units

preferred_rows_per_event <- 1.5

burnout <- 60 # In original units.
numEnvs <- 10

### Load: #####################################################################
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

##### Careful to process as we load: ##########################################

results_specialism <- foreach::foreach(
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

    # PerIslandDistance <- theParameters$Space
    # PerIslandDistance <- 10^(as.numeric(sub(pattern = "_", replacement = "-",
    #                                         x = PerIslandDistance, fixed = TRUE)))
    #
    # # Ring Dynamics
    # DistanceMatrix <- Matrix::bandSparse(
    #   numEnvs, k = c(-1, 1),
    #   diagonals = list(rep(PerIslandDistance, numEnvs - 1),
    #                    rep(PerIslandDistance, numEnvs - 1))
    # )
    # DistanceMatrix[numEnvs, 1] <- PerIslandDistance
    # DistanceMatrix[1, numEnvs] <- PerIslandDistance
    #
    # DispersalMatrices <- RMTRCode2::CreateDispersalMatrix(
    #   EnvironmentDistances = DistanceMatrix,
    #   SpeciesSpeeds = rep(1, nrow(PoolsMats[[1]]$Pool))
    # )

    GeneralismByGuildPatch <-
      apply(fileContents$Abundance[seq(
        from = 1,
        to = nrow(fileContents$Abundance),
        by = round(nrow(fileContents$Abundance) / 100)
      ),],
      MARGIN = 1, function(row) {
        # Row = time step + abundance
        interactions <- t(apply(
          rlang::fn_env(PerCapitaDynamics)$InteractionMatrix,
          MARGIN = 1,
          function(interactionsReceived) interactionsReceived * row[-1]
        )) # Transpose to have rows = species on patch, columns = contributions.

        # if (sum(interactions != 0) > 5) {
        #   print("hi")
        # }

        # Multiplication is safe so long as we consider each guild-action type
        # to be different.
        # Interactions is a matrix: row is recipient, column is contribution.

        # We trust that st has a breakdown of the interaction types.
        # Something different is needed for pred-prey.
        interactions <- abs(interactions)

        # We need to convert into percentage contribution for each recipient and
        # interaction type (and patch, but this is implicit in the structure).
        # Note we are not worrying about dispersal.

        st <- rlang::fn_env(PerCapitaDynamics)$st

        Generalism <- matrix(NA, nrow = length(row) - 1, ncol = length(st) - 1)
        # Rows: Species-Patch, Columns: Guild-Patch.

        for (i in 2:length(st)) {
          lo <- st[i - 1] + 1; hi <- st[i]
          contributionsTotal <- rowSums(interactions[, lo:hi])
          interactions[, lo:hi] <- interactions[, lo:hi] / contributionsTotal
          # Reminder: R is column indexed; each row is divided by same value.

          # Note: those receiving nothing will be 0 / 0. We'll fix at the end.

          # Note: for Entropy calcs, as x -> 0, x log x -> 0.
          interactions[, lo:hi] <- ifelse(
            interactions[, lo:hi] == 0, 0,
            - interactions[, lo:hi] * log(interactions[, lo:hi]) / log(hi - lo)
          )

          Generalism[, i - 1] <- rowSums(interactions[, lo:hi])
          Generalism[contributionsTotal == 0, i - 1] <- NA
        }

        # Can't have generalism/specialism if not present.
        Generalism[row[-1] == 0, ] <- NA

        retval <- list(
          GeneralismScores = Generalism,
          Time = row[1]
        )

        return(retval)
      })

    return(GeneralismByGuildPatch)
  })
}

results_specialismDF <- lapply(
  results_specialism,
  function(fileresults) {
    dplyr::bind_rows(lapply(fileresults, function(timeresults) {
      # Desired output:
      # data.frame(
      #   Time = ,
      #   Species = ,
      #   SpeciesEnvironment =,
      #   Guild = ,
      #   GuildEnvironment =,
      #   Generalism =
      # )
      temp <- timeresults$GeneralismScores %>% melt %>% na.omit
      temp <- temp %>% dplyr::mutate(
        Time = timeresults$Time,
        Species =
          ((Var1 - 1) %% (sum(numSpecies))) + 1,
        SpeciesEnvironment =
          ((Var1 - 1) %/% (sum(numSpecies))) + 1,
        Guild =
          ((Var2 - 1) %% (length(numSpecies))) + 1,
        GuildEnvironment =
          ((Var2 - 1) %/% (length(numSpecies))) + 1
      ) %>% dplyr::rename(
        Generalism = value
      ) %>% dplyr::select(
        -Var1, -Var2
      )
    }))
  })

results_specialismDF <- dplyr::bind_rows(
  lapply(seq_along(results_specialismDF), function(i, df, properties) {
    cbind(df[[i]], properties[i, ])
  }, df = results_specialismDF, properties = Properties)
)

ggplot2::ggplot(
  results_specialismDF %>% dplyr::mutate(
    Space = gsub(pattern = "_", replacement = "-", x = Space, fixed = TRUE),
    Space = 10^as.numeric(Space),
    Dispersal = 1 - exp( -2 / as.numeric(Space) )
  ), ggplot2::aes(
    x = Time, y = Generalism,
    group = interaction(Species, SpeciesEnvironment, Guild, GuildEnvironment),
    color = Species, linetype = factor(Guild)
  )
) + ggplot2::geom_line() + ggplot2::facet_wrap(Guild ~ Dispersal, nrow = 2)
