# Extracted from "FirstAttempt-Doc-Analysis2-Gallery3.Rmd".
# Then modified to comply with the standards in
# "Viking_HandleDiversity_ParametersAndPlots3.R".
# This should result in a 2x2 figure.
# Row 1: Species Presence Absence, by Size.
# Row 2: Diversity Metrics over Time.

library(dplyr)        # Data Manipulation
library(tidyr)        # Data Pivotting
library(ggplot2)      # 2-D Plot
library(gridExtra)

# Problems with X11
options(bitmapType = "cairo")

by_for_thinning <- 100 # time steps
divide_time_by <- 1E4 # time units
burn_in <- 1E4 # time units

# Functions: ###################################################################

load_safe_thin <- function(fname, bythin, divtime, burn) {
  loaded <- tryCatch({load(fname)},
                     error = function(e) {
                       print(fname)
                       print(e)
                       return(NA)
                     })
  if (is.na(loaded)) {
    return(NA)
  } else {
    loaded <- get(loaded)
    loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                             to = nrow(loaded$Abundance),
                                             by = bythin), ]

    toEliminate <-
      loaded$Abundance[, -1] <
      loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
    loaded$Abundance[, -1][toEliminate] <- 0

    loaded$Abundance <- loaded$Abundance[
      loaded$Abundance[, 1] > burn,
      ]

    loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divtime
    return(loaded)
  }
}

load_safe <- function(fname) {
  loaded <- tryCatch({load(fname)},
                     error = function(e) {
                       print(fname)
                       print(e)
                       return(NA)
                     })
  if (all(is.na(loaded))) {
    return(NA)
  } else {
    return(sapply(loaded, get,
                  envir = sys.frame(sys.parent(0)),
                  simplify = FALSE, USE.NAMES = TRUE))
  }
}

Calculate_Diversity <- function(loaded, nspecies) {
  loaded$Abundance[, -1] <-
    loaded$Abundance[, -1] > loaded$Parameters$EliminationThreshold

  ### Alpha Diversity: ##################################################
  diversity_alpha <- lapply(
    1:loaded$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      env_basal <- env[, 1:nspecies[1]]
      env_consumer <- env[, nspecies[1] + 1:nspecies[2]]
      richness <- rowSums(env)
      richness_basal <- rowSums(env_basal)
      richness_consumer <- rowSums(env_consumer)
      species <- apply(
        env, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_basal <- apply(
        env_basal, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_consumer <- apply(
        env_consumer, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0) + nspecies[1]) # + basals doesn't change it, but
          # results in a more intuitive numbering system of species.
        }
      )
      data.frame(Time = time,
                 Richness = richness,
                 Richness_Basal = richness_basal,
                 Richness_Consumer = richness_consumer,
                 Species = species,
                 Species_Basal = species_basal,
                 Species_Consumer = species_consumer,
                 Environment = i,
                 stringsAsFactors = FALSE)
    },
    abund = loaded$Abundance,
    numSpecies = sum(nspecies)
  )

  diversity_alpha <- dplyr::bind_rows(diversity_alpha)

  print("alpha")
  ### Gamma Diversity: ##################################################
  diversity_gamma <- diversity_alpha %>% dplyr::group_by(
    Time
  ) %>% dplyr::summarise(
    Mean = mean(Richness),
    Mean_Basal = mean(Richness_Basal),
    Mean_Consumer = mean(Richness_Consumer),
    Var = var(Richness),
    Var_Basal = var(Richness_Basal),
    Var_Consumer = var(Richness_Consumer),

    SpeciesTotal = toString(sort(unique(unlist(strsplit(paste(
      Species, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Basal = toString(sort(unique(unlist(strsplit(paste(
      Species_Basal, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Consumer = toString(sort(unique(unlist(strsplit(paste(
      Species_Consumer, collapse = ", "), split = ", ", fixed = TRUE))))),

    Richness = unlist(lapply(
      strsplit(
        SpeciesTotal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Richness_Basal = unlist(lapply(
      strsplit(
        SpeciesTotal_Basal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Richness_Consumer = unlist(lapply(
      strsplit(
        SpeciesTotal_Consumer, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    ))
  ) %>% dplyr::select(
    -dplyr::starts_with("Species")
  ) %>% tidyr::pivot_longer(
    cols = !Time,
    names_to = "Measurement",
    values_to = "Value"
  ) %>% dplyr::mutate(
    Environment = "Gamma"
  )

  print("gamma")
  ### Beta Diversity (Jaccard, Space): ##################################
  diversity_beta <- apply(
    loaded$Abundance,
    MARGIN = 1, # Rows
    function(row, envs) {
      time <- row[1]
      # Vegan complains about rows with all 0's.
      # The warning is generic, so we cannot silence it specifically.
      dists <- suppressWarnings(vegan::vegdist(
        method = "jaccard",
        x = matrix(row[-1] > 0, nrow = envs, byrow = TRUE)
      ))

      dataf <- expand.grid(
        Env1 = 1:envs,
        Env2 = 1:envs
      ) %>% dplyr::filter(
        Env1 < Env2
      ) %>% dplyr::mutate(
        Time = time,
        Jaccard = as.numeric(dists)
      )

      return(dataf)
    },
    envs = loaded$NumEnvironments
  )
  diversity_beta <- dplyr::bind_rows(
    diversity_beta
  ) %>% tidyr::unite(
    "Environment", Env1, Env2, sep = " "
  ) %>% dplyr::rename(
    Value = Jaccard
  ) %>% dplyr::mutate(
    Measurement = "Jaccard"
  )

  diversity_beta_avg <- diversity_beta %>% dplyr::group_by(
    Time, Measurement
  ) %>% dplyr::summarise(
    Value = mean(Value)
  ) %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    Environment = "Mean"
  )
  print("beta")

  Diversities <- dplyr::bind_rows(
    diversity_alpha %>% dplyr::select(
      -Species_Basal, -Species_Consumer, -Species
    ) %>% tidyr::pivot_longer(
      c(-Time, -Environment),
      names_to = "Measurement",
      values_to = "Value"
    ) %>% dplyr::mutate(
      Environment = as.character(Environment)
    ),
    diversity_beta,
    diversity_beta_avg,
    diversity_gamma
  )

  ### Return Diversities: ###############################################
  return(Diversities)
}

Calculate_Species <- function(result, bintimes = FALSE) {
  SpeciesPerEnvironment <- lapply(
    1:result$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      # Need to retrieve Position and Value
      species <- apply(
        cbind(time, env), MARGIN = 1,
        FUN = function(x) {
          time <- x[1]
          dat <- x[-1]
          if (any(dat > 0)) {
            positions <- (which(dat > 0))
            values <- dat[positions]
            data.frame(
              Time = time,
              Species = positions,
              Abundance = values,
              row.names = NULL
            )
          } else {NULL}
          # Returns as list
        }
      )
      return(
        dplyr::bind_rows(species) %>% dplyr::mutate(
          Environment = i
        )
      )
    },
    abund = result$Abundance,
    numSpecies = (ncol(result$Abundance) - 1) / result$NumEnvironments
  )

  if (bintimes) {
    # Should equalise time steps.
    SpeciesPerEnvironment <- lapply(
      SpeciesPerEnvironment, function(SPE) {
        SPE %>% dplyr::mutate(
          TimeFloor = floor(Time*10)/10
        ) %>% dplyr::group_by(
          TimeFloor, Species, Environment
        ) %>% dplyr::summarise(
          Abundance = median(Abundance, na.rm = TRUE)
        )
      })
  }

  return(dplyr::bind_rows(SpeciesPerEnvironment))
}

# https://stackoverflow.com/a/46079299
extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}

# Implementation: ##############################################################

# All .RData
files_dat <- dir(
  path = "Data_2022-10-12",
  pattern = "MNA[-]ExampleOutcome[-].+[.]RData$",
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

Results <- sapply(
  files_dat,
  load_safe_thin,
  bythin = by_for_thinning,
  divtime = divide_time_by,
  burn = burn_in,
  simplify = FALSE, USE.NAMES = TRUE
)

names(Results) <- basename(names(Results))

PoolsMats <- sapply(
  files_dat_PM,
  load_safe,
  simplify = FALSE, USE.NAMES = TRUE
)

numSpecies <- table(PoolsMats[[1]]$Pool$Type)

Diversity <- sapply(
  USE.NAMES = TRUE, simplify = FALSE,
  Results, function(result) {
    if (length(result) == 1 && is.na(result)) {
      # Problem case.
      return(NA)
    }
    # print(paste("Calculating", Sys.time()))

    # Calculate the diversity.
    # We will need to extract the system properties from
    # the file names which we carry through using sapply.
    return(Calculate_Diversity(result, numSpecies))
  }
)

SpeciesPresence <-  sapply(
  USE.NAMES = TRUE, simplify = FALSE,
  Results, function(result) {
    if (length(result) == 1 && is.na(result)) {
      # Problem case.
      return(NA)
    }
    # print(paste("Calculating", Sys.time()))

    # Calculate the diversity.
    # We will need to extract the system properties from
    # the file names which we carry through using sapply.
    return(Calculate_Species(result))
  }
)

Properties <- strsplit(names(SpeciesPresence), '-',
                       fixed = TRUE)
# 1st Chunk: Name, Discard
# 2nd Chunk: Iteration + Distance
# 3rd Chunk: Result or Extinction Rate (or Arrival?)
# 4th Chunk: Number of Environments
# 5th Chunk: Space Type + .RData
# Note the mix of Keyword and Location Structure (oops).
# Note also that this strsplit character is a bad decision and should be changed for next time. (D'oh.)
# (E.g. Dates DD-MM-YYYY, Decimals 1.35e-05.)
Properties <- data.frame(
  do.call(rbind, Properties),
  stringsAsFactors = FALSE
)
names(Properties)[1:8] <- c(
  "MNA", "ExampleOutcome", "Result",
  "EnvNum", "Space", "Dist",
  "ImmigrationRate", "ExtirpationRateAND.RData"
)

Properties$FullName <- names(SpeciesPresence)

# Capture the position between the text (first group)
# and the set of numbers (somehow without the +).
# The \\K resets so that we do not capture any text.
patternString <- "((?>[a-zA-Z]+)(?=[0-9eE]))\\K"

# Split strings. Some of the trick will be to introduce
# a character to make the separation around. We use "_".
Properties <- Properties %>% dplyr::mutate(
  ExtirpationRateAND.RData = gsub(pattern = patternString,
                     replacement = "_",
                     x = ExtirpationRateAND.RData, perl = TRUE),
  EnvNum = gsub(pattern = patternString,
                replacement = "_",
                x = EnvNum, perl = TRUE)
) %>% tidyr::separate(
  ExtirpationRateAND.RData, into = c("ExtirpationRate", ".RData"),
  sep = "[.](?=[R])"
) %>% tidyr::separate(
  EnvNum, into = c("Env", "Environments"),
  sep = "[_]"
) %>% dplyr::select(
  -MNA, -ExampleOutcome, -Result, -.RData, -Env
) %>% dplyr::mutate(
  Distance = 10^as.numeric(Dist),
  ImmigrationRate = paste0("Imm.: ", ImmigrationRate),
  ExtirpationRate = paste0("Ext.: ", ExtirpationRate)
)

Diversity <- lapply(1:length(Diversity),
                    function(i, df, nm) {
                      df[[i]] %>% dplyr::mutate(
                        Simulation = nm[i]
                      )
                    },
                    df = Diversity,
                    nm = names(Diversity))

SpeciesPresence <- lapply(1:length(SpeciesPresence),
                          function(i, df, nm) {
                            df[[i]] %>% dplyr::mutate(
                              Simulation = nm[i]
                            )
                          },
                          df = SpeciesPresence,
                          nm = names(SpeciesPresence))

Diversity <- dplyr::left_join(
  dplyr::bind_rows(Diversity),
  Properties,
  by = c("Simulation" = "FullName")
)

SpeciesPresence <- dplyr::left_join(
  dplyr::bind_rows(SpeciesPresence),
  Properties,
  by = c("Simulation" = "FullName")
)

SpeciesPresence$Sizes <-
  PoolsMats[[1]]$Pool$Size[
    SpeciesPresence$Species
    ]
PoolsMats[[1]]$Pool <-
  PoolsMats[[1]]$Pool %>% dplyr::arrange(
    Size
  ) %>% dplyr::mutate(
    SizeID = 1:100
  )
SpeciesPresence$SizeID <-
  PoolsMats[[1]]$Pool$SizeID[
    SpeciesPresence$Species
    ]

SpeciesPresence <- SpeciesPresence %>% dplyr::group_by(
  Environments, Space, Dist, Distance, ImmigrationRate, ExtirpationRate,
  Species, Time, SizeID
) %>% dplyr::summarise(
  Count = n(), .groups = "drop"
) %>% dplyr::mutate(
  # Distance = ifelse(
  #   Space == "Line" | Space == "Ring", 1, Inf
  # ),
  Space = Distance,
  Dispersal = 1 - exp( -9 / as.numeric(Space) ),
  Dispersal = paste0(
    formatC(Dispersal))
)

Diversity <- Diversity %>% dplyr::mutate(
  # Distance = ifelse(
  #   Space == "Line" | Space == "Ring", 1, Inf
  # ),
  Space = Distance,
  Dispersal = 1 - exp( -9 / as.numeric(Space) ),
  Dispersal = paste0(
    formatC(Dispersal))
)

DiversityRibbons <- Diversity %>% dplyr::filter(
  !(Environment %in% c("Mean", "Gamma")),
  Measurement == "Richness" | Measurement == "Jaccard"
) %>%  dplyr::group_by(
  Time, Environments,
  Space, Dispersal,
  ImmigrationRate, ExtirpationRate,
  Measurement

  # Pool, Noise, Neutral, Space
) %>% dplyr::summarise(
  Low = unlist(dplyr::across(dplyr::any_of("Value"),
                             .fns = ~ quantile(.x, p = 0.1, na.rm = TRUE))),
  High = unlist(dplyr::across(dplyr::any_of("Value"),
                              .fns = ~ quantile(.x, p = 0.9, na.rm = TRUE))),
  .groups = "drop"
)

# Plots: #######################################################################

PLOT_TL <- ggplot2::ggplot(
  SpeciesPresence %>% dplyr::filter(
    Space == Inf
  ) %>% dplyr::mutate(
    Dispersal = "No Dispersal"
  ),
  ggplot2::aes(x = Time, y = SizeID, color = Count)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1,
  limits = c(1, 10),
  name = "No. of\nPatches"
) + ggplot2::geom_hline(
  yintercept = 34.5, color = "red"
) + ggplot2::labs(
  y = "Species by Size",
  x = paste0("Time, ", divide_time_by, " units")
) + ggplot2::theme_bw(
) + ggplot2::theme(
  strip.background = ggplot2::element_rect(
    colour = "black",
    fill = "cyan"
  )#,
  #legend.position = "none",
  #axis.text.x = ggplot2::element_blank()
) + ggplot2::scale_y_continuous(
  limits = c(0, 100)
) + ggplot2::facet_grid(
  ImmigrationRate ~ ExtirpationRate
) + ggplot2::ggtitle("No Dispersal")

PLOT_TM <- ggplot2::ggplot(
  SpeciesPresence %>% dplyr::filter(
    Space == 1e+05
  ) %>% dplyr::mutate(
    Dispersal = "Med. Dispersal"
  ),
  ggplot2::aes(x = Time, y = SizeID, color = Count)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1,
  limits = c(1, 10),
  name = "No. of\nPatches"
) + ggplot2::geom_hline(
  yintercept = 34.5, color = "red"
) + ggplot2::labs(
  y = "Species by Size",
  x = paste0("Time, ", divide_time_by, " units")
) + ggplot2::theme_bw(
) + ggplot2::theme(
  strip.background = ggplot2::element_rect(
    colour = "black",
    fill = "plum1"
  )#,
  #legend.position = "none",
  #axis.text.x = ggplot2::element_blank()
) + ggplot2::scale_y_continuous(
  limits = c(0, 100)
) + ggplot2::facet_grid(
  ImmigrationRate ~ ExtirpationRate
) + ggplot2::ggtitle("Med. Dispersal")

PLOT_TR <- ggplot2::ggplot(
  SpeciesPresence %>% dplyr::filter(
    Space == 1
  ) %>% dplyr::mutate(
    Dispersal = "Full Dispersal"
  ),
  ggplot2::aes(x = Time, y = SizeID, color = Count)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1,
  limits = c(1, 10),
  name = "No. of\nPatches"
) + ggplot2::geom_hline(
  yintercept = 34.5, color = "red"
) + ggplot2::labs(
  y = "Species by Size",
  x = paste0("Time, ", divide_time_by, " units")
) + ggplot2::theme_bw(
) + ggplot2::theme(
  strip.background = ggplot2::element_rect(
    colour = "black",
    fill = "darkorange"
  )#,
  #axis.text.x = ggplot2::element_blank()
) + ggplot2::scale_y_continuous(
  limits = c(0, 100)
) + ggplot2::facet_grid(
  ImmigrationRate ~ ExtirpationRate
) + ggplot2::ggtitle("Full Dispersal")

#PLOT_T_Legend <- extractLegend(PLOT_TR)
#PLOT_TR <- PLOT_TR + ggplot2::theme(legend.position = "none")

pasteCustom <- function(x, y) {
  #paste0("(", x, ", ", y, ")")
  ifelse(is.infinite(y), "None",
         ifelse(y == 1, "Full", "Med."))
}
legend_bl_name <- "Dispersal"

PLOT_BL <- ggplot2::ggplot(
  Diversity %>% dplyr::filter(
    Measurement == "Richness",
    Environment != "Gamma",
    Space != "Ring",
    Environment != "Mean"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space),
    group = interaction(ImmigrationRate, ExtirpationRate,
                        Dispersal)
  )
) + ggplot2::geom_line(
  alpha = 0.4,
  mapping = ggplot2::aes(
    group = interaction(ImmigrationRate, ExtirpationRate,
                        Dispersal, Environment)
  )
) + ggplot2::geom_line(
  data = Diversity %>% dplyr::filter(
    Measurement == "Richness",
    Environment != "Gamma",
    Space != "Ring",
    Environment == "Mean"
  ),
  size = 1.5
) + ggplot2::geom_ribbon(
  data = DiversityRibbons %>% dplyr::filter(
    Measurement == "Richness"
  ),
  mapping = ggplot2::aes(
    ymin = Low,
    ymax = High,
    x = Time,
    fill = pasteCustom(Dispersal, Space)
  ),
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::theme_bw(
) + ggplot2::labs(
  y = "Local Richness", # Number of Species",
  x = paste0("Time, ", divide_time_by, " units")
  #x = ""
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::scale_fill_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4")
) + ggplot2::facet_grid(
  ImmigrationRate ~ ExtirpationRate
) + ggplot2::ggtitle("Local Richness")

PLOT_BR <- ggplot2::ggplot(
  Diversity %>% dplyr::filter(
    Measurement == "Richness",
    Environment == "Gamma",
    Space != "Ring"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space),
    group = interaction(ImmigrationRate, ExtirpationRate, Dispersal)
  )
) + ggplot2::geom_line(
) + ggplot2::labs(
  y = "Regional Richness", #Number of Species",
  x = paste0("Time, ", divide_time_by, " units")
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4") #c("darkorange", "plum1", "cyan")
) + ggplot2::theme_bw(
) + ggplot2::coord_cartesian(
  ylim = c(0, 65)
) + ggplot2::facet_grid(
  ImmigrationRate ~ ExtirpationRate
) + ggplot2::ggtitle("Regional Richness")

PLOT_BRR <- ggplot2::ggplot(
  Diversity %>% dplyr::filter(
    Measurement == "Jaccard",
    Environment != "Gamma",
    Space != "Ring",
    Environment != "Mean"
  ),
  mapping = ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space),
    group = interaction(ImmigrationRate, ExtirpationRate,
                        Dispersal, Environment)
  )
) + ggplot2::geom_line(
  alpha = 0.2
) + ggplot2::geom_line(
  data = Diversity %>% dplyr::filter(
    Measurement == "Jaccard",
    Environment != "Gamma",
    Space != "Ring",
    Environment == "Mean"
  ),
  ggplot2::aes(
    group = interaction(ImmigrationRate, ExtirpationRate,
                        Dispersal)
  ),
  size = 1.5
) + ggplot2::geom_ribbon(
  data = DiversityRibbons %>% dplyr::filter(
    Measurement == "Jaccard"
  ),
  mapping = ggplot2::aes(
    ymin = Low,
    ymax = High,
    x = Time,
    fill = pasteCustom(Dispersal, Space)
  ),
  alpha = 0.4,
  inherit.aes = FALSE
) + ggplot2::theme_bw(
) + ggplot2::labs(
  y = "Jaccard Distance", # Number of Species",
  x = paste0("Time, ", divide_time_by, " units")
  # x = ""
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::scale_fill_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4")
#) + ggplot2::theme(
#  legend.position = "none"
) + ggplot2::facet_grid(
  ImmigrationRate ~ ExtirpationRate
) + ggplot2::ggtitle("Spatial Jaccard Turnover")


ggplot2::ggsave(
  filename = "MNA-Image-SupMat-Presence-No.pdf",
  plot = PLOT_TL,
  height = 11, width = 12, dpi = 480, units = "cm"
)
ggplot2::ggsave(
  filename = "MNA-Image-SupMat-Presence-Med.pdf",
  plot = PLOT_TM,
  height = 11, width = 12, dpi = 480, units = "cm"
)
ggplot2::ggsave(
  filename = "MNA-Image-SupMat-Presence-Full.pdf",
  plot = PLOT_TR,
  height = 11, width = 12, dpi = 480, units = "cm"
)
ggplot2::ggsave(
  filename = "MNA-Image-SupMat-Presence-Local.pdf",
  plot = PLOT_BL,
  height = 11, width = 12, dpi = 480, units = "cm"
)
ggplot2::ggsave(
  filename = "MNA-Image-SupMat-Presence-Regional.pdf",
  plot = PLOT_BR,
  height = 11, width = 12, dpi = 480, units = "cm"
)
ggplot2::ggsave(
  filename = "MNA-Image-SupMat-Presence-Jaccard.pdf",
  plot = PLOT_BRR,
  height = 11, width = 12, dpi = 480, units = "cm"
)
