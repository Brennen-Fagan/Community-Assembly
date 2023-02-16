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
  path = "Data_2022-09-16",
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
names(Properties)[1:6] <- c(
  "MNA", "IterANDRate", "Modifier", "EnvNum", "Space", "DistAND.RData"
)

Properties$FullName <- names(SpeciesPresence)

# Capture the position between the text (first group)
# and the set of numbers (somehow without the +).
# The \\K resets so that we do not capture any text.
patternString <- "((?>[a-zA-Z]+)(?=[0-9eE]))\\K"

# Split strings. Some of the trick will be to introduce
# a character to make the separation around. We use "_".
Properties <- Properties %>% dplyr::mutate(
  IterANDRate = gsub(pattern = patternString,
                     replacement = "_",
                     x = IterANDRate, perl = TRUE),
  Modifier = gsub(pattern = patternString,
                  replacement = "_",
                  x = Modifier, perl = TRUE),
  EnvNum = gsub(pattern = patternString,
                replacement = "_",
                x = EnvNum, perl = TRUE)
) %>% tidyr::separate(
  IterANDRate, into = c("Iter", "Rate"),
  sep = "[_]", fill = "right"
) %>% tidyr::separate(
  Modifier, into = c("Modifier", "ModIntensity"),
  sep = "[_]", fill = "right"
) %>% tidyr::separate(
  EnvNum, into = c("Env", "Environments"),
  sep = "[_]"
) %>% tidyr::separate(
  DistAND.RData, into = c("Distance", ".RData"),
  sep = "[.]"
) %>% dplyr::select(
  -MNA, -.RData, -Env
) %>% dplyr::mutate(
  Distance = dplyr::case_when(
    is.na(Distance) ~ "1e+00",
    TRUE ~ Distance
  )
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
  Modifier, ModIntensity, Space, Distance,
  Species, Time, SizeID
) %>% dplyr::summarise(
  Count = n(), .groups = "drop"
) %>% dplyr::mutate(
  # Distance = ifelse(
  #   Space == "Line" | Space == "Ring", 1, Inf
  # ),
  Distance = 10^as.numeric(Distance),
  Space = Distance,
  Dispersal = 1 - exp( -9 / as.numeric(Space) ),
  Dispersal = paste0(
    formatC(Dispersal))
)

Diversity <- Diversity %>% dplyr::mutate(
  # Distance = ifelse(
  #   Space == "Line" | Space == "Ring", 1, Inf
  # ),
  Distance = 10^as.numeric(Distance),
  Space = Distance,
  Dispersal = 1 - exp( -9 / as.numeric(Space) ),
  Dispersal = paste0(
    formatC(Dispersal))
)

DiversityRibbons <- Diversity %>% dplyr::filter(
  !(Environment %in% c("Mean", "Gamma")),
  Measurement == "Richness" | Measurement == "Jaccard"
) %>%  dplyr::group_by(
  Time, Iter, Distance, Modifier, ModIntensity, Environments, Space, Dispersal,
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
  limits = c(1, 10)
) + ggplot2::facet_grid(
  . ~ #Space +
    Dispersal
) + ggplot2::geom_hline(
  yintercept = 34.5, color = "red"
) + ggplot2::labs(
  y = "Species by Size",
  x = "" # paste0("Time, ", divide_time_by, " units")
) + ggplot2::theme_bw(
) + ggplot2::theme(
  strip.background = ggplot2::element_rect(
    colour = "black",
    fill = "cyan"
  ),
  legend.position = "none",
  axis.text.x = ggplot2::element_blank()
) + ggplot2::scale_y_continuous(
  limits = c(0, 100)
)

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
  limits = c(1, 10)
) + ggplot2::facet_grid(
  . ~ #Space +
    Dispersal
) + ggplot2::geom_hline(
  yintercept = 34.5, color = "red"
) + ggplot2::labs(
  y = "Species by Size",
  x = "" # paste0("Time, ", divide_time_by, " units")
) + ggplot2::theme_bw(
) + ggplot2::theme(
  strip.background = ggplot2::element_rect(
    colour = "black",
    fill = "plum1"
  ),
  legend.position = "none",
  axis.text.x = ggplot2::element_blank()
) + ggplot2::scale_y_continuous(
  limits = c(0, 100)
)

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
) + ggplot2::facet_grid(
  . ~ #Space +
    Dispersal
) + ggplot2::geom_hline(
  yintercept = 34.5, color = "red"
) + ggplot2::labs(
  y = "Species by Size",
  x = "" # paste0("Time, ", divide_time_by, " units")
) + ggplot2::theme_bw(
) + ggplot2::theme(
  strip.background = ggplot2::element_rect(
    colour = "black",
    fill = "darkorange"
  ),
  axis.text.x = ggplot2::element_blank()
) + ggplot2::scale_y_continuous(
  limits = c(0, 100)
)

PLOT_T_Legend <- extractLegend(PLOT_TR)
PLOT_TR <- PLOT_TR + ggplot2::theme(legend.position = "none")

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
    color = pasteCustom(Dispersal, Space)
  )
) + ggplot2::geom_line(
  alpha = 0.4,
  mapping = ggplot2::aes(
    group = interaction(Dispersal, Environment)
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
  # x = paste0("Time, ", divide_time_by, " units")
  x = ""
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::scale_fill_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4")
)


PLOT_B_Legend <- extractLegend(PLOT_BL)
PLOT_BL <- PLOT_BL + ggplot2::theme(legend.position = "none")

PLOT_BR <- ggplot2::ggplot(
  Diversity %>% dplyr::filter(
    Measurement == "Richness",
    Environment == "Gamma",
    Space != "Ring"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space)
  )
) + ggplot2::geom_line(
) + ggplot2::labs(
  y = "Regional Richness", #Number of Species",
  x = paste0("Time, ", divide_time_by, " units")
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::theme_bw(
) + ggplot2::coord_cartesian(
  ylim = c(0, 45)
) + ggplot2::theme(
  legend.position = "none"
)

PLOT_BRR <- ggplot2::ggplot(
  Diversity %>% dplyr::filter(
    Measurement == "Jaccard",
    Environment != "Gamma",
    Space != "Ring",
    Environment != "Mean"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space)
  )
) + ggplot2::geom_line(
  alpha = 0.2,
  mapping = ggplot2::aes(
    group = interaction(Dispersal, Environment)
  )
) + ggplot2::geom_line(
  data = Diversity %>% dplyr::filter(
    Measurement == "Jaccard",
    Environment != "Gamma",
    Space != "Ring",
    Environment == "Mean"
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
  y = "Jaccard Turnover", # Number of Species",
  # x = paste0("Time, ", divide_time_by, " units")
  x = ""
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::scale_fill_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4")
) + ggplot2::theme(legend.position = "none")

blankGrob <- ggplot2::ggplot() + ggplot2::theme_void()

PLOT_LegendGrob <- gridExtra::arrangeGrob(
  PLOT_T_Legend, PLOT_B_Legend,
  blankGrob, blankGrob,
  ncol = 2, nrow = 2,
  widths = c(1, 1),
  heights = c(1, 0.2)
)

obj <- gridExtra::arrangeGrob(
  PLOT_TL, PLOT_TM, PLOT_TR, PLOT_T_Legend,
  PLOT_BL, PLOT_BR, PLOT_BRR, PLOT_B_Legend,
  ncol = 4, nrow = 2, widths = c(1, 1, 1, 0.5)
)
# objobj <- gridExtra::arrangeGrob(
#   obj, PLOT_LegendGrob,
#   ncol = 2, widths = c(10, 2)
# )

# Alternate Bottom
Diversity <- Diversity %>% dplyr::mutate(
  Measurement2 = dplyr::case_when(
    Measurement == "Jaccard" ~ "Turnover",
    Measurement == "Richness" & Environment == "Gamma" ~ "Regional Rich.",
    Measurement == "Richness" & Environment == "Mean" ~ "Local Rich.", # Panel
    Measurement == "Richness"  ~ "Local Rich.", # Otherwise
    TRUE ~ Measurement
  )
)

PLOT_TALT_pre <- ggplot2::ggplot(
  SpeciesPresence %>% dplyr::mutate(
    Dispersal = dplyr::case_when(
      Space == 1 ~ "Full Dispersal",
      Space == 1e+05 ~ "Med. Dispersal",
      Space == Inf ~ "No Dispersal"
    )
  ),
  ggplot2::aes(x = Time, y = SizeID, color = Count)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1,
  limits = c(1, 10)
) + ggplot2::facet_grid(
  . ~ #Space +
    factor(Dispersal, ordered = T,
           levels = c("No Dispersal", "Med. Dispersal", "Full Dispersal")),
  scales = "free_y"
) + ggplot2::geom_hline(
  yintercept = 34.5, color = "red"
) + ggplot2::labs(
  y = "Species by Size",
  x = paste0("Time, ", divide_time_by, " units"),
  tag = "a)"
) + ggplot2::theme_bw(
) + ggplot2::theme(
  axis.text.x = ggplot2::element_blank(),
  plot.tag.position = c(0.02, 0.98)
) + ggplot2::scale_y_continuous(
  limits = c(0, 100)
)

# https://stackoverflow.com/a/65024951
PLOT_TALT_post <-
  ggplot2::ggplot_gtable(ggplot2::ggplot_build(PLOT_TALT_pre))
strips <- which(startsWith(PLOT_TALT_post$layout$name, 'strip'))
for (s in seq_along(strips)) {
  PLOT_TALT_post$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- c(
    "cyan", "plum1", "darkorange"
  )[s]
}

PLOT_BALT <- ggplot2::ggplot(
  Diversity %>% dplyr::filter(
    Measurement2 %in% c("Turnover", "Local Rich.", "Regional Rich."),
    Environment != "Mean"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space)
  )
) + ggplot2::geom_line(
  alpha = 0.4,
  mapping = ggplot2::aes(
    group = interaction(Dispersal, Environment)
  )
) + ggplot2::geom_line(
  data = Diversity %>% dplyr::filter(
    Measurement2 %in% c("Turnover", "Local Rich.", "Regional Rich."),
    Environment == "Mean"
  ),
  size = 1.5
) + ggplot2::geom_ribbon(
  data = DiversityRibbons %>% dplyr::mutate(
    Measurement2 = dplyr::case_when(
      Measurement == "Jaccard" ~ "Turnover",
      Measurement == "Richness" ~ "Local Rich.",
      TRUE ~ Measurement
    )
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
  y = "Value", # Number of Species",
  x = paste0("Time, ", divide_time_by, " units"),
  tag = "b)"
  # x = ""
) + ggplot2::theme(
  plot.tag.position = c(0.02, 0.98)
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::scale_fill_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4")
) + ggplot2::facet_wrap(
  . ~ factor(
    Measurement2, ordered = T,
    levels = c("Local Rich.", "Regional Rich.", "Turnover")
  ), nrow = 1, scales = "free_y"
)

objalt <- gridExtra::arrangeGrob(
  # gridExtra::arrangeGrob(
  #   PLOT_TL, PLOT_TM, PLOT_TR, PLOT_T_Legend,
  #   ncol = 4, nrow = 1, widths = c(1, 1, 1, 0.5)
  # ),
  PLOT_TALT_post,
  PLOT_BALT, nrow = 2
)

ggplot2::ggsave(
  filename = "MNA-Image-Example-Presence.png",
  plot = objalt,
  height = 11, width = 12, dpi = 480, units = "cm"
)

# # Extract Sources of Abundance: ###############################################
# # Two primary sources: the Events list and the Abundance matrix.
# # While I expect to throw away all species information at the end of the
# # calculation, the calculation itself makes more sense while it is preserved.
#
#
#
# # For a given time t_n, we calculate the change in abundance due to a
# # particular flux element f as
# #    (f(t_n) + f(t_n + 1)) / 2 * (t_(n + 1) - t_(n)).
# # Note that, since we have all of the results already, we know what each
# # abundance is at each t_n and need only worry about accurate enough f's.
#
# # Note implicitly that I have taken the premise of assigning the result to the
# # time of its resolution, rather than the time of its generation (n+1, not n).
# centralDiff <- function(x, timediff) {
#   # Applied to each column separately.
#   c(NA, timediff * (x[-1] + dplyr::lag(x)[-1]) / 2) * divide_time_by
# }
#
# AbundanceChanges <- lapply(
#   seq_along(Results),
#   function(i, results, pool, mats, distmats) {
#     result <- results[[i]]
#
#
#     # Flux/Continuous Sources/Sinks: ###########################################
#     TimeDiffs <- diff(result$Abundance[, 1])
#
#     # Abundance * Per Capita Type 1 + DispersalMatrix %*% Abundance
#     # PCT1 = rep(ReproductionRate, NumEnvironments) + InteractionMatrix %*% y
#     # DispersalMatrix = CreateDispersalMatrix(EnvironmentDistances, Speeds)
#
#     InteractionMatrix <- Matrix::bdiag(mats$Mats)
#     InteractionMatrixPos <- InteractionMatrix
#     InteractionMatrixPos@x <- ifelse(InteractionMatrix@x > 0, InteractionMatrix@x, 0)
#     InteractionMatrixNeg <- InteractionMatrix
#     InteractionMatrixNeg@x <- ifelse(InteractionMatrix@x < 0, InteractionMatrix@x, 0)
#
#     DispersalMatrix <- RMTRCode2::CreateDispersalMatrix(distmats[[i]],
#                                                         rep(1, nrow(pool)))
#     DispersalMatrixPos <- DispersalMatrix
#     DispersalMatrixPos@x <- ifelse(DispersalMatrix@x > 0, DispersalMatrix@x, 0)
#     DispersalMatrixNeg <- DispersalMatrix
#     DispersalMatrixNeg@x <- ifelse(DispersalMatrix@x < 0, DispersalMatrix@x, 0)
#
#     # Need to t(ranspose) the results, since it applies to each row and then
#     # populates a matrix by columns.
#     BirthDeath <- t(apply(
#       result$Abundance[, -1], 1,
#       function(abundance) abundance * rep(pool$ReproductionRate, 10)
#     ))
#     Birth <- ifelse(BirthDeath > 0, BirthDeath, 0)
#     Death <- ifelse(BirthDeath < 0, BirthDeath, 0)
#
#     InteractionsPos <- do.call(rbind, apply(
#       result$Abundance[, -1], 1,
#       function(abundance) Matrix::t((InteractionMatrixPos %*% abundance) * abundance)
#     ))
#
#     InteractionsNeg <- do.call(rbind, apply(
#       result$Abundance[, -1], 1,
#       function(abundance) Matrix::t((InteractionMatrixNeg %*% abundance) * abundance)
#     ))
#
#     DispersalPos <- do.call(rbind, apply(
#       result$Abundance[, -1], 1,
#       function(abundance) Matrix::t(DispersalMatrixPos %*% abundance)
#     ))
#
#     DispersalNeg <- do.call(rbind, apply(
#       result$Abundance[, -1], 1,
#       function(abundance) Matrix::t(DispersalMatrixNeg %*% abundance)
#     ))
#
#     FluxChanges <- lapply(
#       list(
#         Birth = Birth, Death = Death,
#         InteractionsPos = InteractionsPos, InteractionsNeg = InteractionsNeg,
#         DispersalPos = DispersalPos, DispersalNeg = DispersalNeg
#       ),
#       function(m) {
#         m %>% as.matrix %>% data.frame %>% dplyr::mutate(
#           dplyr::across(.fns = centralDiff, timediff = TimeDiffs)
#         )
#       }
#     )
#
#     # Events/Discrete Sources/Sinks: ##########################################
#     # Note that, since these occur at discrete times, they have precise timing,
#     # compared to the continuous case (which we right assigned instead of left).
#
#     # We need to essentially convert Events into a sparse matrix:
#     # Times -> Row, Species + Environment -> Column, Type + Success -> Matrix
#     # Entries are then drawn from result$Abundance for Extinctions, 0.4 o.w.
#
#     EventChanges <- list(
#       Immigration = NA, Establishment = NA, Extinction = NA
#     )
#     EventChanges <- lapply(EventChanges, function(m, dim) {
#       return(matrix(nrow = dim[1], ncol = dim[2]))
#     },
#     dim = dim(Birth)
#     )
#
#
#     # Rescale time so that they are both operating on the same units.
#     # Discard the burn-in by removing anything that happens before the first
#     # recorded abundance step.
#     result$Events <- result$Events %>% dplyr::mutate(
#       Times = Times / divide_time_by
#     ) %>% dplyr::filter(
#       Times >= result$Abundance[1, 1]
#     )
#
#     for (r in 1:nrow(result$Events)) {
#       event <- result$Events[r, ]
#
#       # Why which.max, not which? During filtering I think some values were
#       # combined, resulting in misplacement of some exact times.
#       time <- which.max(result$Abundance[, 1] >= event$Times)
#       id <- (event$Environment - 1) * 100 + event$Species
#
#       if (event$Type == "Arrival") {
#         EventChanges$Immigration[time, id] <- result$Parameters$ArrivalDensity
#         if (event$Success == TRUE) {
#           EventChanges$Establishment[time, id] <- result$Parameters$ArrivalDensity
#         }
#       } else if (event$Type == "Extinct") {
#         if (event$Success == TRUE && time - 1 != 0) { # + 1 for Time Column.
#           EventChanges$Extinction[time, id] <- -result$Abundance[time - 1, id + 1]
#         }
#       }
#     }
#
#     return(c(
#       FluxChanges, EventChanges
#     ))
#   },
#   results = Results,
#   pool = PoolsMats$`MNA-FirstAttempt-PoolMats-Env10.RData`$Pool,
#   mats = PoolsMats$`MNA-FirstAttempt-PoolMats-Env10.RData`$InteractionMatrices,
#   distmats = list(
#     Matrix::bandSparse(
#       10, k = c(-1, 1),
#       diagonals = list(rep(1, 10 - 1),
#                        rep(1, 10 - 1))
#     ),
#     Matrix::bandSparse(
#       10, k = c(-1, 1),
#       diagonals = list(rep(Inf, 10 - 1),
#                        rep(Inf, 10 - 1))
#     )
#   )
# )
#
# tableChanges <- lapply(AbundanceChanges, function(type) {
#   lapply(type, sum, na.rm = TRUE)
# })
#
# tableChangesPercent <- lapply(tableChanges, function(type) {
#   pos <- 0; neg <- 0
#   for (i in type) {
#     if (i > 0) pos <- pos + i
#     if (i < 0) neg <- neg + i
#   }
#   for (i in seq_along(type)) {
#     if (type[[i]] > 0) type[[i]] <- type[[i]] / pos
#     if (type[[i]] < 0) type[[i]] <- type[[i]] / neg
#   }
#   return(type)
# })
#
# # Test Run for Occupation Times: ###############################################
# # Idea is to look at the Abundance matrix (make sure it is corrected below
# # eliminations!) and look for positive runs, which are then assigned the time
# # difference between the start and end of the run. These run lengths should
# # form a distribution, which we can then present.
#
# occupationTimes <- lapply(
#   Results, function(result) {
#     # Calculate Runs
#     runs <- apply(result$Abundance[, -1], 2, function(x) {
#       positives <- which(x > 0)
#       if (length(positives) == 0) {return(NULL)}
#
#       bigsteps <- which(diff(positives) > 1)
#
#       ends <- c(positives[bigsteps], positives[length(positives)])
#
#       starts <- c(positives[1], positives[bigsteps + 1])
#
#       runs <- data.frame(
#         start = starts,
#         end = ends
#       ) %>% dplyr::filter(
#         start != end
#       )
#
#       return(runs)
#     },
#     simplify = FALSE)
#
#     # Add Timing Information
#     runs <- lapply(runs, function(df, times) {
#       if (is.null(df)) {return(NULL)}
#       df %>% dplyr::mutate(
#         startTime = times[start],
#         endTime = times[end],
#         diffs = endTime - startTime
#       )
#     },
#     times = result$Abundance[, 1])
#
#     # Add Identifying Information
#     # Species, Environment, Basal/Consumer, LastAbundance (=> Neutral/Dynamic)
#     runs <- lapply(
#       seq_along(runs),
#       function(i, dfs, species, environments, abundance) {
#         if (is.null(dfs[[i]])) {return(NULL)}
#         dfs[[i]] %>% dplyr::mutate(
#           Species = ((i - 1) %% species) + 1,
#           Environment = ((i - 1) %/% species) + 1,
#           Type = dplyr::case_when(
#             Species <= 34 ~ "Basal",
#             Species > 34 ~ "Consumer",
#             TRUE ~ "Oops"
#           ),
#           LastAbundance = abundance[end, i]
#         )
#       },
#       dfs = runs, species = 100, environments = 10,
#       abundance = result$Abundance[, -1])
#
#     return(dplyr::bind_rows(runs))
#   }
# )
#
# occupationTimes$`MNA-FirstAttempt-Result-Env10-Line.RData` <-
#   occupationTimes$`MNA-FirstAttempt-Result-Env10-Line.RData` %>% dplyr::mutate(
#     Space = 1
#   )
#
# occupationTimes$`MNA-FirstAttempt-Result-Env10-None.RData` <-
#   occupationTimes$`MNA-FirstAttempt-Result-Env10-None.RData` %>% dplyr::mutate(
#     Space = Inf
#   )
#
# occupationTimes <- dplyr::bind_rows(occupationTimes)
#
# ggplot2::ggsave(
#   plot = ggplot2::ggplot(
#     occupationTimes,
#     ggplot2::aes(
#       x = diffs,
#       y = LastAbundance
#     )
#   ) + ggplot2::geom_bin2d(
#   ) + ggplot2::facet_wrap(
#     Space ~ ., nrow = 2
#   ) + ggplot2::scale_fill_viridis_c(
#     direction = -1, trans = "log10"
#   ) + ggplot2::scale_x_log10(
#   ) + ggplot2::scale_y_log10(
#   ) + ggplot2::theme_bw(
#   ) + ggplot2::labs(
#     x = paste0("Length of Run, ", divide_time_by, " units"),
#     y = "Abundance at End of Run",
#     title = "Example Occupation Times and Final Abundances"
#   ),
#   filename = "MNA-Image-Example-Occupations.png",
#   height = 6, width = 8, dpi = 320, units = "in"
# )
#
# occupationTimesLong <- tidyr::pivot_longer(
#   occupationTimes, cols = c("startTime", "endTime"),
#   names_to = "EventType", values_to = "Time"
# ) %>% dplyr::group_by(
#   Space
# ) %>% dplyr::arrange(
#   Time
# ) %>% dplyr::mutate(
#   Effect = dplyr::case_when(
#     EventType == "startTime" ~ 1,
#     EventType == "endTime" ~ -1,
#     TRUE ~ as.double(NA_integer_)
#   ),
#   Occupancy = cumsum(Effect),
#   AverageTurnoverToDateTotal = cumsum(abs(Effect)) / (Time * divide_time_by)
# ) %>% dplyr::group_by(
#   Space, Effect
# ) %>% dplyr::mutate(
#   AverageTurnoverToDate = cumsum(Effect) / (Time * divide_time_by)
# )
#
# ggplot2::ggsave(
#   plot = ggplot2::ggplot(
#     occupationTimes
#   ) + ggplot2::geom_histogram(
#     # EXITS
#     ggplot2::aes(x = endTime, y = -1*..count..),
#     fill = "red", binwidth = 0.05
#   ) + ggplot2::geom_histogram(
#     # ENTRANCES
#     ggplot2::aes(x = startTime),
#     fill = "blue", binwidth = 0.05
#   ) + ggplot2::geom_histogram(
#     # ABSOLUTE VALUE
#     data = occupationTimesLong,
#     mapping = ggplot2::aes(x = Time),
#     fill = "transparent", color = "black", binwidth = 0.05
#   ) + ggplot2::geom_line(
#     # DIFFERENCE
#     data = occupationTimesLong,
#     mapping = ggplot2::aes(x = Time, y = Occupancy, color = as.factor(Space))
#   ) + ggplot2::facet_wrap(
#     . ~ Space, nrow = 2
#   ) + ggplot2::labs(
#     y = "Species Change",
#     x = paste0("Time, ", divide_time_by, " units"),
#     caption = paste(
#       "Occupancy refers to the number of species-island pairs with positive abundance.",
#       "Bins are 0.05 wide.", sep = "\n"
#     )
#   ) + ggplot2::scale_y_continuous(
#     minor_breaks = seq(from = -100, to = 120, by = 10)
#   ) + ggplot2::theme_bw(
#   ) + ggplot2::scale_color_manual(
#     name = "Color Guide",
#     breaks = c("red", "blue", "transparent", "1", "Inf"),
#     limits = c("red", "blue", "transparent", "1", "Inf"),
#     values = c("red", "blue", "transparent", "darkorange", "cyan"),
#     labels = c("Loss", "Gain", "Total", "Occupancy", "Occupancy"),
#     drop = FALSE
#   ) + ggplot2::theme(
#     legend.key = ggplot2::element_rect(color = "black")
#   ),
#   filename = "MNA-Image-Example-Turnover.png",
#   height = 6, width = 8, dpi = 320, units = "in"
# )
#
# ggplot2::ggsave(
#   plot = ggplot2::ggplot(
#     occupationTimesLong,
#     ggplot2::aes(x = Time, y = AverageTurnoverToDate, color = factor(Effect))
#   ) + ggplot2::geom_line(
#   ) + ggplot2::geom_line(
#     mapping = ggplot2::aes(x = Time, y = AverageTurnoverToDateTotal),
#     color = "black"
#   ) + ggplot2::facet_wrap(
#     . ~ Space, nrow = 2
#   ) + ggplot2::theme_bw(
#   ) + ggplot2::scale_color_manual(
#     name = "Color Guide",
#     breaks = c("-1", "1", "black"),
#     limits = c("-1", "1", "black"),
#     values = c("red", "blue", "black"),
#     labels = c("Loss", "Gain", "Total"),
#     drop = FALSE
#   ) + ggplot2::labs(
#     y = "Average Turnover to Date (Species / 1 Time Unit)",
#     x = paste0("Time, ", divide_time_by, " units"),
#   ),
#   filename = "MNA-Image-Example-AverageTurnover.png",
#   height = 6, width = 8, dpi = 320, units = "in"
# )
