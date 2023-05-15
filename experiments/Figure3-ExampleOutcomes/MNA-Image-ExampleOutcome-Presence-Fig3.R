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
  path = "Data_Figure3", # "Data_2023-03-01", # "Data_2022-09-16",
  pattern = "MNA-Example.+[.]RData$", # "MNA[-]ExampleOutcome[-].+[.]RData$",
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
  Dispersal = 1 - exp( -2 / as.numeric(Space) ),
  Dispersal = paste0(
    formatC(Dispersal))
)

Diversity <- Diversity %>% dplyr::mutate(
  # Distance = ifelse(
  #   Space == "Line" | Space == "Ring", 1, Inf
  # ),
  Distance = 10^as.numeric(Distance),
  Space = Distance,
  Dispersal = 1 - exp( -2 / as.numeric(Space) ),
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

DiversityRibbons_Gamma <- Diversity %>% dplyr::filter(
  (Environment %in% c("Gamma")),
  Measurement == "Richness"
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
) %>% dplyr::mutate(
  Measurement = "Regional Rich."
)

# Plots: #######################################################################

pasteCustom <- function(x, y) {
  #paste0("(", x, ", ", y, ")")
  ifelse(is.infinite(y), "None",
         ifelse(y == 1, "Full", "Med."))
}
legend_bl_name <- "Dispersal"

Diversity <- Diversity %>% dplyr::mutate(
  Measurement2 = dplyr::case_when(
    Measurement == "Jaccard" ~ "Spatial Diss.",
    Measurement == "Richness" & Environment == "Gamma" ~ "Regional Rich.",
    Measurement == "Richness" & Environment == "Mean" ~ "Local Rich.", # Panel
    Measurement == "Richness"  ~ "Local Rich.", # Otherwise
    TRUE ~ Measurement
  )
)

PLOT_T_pre <- ggplot2::ggplot(
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
PLOT_T_post <-
  ggplot2::ggplot_gtable(ggplot2::ggplot_build(PLOT_T_pre))
strips <- which(startsWith(PLOT_T_post$layout$name, 'strip'))
for (s in seq_along(strips)) {
  PLOT_T_post$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- c(
    "cyan", "plum1", "darkorange"
  )[s]
}

PLOT_B <- ggplot2::ggplot(
  Diversity %>% dplyr::filter(
    Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    Environment != "Mean"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = pasteCustom(Dispersal, Space)
  )
) + ggplot2::geom_line(
  # alpha = 0.4,
  mapping = ggplot2::aes(
    group = interaction(Dispersal, Environment),
    alpha = ifelse(Measurement2 == "Regional Rich.", 1, 0.4)
  )
) + ggplot2::geom_line(
  data = Diversity %>% dplyr::filter(
    Measurement2 %in% c("Spatial Diss.", "Local Rich.", "Regional Rich."),
    Environment == "Mean"
  ),
  size = 1.5
) + ggplot2::geom_ribbon(
  data = dplyr::bind_rows(
    DiversityRibbons,
    DiversityRibbons_Gamma
  ) %>% dplyr::mutate(
    Measurement2 = dplyr::case_when(
      Measurement == "Jaccard" ~ "Spatial Diss.",
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
  plot.tag.position = c(0.02, 0.98),
  strip.text.x = ggplot2::element_text(size = 8)
) + ggplot2::scale_color_manual(
  name = legend_bl_name,
  values = c("darkorange", "plum1", "cyan")
) + ggplot2::scale_fill_manual(
  name = legend_bl_name,
  values = c("darkorange4", "plum4", "cyan4")
) + ggplot2::facet_wrap(
  . ~ factor(
    Measurement2, ordered = T,
    levels = c("Local Rich.", "Regional Rich.", "Spatial Diss.")
  ), nrow = 1, scales = "free_y"
) + ggplot2::scale_alpha(guide = "none") + ggplot2::coord_cartesian(
  ylim = c(0, NA)
)

obj <- gridExtra::arrangeGrob(
  PLOT_T_post,
  PLOT_B, nrow = 2
)

ggplot2::ggsave(
  filename = "MNA-Image-Example-Presence.pdf",
  plot = obj,
  height = 11, width = 12, dpi = 480, units = "cm"
)


