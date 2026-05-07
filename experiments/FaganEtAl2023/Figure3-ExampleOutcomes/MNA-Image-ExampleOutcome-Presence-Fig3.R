# Extracted from "FirstAttempt-Doc-Analysis2-Gallery3.Rmd".
# Then modified to comply with the standards in
# "Viking_HandleDiversity_ParametersAndPlots3.R".
# This should result in a 2x2 figure.
# Row 1: Species Presence Absence, by Size.
# Row 2: Diversity Metrics over Time.

library(RMTRCode2)    # Helper Functions
library(dplyr)        # Data Manipulation
library(tidyr)        # Data Pivotting
library(ggplot2)      # 2-D Plot
library(gridExtra)



# Problems with X11
options(bitmapType = "cairo")

by_for_thinning <- 100 # time steps
divide_time_by <- 1E4 # time units
burn_in <- 1E4 # time units

pathRfunctions <- file.path("..", "R")
source(file.path(pathRfunctions, "Calculate_Diversity.R"))
source(file.path(pathRfunctions, "Calculate_Species.R"))
source(file.path(pathRfunctions, "load_safe.R"))
source(file.path(pathRfunctions, "load_safe_thin.R"))

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

Properties <- gsub(pattern = ".RData", replacement = "",
                   x = names(SpeciesPresence), fixed = TRUE)
Properties <- strsplit(Properties, '-',
                       fixed = TRUE)
Properties <- lapply(Properties,
                     function(x) if (length(x) < 10) {c(x, NA)} else {x})
# Note the mix of Keyword and Location Structure (oops).
# (E.g. Dates DD-MM-YYYY, Decimals 1.35e-05.)
Properties <- data.frame(
  do.call(rbind, Properties),
  stringsAsFactors = FALSE
)
names(Properties)[1:10] <- c(
  "PathMNA", "ExampleExtProp", "Result", "EnvNum", "Space",
  "Distance", "Imm", "Ext", "ExtProp", "MaybeSubset"
)

Properties$FullName <- names(SpeciesPresence)

# Capture the position between the text (first group)
# and the set of numbers (somehow without the +).
# The \\K resets so that we do not capture any text.
patternString <- "((?>[a-zA-Z]+)(?=[0-9eE]))\\K"

# Split strings. Some of the trick will be to introduce
# a character to make the separation around. We use "_".
Properties <- Properties %>% dplyr::mutate(
  ExtProp = gsub(pattern = patternString,
                     replacement = "_",
                     x = ExtProp, perl = TRUE),
  EnvNum = gsub(pattern = patternString,
                replacement = "_",
                x = EnvNum, perl = TRUE)
) %>% tidyr::separate(
  ExtProp, into = c("ExtProp", "ExtProportion"),
  sep = "[_]", fill = "right"
) %>% tidyr::separate(
  EnvNum, into = c("Env", "Environments"),
  sep = "[_]"
) %>% dplyr::select(
  -PathMNA, -ExampleExtProp, -Result, -Env, -ExtProp, -MaybeSubset
) %>% dplyr::mutate(
  Distance = dplyr::case_when(
    is.na(Distance) ~ "1e+00",
    TRUE ~ Distance
  )
) %>% dplyr::distinct(
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
  # Modifier, ModIntensity,
  Imm, Ext, ExtProportion,
  Space, Distance,
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
  Time, #Iter,
  Distance,
  Imm, Ext, ExtProportion,
  Environments, Space, Dispersal,
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
  Time, #Iter,
  Distance,
  Imm, Ext, ExtProportion,
  Environments, Space, Dispersal,
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
  tag = "(a)"
) + ggplot2::theme_bw(
) + ggplot2::theme(
  axis.text.x = ggplot2::element_blank(),
  plot.tag.position = c(0.02, 0.98),
  plot.tag = ggplot2::element_text(face = "bold")
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
  tag = "(b)"
  # x = ""
) + ggplot2::theme(
  plot.tag.position = c(0.02, 0.98),
  plot.tag = ggplot2::element_text(face = "bold"),
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


