# Almost a copy of Mutualism-ExampleOutcome-Image up until the plot.
# Difference 0 is removal of functions that are now in package.
# Difference 1 is in the Properties (and number of files targeted).
# We preserve the initial Column which dictates the Simulation Run.
# Difference 2 is that we now have multiple sim types, so we need to adjust the
# Species Typing to reflect this.

library(RMTRCode2)
library(dplyr)        # Data Manipulation
library(tidyr)        # Data Pivotting
library(ggplot2)      # 2-D Plot
library(gridExtra)
library(patchwork)

# Problems with X11
options(bitmapType = "cairo")

by_for_thinning <- 10 # time steps
divide_time_by <- 1E1 # time units
burn_in <-0 # 1E1 # time units

# Implementation: ##############################################################

# All .RData
files_dat <- dir(
  path = c("Data_2023-06-26", "Data_2023-06-27", "Data_2023-07-06",
           "Data2_2023-06-30"), # "Data_2022-09-16",
  pattern = ".+Example.+[.]RData$", # "MNA[-]ExampleOutcome[-].+[.]RData$",
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
  RMTRCode2::load_safe_thin,
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

# In principle, we should change this to a lapply:
#   lapply(PoolsMats, function(x) table(x$Pool$Type))
# but for now, note that 68 on first layer, 132 on second.
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
Properties <- data.frame(
  do.call(rbind, Properties),
  stringsAsFactors = FALSE
)
names(Properties)[1:9] <- c(
  "SimType", "ExampleExtProp", "Result", "EnvNum",
  "Space", "Distance", "Imm", "Ext", "ExtPropAND.RData"
)

Properties$FullName <- names(SpeciesPresence)

# Capture the position between the text (first group)
# and the set of numbers (somehow without the +).
# The \\K resets so that we do not capture any text.
patternString <- "((?>[a-zA-Z]+)(?=[0-9eE]))\\K"

# Split strings. Some of the trick will be to introduce
# a character to make the separation around. We use "_".
Properties <- Properties %>% dplyr::mutate(
  ExtPropAND.RData = gsub(pattern = patternString,
                          replacement = "_",
                          x = ExtPropAND.RData, perl = TRUE),
  EnvNum = gsub(pattern = patternString,
                replacement = "_",
                x = EnvNum, perl = TRUE)
) %>% tidyr::separate(
  ExtPropAND.RData, into = c("ExtProp", "ExtinctionProportionAND.RData"),
  sep = "[_]", fill = "right"
) %>% tidyr::separate(
  EnvNum, into = c("Env", "Environments"),
  sep = "[_]"
) %>% tidyr::separate(
  ExtinctionProportionAND.RData, into = c("ExtinctionProportion", ".RData"),
  sep = "[.]"
) %>% dplyr::select(
  -ExampleExtProp, -Result, -.RData, -Env, -ExtProp
) %>% dplyr::mutate(
  Distance = gsub(Distance, pattern = "_", replacement = "-", fixed = TRUE),
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

SpeciesTypes <- lapply(
  seq_along(PoolsMats), function(i, x, n) {
    x[[i]]$Pool %>% dplyr::arrange(
      ReproductionRate
      ) %>% dplyr::mutate(
        SimType = n[i],
      RepRateID = 1:nrow(x[[i]]$Pool)
      )
  }, x = PoolsMats,
  n = sub(pattern = "[-].+", replacement = "", x = basename(names(PoolsMats)))
) %>% dplyr::bind_rows()

SpeciesPresence <- dplyr::left_join(
  SpeciesPresence,
  SpeciesTypes,
  by = c("SimType", "Species" = "ID")
)


SpeciesPresence <- SpeciesPresence %>% dplyr::group_by(
  SimType, Imm, Ext, ExtinctionProportion, Space, Distance,
  Species, Time, RepRateID
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
  SimType, Time, Distance, Imm, Ext, ExtinctionProportion, Environments,
  Space, Dispersal,
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
  SimType, Time, Distance, Imm, Ext, ExtinctionProportion, Environments,
  Space, Dispersal,
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
  factor(x, levels = c(
    "0", "2e-09", "2e-08", "2e-07", "2e-06", "2e-05",
    "0.0002", "0.001998", "0.0198", "0.1813", "0.8647"
  ), ordered = TRUE)
  # paste0("(", x, ", ", y, ")")
  # ifelse(y == max(y), "None",
  #        ifelse(y == min(y), "Full", "Med."))
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

PLOT_Presence <- lapply(
  SpeciesPresence$SimType %>% unique(),
  function(SType) {
    ggplot2::ggplot(
      SpeciesPresence %>% dplyr::filter(SimType == SType, Dispersal != 1),
      ggplot2::aes(x = Time, y = RepRateID, color = Count)
    ) + ggplot2::geom_point(
      shape = '.'
    ) + ggplot2::scale_color_viridis_c(
      direction = -1,
      limits = c(1, 10)
    ) + ggplot2::facet_grid(
      SimType ~ #Space +
        factor(Dispersal, levels = c(
          "0", "2e-09", "2e-08", "2e-07", "2e-06", "2e-05",
          "0.0002", "0.001998", "0.0198", "0.1813", "0.8647"
        ), ordered = TRUE),
      scales = "free_y"
    ) + ggplot2::geom_hline(
      yintercept = 1/3 * nrow(PoolsMats[[1]]$Pool), color = "red"
    ) + ggplot2::labs(
      y = "Species by Size",
      x = paste0("Time, ", divide_time_by, " units"),
      tag = "a)"
    ) + ggplot2::theme_bw(
    ) + ggplot2::theme(
      #axis.text.x = ggplot2::element_blank(),
      plot.tag.position = c(0.02, 0.98)
    )
  }
)

# ggplot2::ggsave(
#   patchwork::wrap_plots(PLOT_Presence, nrow = length(PLOT_Presence), ncol = 1),
#   "PresentPlots.png", units = "in", width = 12, height = 9
# )

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
# ) + ggplot2::scale_color_manual(
#   name = legend_bl_name,
#   values = c("darkorange", "plum1", "cyan")
# ) + ggplot2::scale_fill_manual(
#   name = legend_bl_name,
#   values = c("darkorange4", "plum4", "cyan4")
) + ggplot2::facet_wrap(
  SimType ~ factor(
    Measurement2, ordered = T,
    levels = c("Local Rich.", "Regional Rich.", "Spatial Diss.")
  ), #nrow = unique(Diversity$SimType),
  ncol = 3,
  scales = "free"
) + ggplot2::scale_alpha(guide = "none") + ggplot2::coord_cartesian(
  ylim = c(0, NA)
)


# ggplot2::ggsave(
#     PLOT_B,
#   "DiversityPlots.png", units = "in", width = 12, height = 9
# )
