# Introduction: ###############################################################
# Follows TimeSpaceAndTimeSeries-4b-Sampling.R
# Composed roughly of two sets of plots.
#   1) Plots corresponding to our old figure 3 of diversities.
#   2) Plots corresponding to our samplings, e.g. esp. File 1b.

# Warning: care with combining dat folders.
#          At minimum, the first digit should agree (same simulations).
#          The second and third digits correspond to the intervention types.
#          In particular, 2nd = 1 -> no intervention, so 3rd is meaningless.
datfolders <- c(
  "TSTS_Simulations_1-1-1_2024-01-16",
  "TSTS_Simulations_1-2-1_2024-01-10",
  "TSTS_Simulations_1-2-2_2024-01-10",
  "TSTS_Simulations_1-2-3_2024-01-10",
  "TSTS_Simulations_1-3-1_2024-01-12",
  "TSTS_Simulations_1-3-2_2024-01-12",
  "TSTS_Simulations_1-3-3_2024-01-15"#,
  # "TSTS_Simulations_2-1-6_2024-01-19",
  # "TSTS_Simulations_2-2-6_2024-01-19",
  # "TSTS_Simulations_2-3-6_2024-01-19"#,
  # "TSTS_Simulations_3-1-7_2024-01-22",
  # "TSTS_Simulations_3-2-7_2024-01-22",
  # "TSTS_Simulations_3-3-7_2024-01-22"
)

# Libraries: ##################################################################
library(dplyr)

library(ggplot2)
library(patchwork)

library(RMTRCode2)

# Load Data: ##################################################################

results <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Simulation"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

samples <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Sampling"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

stopifnot(length(results) >= length(datfolders),
          length(samples) == length(datfolders))

# Convert to diversity metrics: ###############################################





