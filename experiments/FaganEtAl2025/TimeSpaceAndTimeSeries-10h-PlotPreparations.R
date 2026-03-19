# Settings / Parameters: ######################################################
# Try to make this more efficient when running multiple images in a single
# session. This function only assigns if the lhs does not exist yet.
`%assign%` <- function(lhs, rhs) {
  lhs <- substitute(lhs)
  env <- parent.frame()
  if (!exists(lhs, envir = env)) {
    assign(deparse(lhs), rhs, envir = env)
  }
}

# Required for interventionStrings and endTimes, which have downstream deps.
# WARNING: This does not respect the dateSwitch below.
# The detected endtimes and intervention strings may therefore be incorrect.
if (!exists("ColExt")) {
  load("ColExt10a1_flat.RData")
}

# Problems with X11
options(bitmapType = "cairo")

dateSwitch <- switch(
  c("Original", "Carrying Capacity", "Timescales")[1],
  "Original" = "2025-07-30",
  "Carrying Capacity" = "2025-09-04",
  "Timescales" = "2026-03-11"
)

# Grey interval that we compute over, usually after intervention (~50%)
# If second number is less than 1, we lose persistent species.
endprops %assign% c(0.602, 0.9045) # Aiming for 20000 - 30000. ~0.0003% away.

defaultNicheDistance %assign% "5" # "3"::2, "5"::5, "7"::10

richnessYMax %assign% 42
basePoolPatchSeeds %assign% as.character(1:44)

dirImages %assign% paste0("TSTS_Images_10j_", dateSwitch)

# Common Libraries / Functions: ###############################################
library(RMTRCode2) # Personal
library(ggplot2)   # Plotting
library(ggpubr)    # Combining Plots
library(tidytable) # Data Manipulation
library(tidygraph) # Graph Manipulations
library(ggraph)    # Graph Plotting

# WISOTT: What it says on the tin.
source("TimeSpaceAndTimeSeries-10-Dictionaries.R") # Defines IDs
source(file.path("R", "changeInterventionLevels.R")) # Land-use names.
source(file.path("R", "changePreferencesLevels.R")) # More human readable WISOTT.
source(file.path("R", "colorPaletteAlg.R")) # Color scheme.
source(file.path("R", "geom_flat_violin.R")) # Half-violin plots.
source(file.path("R", "interventionNamingScheme.R")) # WISOTT.
source(file.path("R", "plotGraph.R")) # WISOTT.
source(file.path("R", "plotMeanAndInner.R")) # WISOTT for Value over Time.
source(file.path("R", "unifyAffinityBins.R")) # Make intervals consistent.

# Resources: ##################################################################
### Folder creation: ##########################################################
if (!dir.exists(dirImages)) {
  dir.create(dirImages, showWarnings = FALSE)
}

### Strings: ##################################################################
# Enhance readability, from 9g TablePlots
interventionStrings %assign% (ColExt |> tidytable::select(
  PatchAffinity, PoolPatch, InterventionPatchType
) |> tidytable::distinct(
) |> tidytable::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    PatchAffinity, PoolPatch, InterventionPatchType
  ))
))

### Standardised intervention names: ##########################################
# WARNING: This should probably be derived from interventionNamingScheme rather
# than be hard-coded (something like outer(initial state, final state)).
interventionMatrix %assign% matrix(c(
  "(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
  "(0.25)->(0)", "(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
  "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
  "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)", "(0.75)->(1)",
  "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)"
),
byrow = TRUE, nrow = 5)

### colors: ###################################################################
colorPalette %assign% sapply(interventionMatrix, colorPaletteAlg)

linetypePalette %assign% c(
  # Internal Aliases
  "100% 0" = "solid",
  "50% 0, 50% 1" = "longdash",
  "Uniform(0, 1)" = "dotdash",
  # Paper Aliases
  "1 Adaptation Type" = "solid",
  "2 Adaptation Types" = "longdash",
  "Multiple Adaptation Types" = "dotdash"
)

### renames: ##################################################################
# For presentation (i.e., "Arrival" is a working term, but not 100% accurate.)
externalNames %assign% c(
  "Arrival"         = "Colonisation",
  "Failed Arrival"  = "Failure",
  "Present"         = "Present",
  "Dispersal"       = "Adjacent",
  "Extinct"         = "Neutral Ext.",
  "Dynamic Loss"    = "Dynamic Ext.",
  "EndOfSimulation" = "Persistent",
  "NA"              = "NA"
)

### End times: #################################################################
# Work out the end times so we can truncate the simulations
# so that we are making sure our comparisons are equivalent.
endTimes %assign% (ColExt |> tidytable::rename(
  DispersalParam = Dispersal
) |> tidytable::filter(
  EventType == "EndOfSimulation"
) |> tidytable::select(
  Times, PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance,
  SpeciesAffinity, SpeciesAffinitySeed, PatchAffinity, PatchAffinitySeed
) |> tidytable::distinct(
) |> tidytable::group_by(
  PoolPatch, PoolPatchSeed, Interactions, InteractionsSeed, Events,
  EventsSeed, InitialConditions, InitialConditionsSeed, DispersalParam,
  NicheDistance,
  SpeciesAffinity, SpeciesAffinitySeed, PatchAffinity, PatchAffinitySeed
) |> tidytable::summarise(
  Times = max(Times),
  ExtraMn = min(Times),
  ExtraS = diff(Times),
  ExtraN = tidytable::n(),
  .groups = "drop"
) |> tidytable::mutate( # In the plots:
  Start = endprops[1] * Times, # Neglect anything with an out time before this.
  Stop = endprops[2] * Times # Neglect anything with an in time after this.
))

if ("ExtraS" %in% colnames(endTimes)) {# Prevent from triggering on repeat
  stopifnot(all(abs(na.omit(endTimes$ExtraS)) < 0.1)) # Should be ~ 1e-11.
  endTimes <- endTimes |> tidytable::select(
    -ExtraMn, -ExtraS, -ExtraN
  )
}
# Example Debugged if there is a problem:
# ColExt |> tidytable::rename(
#   DispersalParam = Dispersal
# ) |> tidytable::filter(
#   EventType == "EndOfSimulation",
#   NicheDistance == "7",
#   SpeciesAffinity == "1",
#   PoolPatchSeed %in% c(26, 29)
# ) |> tidytable::select(
#   -Species, -Environment, -Success, -SpeciesType, -Size,
#   -ReproductionRate, -Affinity, -AffinityBins, -PostIntervention
# ) |> tidytable::distinct(
# ) |> tidytable::arrange(
#   Times
# )


# Verify as we load that the intervention times are calc'd correctedly.
