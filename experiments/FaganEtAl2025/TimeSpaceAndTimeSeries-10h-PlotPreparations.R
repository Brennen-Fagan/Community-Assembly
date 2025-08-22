# Settings / Parameters: ######################################################
# Required for interventionStrings and endTimes, which have downstream deps.
if (!exists("ColExt")) {
  load("ColExt10a1_flat.RData")
}

# Problems with X11
options(bitmapType = "cairo")

# Grey interval that we compute over, usually after intervention (~50%)
# If second number is less than 1, we lose persistent species.
end <- c(0.602, 0.9045) # Aiming for 20000 - 30000. These go ~0.0003% away.

defaultNicheDistance <- "5" # "3"::2, "5"::5, "7"::10

richnessYMax <- 42
basePoolPatchSeeds <- as.character(1:44)

# Common Libraries / Functions: ###############################################
library(RMTRCode2) # Personal
library(ggplot2)   # Plotting
library(ggpubr)    # Combining Plots
library(tidytable) # Data Manipulation
# library(tidygraph)
# library(ggraph)
# library(ggforce)

# WISOTT: What it says on the tin.
source("TimeSpaceAndTimeSeries-10-Dictionaries.R") # Defines IDs
source(file.path("R", "changeInterventionLevels.R")) # Land-use names.
source(file.path("R", "changePreferencesLevels.R")) # More human readable WISOTT.
source(file.path("R", "colorPaletteAlg.R")) # Color scheme.
source(file.path("R", "interventionNamingScheme.R")) # WISOTT.
source(file.path("R", "plotGraph.R")) # WISOTT.
source(file.path("R", "plotMeanAndInner.R")) # WISOTT for Value over Time.
source(file.path("R", "unifyAffinityBins.R")) # Make intervals consistent.

# Resources: ##################################################################
interventionMatrix <- matrix(c(
  "(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
  "(0.25)->(0)", "(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
  "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
  "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)", "(0.75)->(1)",
  "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)"
),
byrow = TRUE, nrow = 5)

### colors: ###################################################################
colorPalette <- sapply(interventionMatrix, colorPaletteAlg)

linetypePalette <- c(
  "100% 0" = "solid",
  "50% 0, 50% 1" = "longdash",
  "Uniform(0, 1)" = "dotdash"
)

### renames: ##################################################################
# For presentation (i.e., "Arrival" is a working term, but not 100% accurate.)
externalNames <- c(
  "Arrival"         = "Colonisation",
  "Failed Arrival"  = "Failure",
  "Present"         = "Present",
  "Dispersal"       = "Adjacent",
  "Extinct"         = "Neutral Ext.",
  "Dynamic Loss"    = "Dynamic Ext.",
  "EndOfSimulation" = "Persistent",
  "NA"              = "NA"
)

### Strings: ##################################################################
# Enhance readability, from 9g TablePlots
interventionStrings <- ColExt |> tidytable::select(
  PatchAffinity, PoolPatch, InterventionPatchType
) |> tidytable::distinct(
) |> tidytable::mutate(
  Intervention = unlist(mapply(
    FUN = interventionNamingScheme,
    PatchAffinity, PoolPatch, InterventionPatchType
  ))
)

### End times: #################################################################
# Work out the end times so we can truncate the simulations
# so that we are making sure our comparisons are equivalent.
endTimes <- ColExt |> tidytable::rename(
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
  Start = end[1] * Times, # Neglect anything with an out time before this.
  Stop = end[2] * Times # Neglect anything with an in time after this.
)

stopifnot(all(abs(na.omit(endTimes$ExtraS)) < 0.1)) # Should be ~ 1e-11.
# Example Debugged:
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

endTimes <- endTimes |> tidytable::select(
  -ExtraMn, ExtraS, ExtraN
)

# Verify as we load that the intervention times are calc'd correctedly.
