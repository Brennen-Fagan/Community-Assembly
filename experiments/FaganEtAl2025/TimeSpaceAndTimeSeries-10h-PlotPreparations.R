# Settings / Parameters: ######################################################
load("diversitiesFlattened10a1_subsetRichness.RData")
load("ColExt10a1_flat.RData")
load("TSTS_Interventions_10a1.RData")

# Problems with X11
options(bitmapType = "cairo")

# Grey interval that we compute over, usually after intervention (~50%)
# If second number is less than 1, we lose persistent species.
end <- c(0.602, 0.9045) # Aiming for 20000 - 30000. These go ~0.0003% away.

defaultNicheDistance <- "5" # "3"::2, "5"::5, "7"::10

# Libraries: ##################################################################
library(RMTRCode2)
library(ggplot2)
library(ggpubr)
library(tidytable)
library(tidygraph)
library(ggraph)
library(ggforce)

source("TimeSpaceAndTimeSeries-10-Dictionaries.R")
# source('TimeSpaceAndTimeSeries-0-Functions.R')
# source("flattenDiversity.R") # Flatten the data for the networks.
# source("CalculateTrophicStructure.R") # Calculator creator.
# source("toCheddar.R") # Updated function.
# source("generateNetworks.R")

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
# Algorithmic: first index = 100%, second index = 50%
#              (0) => 1,0,0, (0.5) => 0,1,0, (1) => 0,0,1
#              (0.25) => 0.5,0.5,0, (0.75) => 0,0.5,0.5
colorPaletteAlg <- function(intervention) {
  intervention <- as.numeric(strsplit(
    gsub(intervention, pattern = "[(|)]", replacement = ""),
    split = "->")[[1]])
  x <- intervention[1]
  y <- if(length(intervention) == 2) {
    intervention[2]
  } else {
    intervention[1]
  }
  DescTools::CmykToRgb(
    min(max(0, (0.5-x)/0.5) + 0.5*max(0, (0.5-y)/0.5), 1),
    min(max(0, (0.5 - abs(x - 0.5))/0.5)
        + 0.5*max(0, (0.5 - abs(y - 0.5))/0.5), 1),
    min(max(0, (x-0.5)/0.5)+ 0.5*max(0, (y-0.5)/0.5), 1),
    0.25
  )
}

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
### Functions: ################################################################
changeAffinityLevels <- function(df) {
  df |> tidytable::mutate(
    SpeciesAffinity = tidytable::case_when(
      SpeciesAffinity == "rep_0" ~ "100% 0",
      SpeciesAffinity == "evensplit_01" ~ "50% 0, 50% 1",
      SpeciesAffinity == "runif" ~ "Uniform(0, 1)",
      TRUE ~ SpeciesAffinity
    ),
    SpeciesAffinity = factor(SpeciesAffinity, levels = c(
      "100% 0", "50% 0, 50% 1", "Uniform(0, 1)"
    ), ordered = TRUE)
  )
}

changeInterventionLevels <- function(df) {
  df |> tidytable::mutate(
    Intervention = factor(
      Intervention,
      levels = t(interventionMatrix)[1:prod(dim(interventionMatrix))],
      ordered = TRUE
    ),
    InterventionInitial = tidytable::case_when(
      Intervention %in% interventionMatrix[1, ] ~ diag(interventionMatrix)[1],
      Intervention %in% interventionMatrix[2, ] ~ diag(interventionMatrix)[2],
      Intervention %in% interventionMatrix[3, ] ~ diag(interventionMatrix)[3],
      Intervention %in% interventionMatrix[4, ] ~ diag(interventionMatrix)[4],
      Intervention %in% interventionMatrix[5, ] ~ diag(interventionMatrix)[5],
      TRUE ~ NA_character_
    ),
    InterventionInitial = factor(
      InterventionInitial,
      levels = c(
        diag(interventionMatrix)
      ), ordered = TRUE
    ),
    InterventionFinal = tidytable::case_when(
      Intervention %in% interventionMatrix[, 1] ~ diag(interventionMatrix)[1],
      Intervention %in% interventionMatrix[, 2] ~ diag(interventionMatrix)[2],
      Intervention %in% interventionMatrix[, 3] ~ diag(interventionMatrix)[3],
      Intervention %in% interventionMatrix[, 4] ~ diag(interventionMatrix)[4],
      Intervention %in% interventionMatrix[, 5] ~ diag(interventionMatrix)[5],
      TRUE ~ NA_character_
    ),
    InterventionFinal = factor(
      InterventionFinal,
      levels = c(
        diag(interventionMatrix)
      ), ordered = TRUE
    )
  )
}

unifyAffinityBins <- function(., n = 5) {
  tidytable::separate_wider_delim(
    .,
    col = "AffinityBins", names = c("Left", "Right"), delim = ","
  ) |> tidytable::mutate(
    Left =
      round(as.numeric(gsub(pattern = "^[(]", replacement = "", x = Left))*n)/n,
    Right =
      round(as.numeric(gsub(pattern = "\\]$", replacement = "", x = Right))*n)/n,
    AffinityBins = ifelse(
      is.na(Right), as.character(Left),
      paste0("(", Left, ", ", Right, "]")
    )
  )
}

# Verify as we load that the intervention times are calc'd correctedly.