# Plots used for BES Macro Presentation.
options(bitmapType = "cairo")

# Currently only set up for one at a time...
datfolders <- c(
  "TSTS_Simulations_18-1_6-6_2024-05-23"
)

targetFiles <- c(
  "18-1-4-15-2_6-6-6.RData", "18-1-4-15-2_6-6-6_12-1-p-3_1-1.RData"
)

targettimes <- c(0,
                 200) # 0 is intervention, ... is timespan, 50% is symmetry.

timeType <- " Time Series " # "Time\nFor Time"
spaceType <- "Space\nFor Time"

histFillColor <- "#70ad47" # green
lineNullColor <- "#5b9bd5" # blue
linePertColor <-  "#ed7d31" # orange

# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer) # Shading: stackoverflow.com/a/24436825
library(ggpubr) # replace patchwork with ggarrange.

library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R")

# Load Data: ##################################################################

samples <- unlist(lapply(targetFiles, function(tf)
  dir(datfolders, full.names = TRUE,
               pattern = paste0("Sampling_", tf))
))
samples <- lapply(
  samples, function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

diversities <- unlist(lapply(targetFiles, function(tf)
  dir(datfolders, full.names = TRUE,
      pattern = paste0("Diversity_", tf))
))
diversities <- lapply(
  diversities, function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = x))
  })

# poolmats <- lapply(
#   dir(datfolders, full.names = TRUE, pattern = "PoolPatchDynamics"), function(x) {
#     names <- load(x)
#     return(c(mget(names), "Dir" = x))
#   })

# Convert formats for plotting: ###############################################

#

# Plotting: ###################################################################
# Parameters
targetpatches <- c(2, 1) # first patch is focal, second patch is control

### Plot of maximally distant true diversities, not separated by niche. #######
# We use the precomputed diversities. We expect experimental patch above,
# control patch below. Each panel should have two lines, one for the no
# intervention case and one for the pool swap case. We make three plots in
# total, one for each distance.

Div_rounding <- 1

diversitiesPlottable <- do.call(
  rbind,
  diversities %>% lapply(function(d) {
    fileproperties <- strsplit(basename(d$Dir), split = "_", fixed = TRUE)[[1]]

    filepropertiesInit <- strsplit(fileproperties[3], split = "-",
                                   fixed = TRUE)[[1]]

    filepropertiesExt <- filepropertiesInit[3]
    filepropertiesExt <- ifelse(
      filepropertiesExt == "2", 1,
      ifelse(filepropertiesExt == "4", 0.9,
             ifelse(filepropertiesExt == "6", 0,
                    NA))
    )

    filepropertiesSpace <- filepropertiesInit[4]
    filepropertiesSpace <- ifelse(
      filepropertiesSpace == "10", 0,
      ifelse(filepropertiesSpace == "15", 5,
             ifelse(filepropertiesSpace == "NA", "Inf", NA))
    )

    if(length(fileproperties) == 6) {
      filepropertiesIntervention <- strsplit(fileproperties[5], split = "-",
                                             fixed = TRUE)[[1]]
      filepropertiesInterventionType <- filepropertiesIntervention[1]
      filepropertiesIntervention <- filepropertiesIntervention[4]
      filepropertiesIntervention <- ifelse(
        filepropertiesIntervention == "2", "2 ^ (1 - 2 euclid(m, n))",
        ifelse(filepropertiesIntervention == "3", "10 ^ (1 - 2 euclid(m, n))",
               "NA")
      )
    } else {
      filepropertiesIntervention <- "NA"
      filepropertiesInterventionType <- "NA"
    }

    # Need to unify Diversities Format
    d$Diversities$alpha <- d$Diversities$alpha %>% dplyr::select(
      -Species
    ) %>% tidyr::pivot_longer(
      cols = Richness:Richness_Consumer,
      names_to = "Aggregation", values_to = "Measurement"
    ) %>% dplyr::mutate(
      Environment = as.character(Environment)
    )

    d$Diversities$beta <- dplyr::bind_rows(
      d$Diversities$beta
    ) %>% dplyr::mutate(
      Environment = paste(Env1, Env2),
      Aggregation = "Jaccard",
      Jaccard = as.double(Jaccard)
    ) %>% dplyr::rename(
      Measurement = Jaccard
    ) %>% dplyr::select(
      -Env1, -Env2
    )

    d$Diversities$gamma <- d$Diversities$gamma %>% dplyr::mutate(
      Aggregation = ifelse(Aggregation == "Gamma", "Richness", Aggregation)
    ) %>% dplyr::rename(
      Agg = Aggregation,
      `_Gamma` = Richness,
      `_Gamma_Basal` = Basals,
      `_Gamma_Consumer` = Consumers
    ) %>% tidyr::pivot_longer(
      cols = `_Gamma`:`_Gamma_Consumer`,
      names_to = "Agg2", values_to = "Measurement"
    ) %>% dplyr::mutate(
      Environment = NA,
      Aggregation = paste0(Agg, Agg2)
    ) %>% dplyr::select(
      -Agg, -Agg2
    )

    d$Diversities <- dplyr::bind_rows(
      d$Diversities
    ) %>% dplyr::rename(
      Value = Measurement,
      Measurement = Aggregation
    )

    if (d$Ellipsis$Timescale == "Simulation") {
      d$Diversities$Time <- d$Diversities$Time / d$Ellipsis$ReactionTime
    }

    d$Diversities %>% dplyr::mutate(
      Time = round(Time / Div_rounding) * Div_rounding
    ) %>% dplyr::group_by(
      Time, Environment, Measurement
    ) %>% dplyr::summarise(
      Value = median(Value),
      Frequency = dplyr::n(),
      .groups = "drop"
    ) %>% dplyr::mutate(
      Space = filepropertiesSpace,
      ExtirpationProportion = filepropertiesExt,
      Intervention = ifelse(
        length(fileproperties) == 6,
        ifelse(filepropertiesInterventionType=="40",
               "(0.5, 0.5) -> (0, 1)",
               ifelse(filepropertiesInterventionType=="12",
                      "(0.5, 0.5) -> (0.5, 1)",
                      "UNLISTED")),
        "No Intervention"),
      InterventionIntensity = filepropertiesIntervention,
      ID = d$Dir
    )
  })
)


yrange <- range(diversitiesPlottable %>% dplyr::filter(
  # Space == Div_spacechoice,
  Environment != "Gamma",
  Measurement == "Richness"
) %>% dplyr::pull(Value))
plotDiversityInPatch <- function(data, yrange) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Time,
      y = Value,
      color = Intervention,
      linetype = factor(
        interaction(Intervention, InterventionIntensity),
        levels = rev(unique(interaction(Intervention, InterventionIntensity)))),
      group = ID
    )
  ) + ggplot2::geom_line(
    show.legend = FALSE, linewidth = 2
  ) + ggplot2::scale_color_manual(
    values = c(
      "No Intervention" = lineNullColor,
      "(0.5, 0.5) -> (0, 1)" = linePertColor,
      "(0.5, 0.5) -> (0.5, 1)" = linePertColor
    )
  ) + ggplot2::ylab(
    "Richness"
  ) + ggplot2::coord_cartesian(
    ylim =
    # ) + ggplot2::facet_wrap(
    #   . ~ ExtirpationProportion
  ) + ggplot2::scale_linetype(
    name = "Intervention+\nIntensity"
  ) + ggplot2::theme_bw(
  ) + ggplot2::theme(
    text = ggplot2::element_text(size = 30)
  )
}

# Richness (local) and Space Specific plot in the first patch.
# Assumes that all data are drawn from the same time/history/system.
Div_experiment_before <- plotDiversityInPatch(
  diversitiesPlottable %>% dplyr::filter(
    # Space == Div_spacechoice,
    Environment %in% as.character(targetpatches[1]),
    Measurement == "Richness",
    Intervention %in% "No Intervention"#,
    # InterventionIntensity %in% Div_intenschoice[1],
    # ExtirpationProportion %in% Div_extprpchoice
  ), yrange)

Div_experiment <- plotDiversityInPatch(
  diversitiesPlottable %>% dplyr::filter(
    # Space == Div_spacechoice,
    Environment %in% as.character(targetpatches[1]),
    Measurement == "Richness"#,
    # Intervention %in% Div_interchoice,
    # InterventionIntensity %in% Div_intenschoice,
    # ExtirpationProportion %in% Div_extprpchoice
  ), yrange)

Div_control_before <- plotDiversityInPatch(
  diversitiesPlottable %>% dplyr::filter(
    # Space == Div_spacechoice,
    Environment %in% as.character(targetpatches[2]),
    Measurement == "Richness",
    Intervention %in% "No Intervention"#,
    # InterventionIntensity %in% Div_intenschoice[1],
    # ExtirpationProportion %in% Div_extprpchoice
  ), yrange)

Div_control <- plotDiversityInPatch(
  diversitiesPlottable %>% dplyr::filter(
    # Space == Div_spacechoice,
    Environment %in% as.character(targetpatches[2]),
    Measurement == "Richness"#,
    # Intervention %in% Div_interchoice,
    # InterventionIntensity %in% Div_intenschoice,
    # ExtirpationProportion %in% Div_extprpchoice
  ), yrange)

ggplot2::ggsave(
  ggpubr::ggarrange(
    Div_experiment_before,
    Div_control_before,
    common.legend = TRUE, ncol = 2, legend = "right"),
  filename = paste0("Image-BESMacro-Div-",
                    targetFiles[1],
                    ".png"),
  width = 16, height = 4, units = "in"
)
ggplot2::ggsave(
  ggpubr::ggarrange(
    Div_experiment,
    Div_control,
    common.legend = TRUE, ncol = 2, legend = "right"),
  filename = paste0("Image-BESMacro-Div-",
                    targetFiles[2],
                    ".png"),
  width = 16, height = 4, units = "in"
)

### Plot of maximally distant true diversities, separated by species type. ####
### and species patch type.
# We can't use the precomputed diversities, which separated by basal consumer.

### Plot of Sampling for Null Intervention. ###################################
##### Preparation: ############################################################
# Originally developed to show contrast between a control patch and an
# intervention patch.
# The focal patch is targetpatches[1] and focal time is targettimes[2].
# Hence the controls are at targetpatches[1] targettimes[1] or 2 and 2.
# Parameters
# interventchoice <- c("No Intervention", "(0.5, 0.5) -> (0.5, 1)")
# intensitychoice <- c("No Intervention",
#                      # "2 ^ (1 - 2 euclid(m, n))")
#                      "10 ^ (1 - 2 euclid(m, n))")

samplesByRun <- dplyr::bind_rows(samples) %>% dplyr::group_by(
  ParentRun
  # samples[[1]][[1]] %>% dplyr::ungroup()%>% dplyr::group_by(ParentRun # debug
) %>% dplyr::group_map(
  .f = function(.x, .y) {

    .yProperties <- strsplit(unlist(strsplit(
      basename(as.character(.y)), "_")), "-"
    )

    space <- if((length(.yProperties) > 2 && .yProperties[[3]][3] == "p")
                || length(.yProperties) == 2) {
      .yProperties[[1]][4]
    } else {
      .yProperties[[3]][3] # More accurately, this means a spatial intervention
    }
    space <- ifelse(
      space == "10", 0,
      ifelse(space == "15", 5,
             ifelse(space == "NA", Inf,
                    NA))
    )

    intervention <-
      ifelse(length(.yProperties) == 2, "No Intervention",
             ifelse(
               (length(.yProperties) > 2 && .yProperties[[3]][1] == "40"),
               "(0.5, 0.5) -> (0, 1)",
               ifelse(length(.yProperties) > 2 && .yProperties[[3]][1] == "12",
                      "(0.5, 0.5) -> (0.5, 1)",
                      NA)) # TECHDEBT (Being lazy...)
      )

    interventionIntensity <-
      ifelse(length(.yProperties) == 2, "No Intervention",
             ifelse(.yProperties[[3]][4] == "2", "2 ^ (1 - 2 euclid(m, n))",
                    ifelse(.yProperties[[3]][4] == "3",  "10 ^ (1 - 2 euclid(m, n))",
                           "NA"
                    ))
      )

    extirpation <-
      ifelse(.yProperties[[1]][3] == "2", 1,
             ifelse(.yProperties[[1]][3] == "4", 0.9,
                    ifelse(.yProperties[[1]][3] == "6", 0,
                           NA)))


    # Subset:
    controls <- .x %>% dplyr::filter(
      Patch %in% targetpatches,
      TimeBase %in% targettimes,
      !(Patch == targetpatches[2] & # Not [[diff time] and [diff place]].
          TimeBase == targettimes[1]),
      !(Patch == targetpatches[1] & # Not [[same time] and [same place]].
          TimeBase == targettimes[2])
    ) %>% dplyr::mutate(
      Control = "Control",
      Type = ifelse(
        Patch == targetpatches[1], timeType, spaceType
      )
    )

    experiments <- .x %>% dplyr::filter(
      Patch == targetpatches[1] & # same time and place.
        TimeBase == targettimes[2]
    ) %>% dplyr::mutate(
      Control = "Experiment"
    )

    experiments <- rbind(
      experiments %>% dplyr::mutate(Type = timeType),
      experiments %>% dplyr::mutate(Type = spaceType)
    )

    combined <- rbind(experiments, controls) %>% dplyr::rowwise(
    ) %>% dplyr::ungroup(
    ) %>% dplyr::group_by(
      Type, SamplingRun
    ) %>% dplyr::group_modify(
      .f = ~ computeSpeciesInControl(.x, Time = "TimeBase")
    ) %>% dplyr::group_modify(
      .f = ~ computeSpeciesInControl(.x, Time = "TimeBase",
                                     IDColumn = "SamplingNonZeroSpecies",
                                     OutPrefix = "True")
    ) %>% dplyr::ungroup(
    ) %>% dplyr::group_by(
      Type, SamplingRun
    ) %>% dplyr::mutate(
      nControl = sum(Control == "Control"),
      TrueAlpha = TrueAlphaNative + TrueAlphaInvasive,
      dplyr::across(
        .cols = c(SamplingAlpha, SamplingAlphaNative, SamplingAlphaInvasive,
                  TrueAlpha, TrueAlphaNative, TrueAlphaInvasive),
        .fns = ~ ifelse(Control == "Control", -.x, .x)
      )
    ) %>% dplyr::summarise(
      DeltaSamplingAlpha = sum(SamplingAlpha),
      DeltaSamplingAlphaNative = sum(SamplingAlphaNative),
      DeltaSamplingAlphaInvasive = sum(SamplingAlphaInvasive),
      DeltaTrueAlpha = sum(TrueAlpha),
      DeltaTrueAlphaNative = sum(TrueAlphaNative),
      DeltaTrueAlphaInvasive = sum(TrueAlphaInvasive),
      .groups = "drop"
    ) %>% dplyr::rename(
      `Overall` = DeltaSamplingAlpha,
      `Detected in Control` = DeltaSamplingAlphaNative,
      `Not Detected in Control` = DeltaSamplingAlphaInvasive,
      `Overall (T)` = DeltaTrueAlpha,
      `Detected in Control (T)` = DeltaTrueAlphaNative,
      `Not Detected in Control (T)` = DeltaTrueAlphaInvasive
    ) %>% tidyr::pivot_longer(
      cols = c(`Overall`, `Detected in Control`, `Not Detected in Control`,
               `Overall (T)`, `Detected in Control (T)`, `Not Detected in Control (T)`),
      names_to = "Species Subset",
      values_to = "Local Species Richness Change (vs. Control)"
    ) %>% dplyr::mutate(
      Sampled = !grepl(pattern = "(T)", fixed = TRUE, x = `Species Subset`),
      `Species Subset (Base)` = gsub(pattern = " (T)", replacement = "",
                                     x = `Species Subset`, fixed = TRUE)
    )

    combined$Space <- space
    combined$Intervention <- intervention
    combined$InterventionIntensity <- interventionIntensity
    combined$Extirpation <- extirpation

    return(combined)
  }
)

# Ok, this gives me per file:
#   Space - For - Time, With Sampling Error
#   Space - For - Time, Without Sampling Error
#   Time - For - Time, With Sampling Error
#   Time - For - Time, Without Sampling Error
##### But now we want to add on the counterfactual ############################
##### which requires comparing between files.

# PREVIOUSLY:
# We want the true niche1 species still present, the true niche2 species that
# have appeared, and the overall change in the number of species in comparison
# to the same patch at the same time but with the alternative treatment.
# (These are what detected in control, not detected in control, and overall
#  are attempting to approximate.)

# BUT NOW:
# We're NOT using niches in this way.
# Instead, we're adopting the no intervention case as the counterfactual
# baseline. Hence, there is 0 change overall, 0 change in species detected
# in the control, and 0 species detected out of control when looking at the
# no intervention case.
# Instead, we want to know the effect of the intervention on each metric.
# The sampling methods are then trying to capture this (previous sentence).

# Make sure to associate each run with the correct no intervention run.
samplesTRUTH <- dplyr::bind_rows(samples) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  ParentRun
  # samples[[1]][[1]] %>% dplyr::ungroup()%>% dplyr::group_by(ParentRun # debug
) %>% dplyr::group_modify(
  .f = function(.x, .y) {

    .yProperties <- strsplit(basename(as.character(.y)), "_")

    # ParentRun is a bit of a misnomer here, since it's actually the
    # Intervention or Simulation run, but we strictly want the Simulations.
    # This is stored in the first pair of .yProperties at this stage.
    initialRun <- paste0(.yProperties[[1]][1:2], collapse = "_")

    .yProperties <- strsplit(unlist(.yProperties), "-")

    space <- if((length(.yProperties) > 2 && .yProperties[[3]][3] == "p")
                || length(.yProperties) == 2) {
      .yProperties[[1]][4]
    } else {
      .yProperties[[3]][3] # More accurately, this means a spatial intervention
    }
    space <- ifelse(
      space == "10", 0,
      ifelse(space == "15", 5,
             ifelse(space == "NA", Inf,
                    NA))
    )

    intervention <-
      ifelse(length(.yProperties) == 2, "No Intervention",
             ifelse(
               (length(.yProperties) > 2 && .yProperties[[3]][1] == "40"),
               "(0.5, 0.5) -> (0, 1)",
               ifelse(length(.yProperties) > 2 && .yProperties[[3]][1] == "12",
                      "(0.5, 0.5) -> (0.5, 1)",
                      NA)) # TECHDEBT (Being lazy...)
      )

    interventionIntensity <-
      ifelse(length(.yProperties) == 2, "No Intervention",
             ifelse(.yProperties[[3]][4] == "2", "2 ^ (1 - 2 euclid(m, n))",
                    ifelse(.yProperties[[3]][4] == "3",  "10 ^ (1 - 2 euclid(m, n))",
                           "NA"
                    ))
      )

    extirpation <-
      ifelse(.yProperties[[1]][3] == "2", 1,
             ifelse(.yProperties[[1]][3] == "4", 0.9,
                    ifelse(.yProperties[[1]][3] == "6", 0,
                           NA)))

    truevalues <- .x %>% dplyr::filter(
      Patch == targetpatches[1] & # same time and place.
        TimeBase == targettimes[2]
    ) %>% dplyr::distinct(
      # Stuff to Keep
      Time, Patch, PatchType, TimeBase, Control, TimeActualRow, TimeActual,
      # The bit we're going to operate on.
      SamplingNonZeroSpecies, SamplingNonZeroAbundances
    )

    truevalues$Space <- space
    truevalues$Intervention <- intervention
    truevalues$InterventionIntensity <- interventionIntensity
    truevalues$Extirpation <- extirpation
    truevalues$InitialRun <- initialRun

    return(truevalues)
  }
)

# Will be useful if we perform "niche" (trait, w/e) evaluation as well.
# samplesTRUTH <- samplesTRUTH %>% dplyr::ungroup(
# ) %>% dplyr::group_by(
#   InitialRun
# ) %>% dplyr::group_modify(
#   .f = function(.x, .y) {
#     speciesIDs <- strsplit(.x$SamplingNonZeroSpecies,
#                            split = ", ", fixed = TRUE)
#     speciesIDs <- lapply(speciesIDs, as.numeric)
#
#     # speciesAbundances <- strsplit(.x$SamplingNonZeroAbundances,
#     #                               split = ", ", fixed = TRUE)
#     # speciesAbundances <- lapply(speciesAbundances, as.numeric)
#     #
#     # speciesIDs <- lapply(seq_along(speciesIDs),
#     #                      function(i, id, ab) id[[i]][ab[[i]] > 1e-4],
#     #                      id = speciesIDs, # TODO TECHDEBT, USE RESULTS
#     #                      ab = speciesAbundances)
#
#     # pool <- which(dirname(.y$InitialRun) ==
#     #                 unlist(lapply(poolmats, function(p) dirname(p$Dir))))
#     # pool <- poolmats[[pool]]$Pool
#     pool <- poolmats[[1]]$Pool
#
#     # IN OUR SPECIFIC CASE: (The column especially, but also single row)
#     # speciesIDsAll <- sort(unique(unlist(speciesIDs)))
#     # niches <- pool$Affinity[match(speciesIDs, pool$ID)]
#     niches <- lapply(lapply(speciesIDs, match, table = pool$ID), # get position
#                      function(position) pool$Affinity[position]) # get affinity
#     nichesPool <- sort(unique(pool$Affinity))
#
#     # Split the speciesIDs up into their types.
#     niches_split <- lapply(seq_along(niches), function(nid) {
#       splits <- lapply(nichesPool, function(niche) {
#         toString(speciesIDs[[nid]][niches[[nid]] == niche])
#       })
#       splits <- data.frame(splits)
#       colnames(splits) <- paste0("Niche", nichesPool)
#       return(splits)
#     }) %>% dplyr::bind_rows()
#
#     # .x$SamplingNonZeroSpecies <- toString(speciesIDs)
#
#     return(cbind(.x, niches_split))
#   }
# )

# samplesTRUTH <- samplesTRUTH %>% dplyr::ungroup(
# ) %>% dplyr::mutate(
#   dplyr::across(
#     .cols = c(SamplingNonZeroSpecies, tidyr::starts_with("Niche")),
#     .fns = list(Richness = function(.x) {
#       unlist(lapply(strsplit(.x, split = ", ", fixed = TRUE), length))
#     })
#   ),
#   dplyr::across(
#     .cols = dplyr::ends_with("Richness"),
#     .fns = ~ ifelse(Intervention == "No Intervention", -.x, .x)
#   )
# ) %>% dplyr::group_by(
#   Space, Patch, Time
# ) %>% dplyr::summarise(
#   `Overall` = sum(SamplingNonZeroSpecies_Richness),
#   `Detected in Control` = sum(Niche1_Richness),
#   `Not Detected in Control` = sum(Niche2_Richness),
#   .groups = "drop"
# )

samplesTRUTH <- samplesTRUTH %>% dplyr::ungroup(
) %>% dplyr::group_by(
  InitialRun
) %>% dplyr::group_modify(
  .f = function(.x, .y) {
    speciesIDs <- strsplit(.x$SamplingNonZeroSpecies,
                           split = ", ", fixed = TRUE)
    speciesIDs <- lapply(speciesIDs, as.numeric)

    baseline <- which(.y$InitialRun == .x$ParentRun)
    if(length(baseline) == 0) {
      stop(paste0("Missing Baseline:", .y$InitialRun))
    }

    inControl <- lapply(speciesIDs, function(ids) {
      length(ids[ids %in% speciesIDs[[baseline]]]) - length(speciesIDs[[baseline]])
    })
    notInControl <- lapply(speciesIDs, function(ids) {
      length(ids[!ids %in% speciesIDs[[baseline]]]) - 0
    })
    overallChange <- lapply(speciesIDs, function(ids) {
      length(ids) - length(speciesIDs[[baseline]])
    })

    differences <- data.frame(
      `Overall` = unlist(overallChange),
      `Detected in Control` = unlist(inControl),
      `Not Detected in Control` = unlist(notInControl)
    )

    return(cbind(.x, differences))
  }
)

##### Plotting: ###############################################################
samplesTRUTH <- samplesTRUTH %>% dplyr::rename(
  `Detected in Control` = Detected.in.Control,
  `Not Detected in Control` = Not.Detected.in.Control
)

###### Master Plot: ###########################################################

histogramData <- dplyr::bind_rows(samplesByRun)  %>% dplyr::filter(
  Sampled, `Species Subset (Base)` %in% c("Overall"),
  !InterventionIntensity %in% c("2 ^ (1 - 2 euclid(m, n))")
  # "No Intervention" "2 ^ (1 - 2 euclid(m, n))" "10 ^ (1 - 2 euclid(m, n))"
) %>% dplyr::group_by(
  `Local Species Richness Change (vs. Control)`,
  Space, Extirpation, Type, `Species Subset (Base)`,
  Intervention, InterventionIntensity
) %>% dplyr::summarise(
  Frequency = dplyr::n(),
  .groups = "drop"
)

factualData <- dplyr::bind_rows(samplesByRun) %>% dplyr::filter(
  SamplingRun > 20, # Temporary patch to deal with my duplicating things.
  # 1st 20 samples are duplicated only and only in some cases...
  Sampled == FALSE, `Species Subset (Base)` %in% c("Overall"),
  !InterventionIntensity %in% c("2 ^ (1 - 2 euclid(m, n))")
) %>% dplyr::select(-SamplingRun) %>% dplyr::distinct()

counterData <- samplesTRUTH %>% dplyr::filter(
) %>% tidyr::pivot_longer(
  cols = `Overall`:`Not Detected in Control`,
  names_to = "Species Subset (Base)",
  values_to = "Local Species Richness Change (vs. Control)"
) %>% dplyr::filter(
  `Species Subset (Base)` %in% c("Overall"),
  !InterventionIntensity %in% c("2 ^ (1 - 2 euclid(m, n))")
)

limitsX <- range(
  c(histogramData$`Local Species Richness Change (vs. Control)`,
    factualData$`Local Species Richness Change (vs. Control)`,
    counterData$`Local Species Richness Change (vs. Control)`)
)
limitsY <- range(histogramData$Frequency)

baseplot <- ggplot2::ggplot(
  histogramData %>% dplyr::filter(
    InterventionIntensity %in% "No Intervention"
  ),
  ggplot2::aes(x = `Local Species Richness Change (vs. Control)`,
               group = `Local Species Richness Change (vs. Control)`,
               y = Frequency)
) + ggplot2::geom_vline(
  xintercept = 0
) + ggplot2::coord_flip(
  xlim = limitsX, ylim = limitsY
) + ggplot2::facet_grid(
  # interaction(Intervention, InterventionIntensity) + Extirpation ~
  # factor(`Species Subset (Base)`, levels = c(
  #   "Overall", "Detected in Control", "Not Detected in Control"),
  #   ordered = TRUE) #+ Type
  .
  ~ Type, switch = "x"
) + ggplot2::theme_bw(
  base_size = 17
) + ggplot2::scale_x_continuous(
  breaks = limitsX[1]:limitsX[2]
) + ggplot2::labs(
  x = "Local Species Richness Change\n(vs. Control)"
) + ggplot2::theme(
  panel.grid = ggplot2::element_blank()
)

baseplot_factual <- baseplot + ggplot2::geom_vline(
  data = factualData %>% dplyr::filter(
    InterventionIntensity %in% "No Intervention"
  ),
  ggplot2::aes(xintercept = `Local Species Richness Change (vs. Control)`),
  color = lineNullColor,
  linetype = "dashed")

baseplot_sampled <- baseplot_factual + ggplot2::geom_bar(
  stat = "identity", position = "dodge", fill = histFillColor
)

ggplot2::ggsave(
  baseplot,
  filename = paste0("histogram_base_TSTS_",targetFiles[1],".png"),
  units = "cm", height = 11, width = 23
)
ggplot2::ggsave(
  baseplot_factual,
  filename = paste0("histogram_factual_TSTS_",targetFiles[1],".png"),
  units = "cm", height = 11, width = 23
)
ggplot2::ggsave(
  baseplot_sampled,
  filename = paste0("histogram_sampled_TSTS_",targetFiles[1],".png"),
  units = "cm", height = 11, width = 23
)

interventionplot <- ggplot2::ggplot(
  histogramData %>% dplyr::filter(
    !InterventionIntensity %in% "No Intervention"
  ),
  ggplot2::aes(x = `Local Species Richness Change (vs. Control)`,
               group = `Local Species Richness Change (vs. Control)`,
               y = Frequency)
) + ggplot2::geom_vline(
  xintercept = 0
) + ggplot2::coord_flip(
  xlim = limitsX, ylim = limitsY
) + ggplot2::facet_grid(
  # interaction(Intervention, InterventionIntensity) + Extirpation ~
  # factor(`Species Subset (Base)`, levels = c(
  #   "Overall", "Detected in Control", "Not Detected in Control"),
  #   ordered = TRUE) #+ Type
  .
  ~ Type, switch = "x"
) + ggplot2::theme_bw(
  base_size = 17
) + ggplot2::scale_x_continuous(
  breaks = limitsX[1]:limitsX[2]
) + ggplot2::labs(
  x = "Local Species Richness Change\n(vs. Control)"
) + ggplot2::theme(
  panel.grid = ggplot2::element_blank()
)

interventionplot_factual <- interventionplot + ggplot2::geom_vline(
  data = factualData %>% dplyr::filter(
    !InterventionIntensity %in% "No Intervention"
  ),
  ggplot2::aes(xintercept = `Local Species Richness Change (vs. Control)`),
  color = lineNullColor,
  linetype = "dashed"
)

interventionplot_sampled <- interventionplot_factual + ggplot2::geom_bar(
  stat = "identity", position = "dodge", fill = histFillColor
)

interventionplot_counter <- interventionplot_sampled + ggplot2::geom_vline(
  data = counterData %>% dplyr::filter(
    !InterventionIntensity %in% "No Intervention"
  ),
  ggplot2::aes(xintercept = `Local Species Richness Change (vs. Control)`),
  color = linePertColor
)

ggplot2::ggsave(
  interventionplot,
  filename = paste0("histogram_base_TSTS_", targetFiles[2], ".png"),
  units = "cm", height = 11, width = 23
)
ggplot2::ggsave(
  interventionplot_factual,
  filename = paste0("histogram_factual_TSTS_", targetFiles[2], ".png"),
  units = "cm", height = 11, width = 23
)
ggplot2::ggsave(
  interventionplot_sampled,
  filename = paste0("histogram_sampled_TSTS_", targetFiles[2], ".png"),
  units = "cm", height = 11, width = 23
)
ggplot2::ggsave(
  interventionplot_counter,
  filename = paste0("histogram_counter_TSTS_", targetFiles[2], ".png"),
  units = "cm", height = 11, width = 23
)

