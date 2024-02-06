# Plots used for ABC Presentation.
options(bitmapType = "cairo")

datfolders <- c(
  "TSTS_Simulations_2-1-6_2024-01-19",
  "TSTS_Simulations_2-2-6_2024-01-19",
  "TSTS_Simulations_2-3-6_2024-01-19",
  "TSTS_Simulations_3-1-7_2024-01-22",
  "TSTS_Simulations_3-2-7_2024-01-22",
  "TSTS_Simulations_3-3-7_2024-01-22"
)
poolfolders <- c(
  "Data_2024-01-17",
  "Data_2024-01-18"
)
# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer) # Shading: stackoverflow.com/a/24436825
library(patchwork)

library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R")

# Load Data: ##################################################################

results <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "(Simulation|Result)"), function(x) {
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

diversities <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Diversity"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(c(get(names), "Dir" = x))
  })

presences <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Presence"), function(x) {
    names <- load(x)
    stopifnot(length(names) == 1)
    return(get(names))
  })

poolmats <- lapply(
  dir(poolfolders, full.names = TRUE, pattern = "PoolMats"), function(x) {
    names <- load(x)
    return(c(mget(names), "Dir" = x))
  })

# Convert formats for plotting: ###############################################

#

# Plotting: ###################################################################
# Parameters
# targetpatches <- c(2, 7) # 2 is experimental

### Plot of maximally distant true diversities, not separated by niche. #######
# We use the precomputed diversities. We expect experimental patch above,
# control patch below. Each panel should have two lines, one for the no
# intervention case and one for the pool swap case. We make three plots in
# total, one for each distance.

Div_rounding <- 1

diversitiesPlottable <- do.call(
  rbind,
  diversities %>% lapply(function(d) {
    targetpatches <- d$Ellipsis$Intervention$Patches[3]
    targetpatches <- c(targetpatches, ((targetpatches + 4) %% 10) + 1)

    fileproperties <- strsplit(basename(d$Dir), split = "([_]|[.])")[[1]][3]
    fileproperties <- strsplit(fileproperties, split = "-", fixed = TRUE)[[1]]

    # In this case, timescale is characteristic, not simulation.
    # if (d$Ellipsis$Timescale == "Simulation") {
      # d$Diversities$Time <- d$Diversities$Time / d$Ellipsis$ReactionTime
      # d$Ellipsis$Intervention$Time <-
      #   d$Ellipsis$Intervention$Time / d$Ellipsis$ReactionTime
    # }

    d$Diversities %>% dplyr::filter(
      Environment %in% targetpatches
    ) %>% dplyr::mutate(
      Time = round(Time / Div_rounding) * Div_rounding
    ) %>% dplyr::group_by(
      Time, Environment, Measurement
    ) %>% dplyr::summarise(
      Value = median(Value),
      Count = dplyr::n(),
      .groups = "drop"
    ) %>% dplyr::mutate(
      Space = "Inf",
      Intervention = paste(
        dplyr::case_when(
          fileproperties[1] == 2 ~ "Low",
          fileproperties[1] == 3 ~ "High"
        ),
        dplyr::case_when(
          fileproperties[3] == 6 ~ "High",
          fileproperties[3] == 7 ~ "Low"
        ),
        sep = "_"
      ),
      InterventionSpeed =
        dplyr::case_when(
          fileproperties[2] == 1 ~ "None",
          fileproperties[2] == 2 ~ "Fast",
          fileproperties[2] == 3 ~ "Slow"
        ),
      Control = Environment == targetpatches[2]
    )
  })
)

# Div_interventionchoice <- c("Low_High", "High_Low")[2]
for(Div_interventionchoice in c("Low_High", "High_Low")) {

Div_experiment <- ggplot2::ggplot(
  diversitiesPlottable %>% dplyr::filter(
    Intervention == Div_interventionchoice,
    # Environment %in% as.character(targetpatches[1]),
    Control == FALSE,
    Measurement == "Richness"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = InterventionSpeed,
    linetype = InterventionSpeed
  )
) + ggplot2::geom_line(
  show.legend = FALSE
) + ggplot2::scale_color_manual(
  values = c(
    "None" = "#44546a",
    "Fast" = "#5b9bd5",
    "Slow" = "#ed7d31"
  )
) + ggplot2::ylab(
  paste0("Richness\n(", "Experiment", " Patch)")
) + ggplot2::coord_cartesian(
  ylim = range(diversitiesPlottable %>% dplyr::filter(
    Intervention == Div_interventionchoice, Environment != "Gamma",
    Measurement == "Richness"
  ) %>% dplyr::pull(Value))
# ) + ggplot2::geom_vline( # All values already after intervention...
#   xintercept = # SAFE ONLY IN OUR SPECIFIC CASE
#     (diversities[[1]]$Ellipsis$Intervention$Time /
#        diversities[[1]]$Ellipsis$ReactionTime)
)

Div_control <- ggplot2::ggplot(
  diversitiesPlottable %>% dplyr::filter(
    Intervention == Div_interventionchoice,
    # Environment %in% as.character(targetpatches[2]),
    Control == TRUE,
    Measurement == "Richness"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = InterventionSpeed,
    linetype = InterventionSpeed
  )
) + ggplot2::geom_line(
  show.legend = FALSE
) + ggplot2::scale_color_manual(
  values = c(
    "None" = "#44546a",
    "Fast" = "#5b9bd5",
    "Slow" = "#ed7d31"
  )
) + ggplot2::ylab(
  paste0("Richness\n(", "Control", " Patch)")
) + ggplot2::coord_cartesian(
  ylim = range(diversitiesPlottable %>% dplyr::filter(
    Intervention == Div_interventionchoice, Environment != "Gamma",
    Measurement == "Richness"
  ) %>% dplyr::pull(Value))
# ) + ggplot2::geom_vline(
#   xintercept = # SAFE ONLY IN OUR SPECIFIC CASE
#     (diversities[[1]]$Ellipsis$Intervention$Time /
#        diversities[[1]]$Ellipsis$ReactionTime)
)

ggplot2::ggsave(
  Div_experiment / Div_control,
  filename = paste0("Image-PrConstrain-Div-",
                    Div_interventionchoice,
                    ".png"),
  width = 12, height = 8, units = "cm"
)
}


### Plot of maximally distant true diversities, separated by species type. ####
### and species patch type.
# We can't use the precomputed diversities, which separated by basal consumer.


### Plot of Sampling for Null Intervention. ###################################
##### Preparation: ############################################################
# Parameters
targettimes <- c(min(dplyr::bind_rows(samples)$TimeBase), # Pre-intervent
                 400) # 100 is intervention
timeType <- "Time\nFor Time"
spaceType <- "Space\nFor Time"

samplesByRun <- dplyr::bind_rows(samples) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  ParentRun
  # samples[[1]][[1]] %>% dplyr::ungroup()%>% dplyr::group_by(ParentRun
) %>% dplyr::group_map(
  .f = function(.x, .y) {
    controlPatches <- unique(.x$Patch[.x$PatchType == "Control"])
    if (1 %in% controlPatches && 10 %in% controlPatches) {
      experimentPatches <- unique(.x$Patch[.x$PatchType == "Experiment"])
      experimentPatches <- median(experimentPatches)
      controlPatches <- ((experimentPatches + 4) %% 10) + 1
    } else {
      controlPatches <- median(controlPatches)
      experimentPatches <- ((controlPatches + 4) %% 10) + 1
    }
    targetpatches <- c(experimentPatches, controlPatches)

    fileproperties <- strsplit(.y$ParentRun, split = "-", fixed = TRUE)[[1]]

    # Subset:
    controls <- .x %>% dplyr::filter(
      Patch %in% targetpatches,
      TimeBase %in% targettimes,
      !(Patch == targetpatches[2] & # Not diff time and place.
          TimeBase == targettimes[1]),
      !(Patch == targetpatches[1] & # Not same time and place.
          TimeBase == targettimes[2])
    ) %>% dplyr::mutate(
      Control = "Control",
      Type = ifelse(
        Patch == targetpatches[1], timeType, spaceType
      )#,
      # dplyr::across(
      #   .cols = SamplingAlpha : SamplingAlphaInvasive,
      #   .fns = ~ -.x
      # )
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
    ) %>% dplyr::mutate(
      SamplingNonZeroSpecies =
        paste0(
          strsplit(SamplingNonZeroSpecies, ", ", fixed = TRUE)[[1]][
            as.numeric(strsplit(SamplingNonZeroAbundances, ", ",
                                fixed = TRUE)[[1]]) > 1e-4 # TODO TECHDEBT
            ], collapse = ", "
        )
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
      values_to = "Local Species Richness Gain (vs. Control)"
    ) %>% dplyr::mutate(
      Sampled = !grepl(pattern = "(T)", fixed = TRUE, x = `Species Subset`),
      `Species Subset (Base)` = gsub(pattern = " (T)", replacement = "",
                                     x = `Species Subset`, fixed = TRUE),
      Space = "Inf",
      Intervention = paste(
        dplyr::case_when(
          fileproperties[1] == 2 ~ "Low",
          fileproperties[1] == 3 ~ "High"
        ),
        dplyr::case_when(
          fileproperties[3] == 6 ~ "High",
          fileproperties[3] == 7 ~ "Low"
        ),
        sep = "_"
      ),
      InterventionSpeed =
        dplyr::case_when(
          fileproperties[2] == 1 ~ "None",
          fileproperties[2] == 2 ~ "Fast",
          fileproperties[2] == 3 ~ "Slow"
        )
    )


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

# We want the true niche1 species still present, the true niche2 species that
# have appeared, and the overall change in the number of species in comparison
# to the same patch at the same type but with the alternative treatment.
# (These are what detected in control, not detected in control, and overall
#  are attempting to approximate.)

samplesTRUTH <- dplyr::bind_rows(samples) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  ParentRun
  # samples[[1]][[1]] %>% dplyr::ungroup(
) %>% dplyr::group_modify(
  .f = function(.x, .y) {
    controlPatches <- unique(.x$Patch[.x$PatchType == "Control"])
    if (1 %in% controlPatches && 10 %in% controlPatches) {
      experimentPatches <- unique(.x$Patch[.x$PatchType == "Experiment"])
      experimentPatches <- median(experimentPatches)
      controlPatches <- ((experimentPatches + 4) %% 10) + 1
    } else {
      controlPatches <- median(controlPatches)
      experimentPatches <- ((controlPatches + 4) %% 10) + 1
    }
    targetpatches <- c(experimentPatches, controlPatches)

    fileproperties <- strsplit(.y$ParentRun, split = "-", fixed = TRUE)[[1]]

    truevalues <- .x %>% dplyr::filter(
      Patch == targetpatches[1] & # same time and place.
        TimeBase == targettimes[2]
    ) %>% dplyr::distinct(
      # Stuff to Keep
      Time, Patch, PatchType, TimeBase, Control, TimeActualRow, TimeActual,
      # The bit we're going to operate on.
      SamplingNonZeroSpecies, SamplingNonZeroAbundances
    ) %>% dplyr::mutate(
      Space = "Inf",
      Intervention = paste(
        dplyr::case_when(
          fileproperties[1] == 2 ~ "Low",
          fileproperties[1] == 3 ~ "High"
        ),
        dplyr::case_when(
          fileproperties[3] == 6 ~ "High",
          fileproperties[3] == 7 ~ "Low"
        ),
        sep = "_"
      ),
      InterventionSpeed =
        dplyr::case_when(
          fileproperties[2] == 1 ~ "None",
          fileproperties[2] == 2 ~ "Fast",
          fileproperties[2] == 3 ~ "Slow"
        )
    )


    return(truevalues)
  }
)

# Evaluate the truth by breaking species up by niche.
# Then rearrange by space, and evaluate differences between interventions.
samplesTRUTH <- samplesTRUTH %>% dplyr::ungroup(
) %>% dplyr::group_by(
  ParentRun
) %>% dplyr::group_modify(
  .f = function(.x, .y) {
    speciesIDs <- strsplit(.x$SamplingNonZeroSpecies,
                           split = ", ", fixed = TRUE)
    speciesIDs <- lapply(speciesIDs, as.numeric)

    speciesAbundances <- strsplit(.x$SamplingNonZeroAbundances,
                                  split = ", ", fixed = TRUE)
    speciesAbundances <- lapply(speciesAbundances, as.numeric)

    speciesIDs <- lapply(seq_along(speciesIDs),
                         function(i, id, ab) id[[i]][ab[[i]] > 1e-4],
                         id = speciesIDs, # TODO TECHDEBT, USE RESULTS
                         ab = speciesAbundances)

    # pool <- which(dirname(.y$ParentRun) ==
    #                 unlist(lapply(poolmats, function(p) dirname(p$Dir))))
    # pool <- poolmats[[pool]]$Pool

    # # IN OUR SPECIFIC CASE: (The column especially, but also single row)
    # speciesIDs <- speciesIDs[[1]]
    # niches <- pool$Niche_Cat[match(speciesIDs, pool$ID)]
    # nichesPool <- sort(unique(pool$Niche_Cat))

    # Split the speciesIDs up into their types.
    # niches_split <- lapply(nichesPool, function(niche) {
    #   toString(speciesIDs[niches == niche])
    # })
    # niches_split <- data.frame(niches_split)
    # names(niches_split) <- nichesPool

    .x$SamplingNonZeroSpecies <- toString(speciesIDs[[1]])

    return(cbind(.x#,
                 #niches_split
                 ))

  }
)

samplesTRUTH <- samplesTRUTH %>% dplyr::ungroup(
) %>% dplyr::group_by(
  Intervention
) %>% dplyr::group_modify(
  .f = function(.x, .y) {
    Control <- .x %>% dplyr::filter(
      InterventionSpeed == "None"
    ) %>% dplyr::pull(
      SamplingNonZeroSpecies
      )
    Control <- strsplit(Control, ", ", fixed = TRUE)[[1]]
    Truth <- .x$SamplingNonZeroSpecies

    DetectedInControl <- unlist(lapply(
      strsplit(Truth, ", ", fixed = TRUE), function(ids) {
        toString(ids[ids %in% Control])
      }))
    NotDetectedInControl <- unlist(lapply(
      strsplit(Truth, ", ", fixed = TRUE), function(ids) {
        toString(ids[!ids %in% Control])
      }))

    ControlRichness <- length(Control)

    OverallRichness <- unlist(lapply(
      strsplit(Truth, ", ", fixed = TRUE), function(ids) {
        length(ids) - ControlRichness
      }))
    DetectedInControlRichness <- unlist(lapply(
      strsplit(DetectedInControl, ", ", fixed = TRUE), function(ids) {
        length(ids) - ControlRichness
      }))
    NotDetectedInControlRichness <- unlist(lapply(
      strsplit(NotDetectedInControl, ", ", fixed = TRUE), function(ids) {
        length(ids) - 0
      }))

    return(cbind(.x,
                 "Overall" = OverallRichness,
                 "Detected in Control" = DetectedInControlRichness,
                 "Not Detected in Control" = NotDetectedInControlRichness))
  }
# ) %>% dplyr::mutate(
#   dplyr::across(
#     .cols = c(SamplingNonZeroSpecies, DetectedInControl, NotDetectedInControl),
#     .fns = list(Richness = function(.x) {
#       unlist(lapply(strsplit(.x, split = ", ", fixed = TRUE), length))
#     })
#   ),
#   dplyr::across(
#     .cols = dplyr::ends_with("Richness"),
#     .fns = ~ ifelse(InterventionSpeed == "None", -.x, .x)
#   )
# ) %>% dplyr::group_by(
#   Space, Patch, Time, Intervention, InterventionSpeed
# ) %>% dplyr::summarise(
#   `Overall` = sum(SamplingNonZeroSpecies_Richness),
#   `Detected in Control` = sum(DetectedInControl_Richness),
#   `Not Detected in Control` = sum(NotDetectedInControl_Richness),
#   .groups = "drop"
)

##### Plotting: ###############################################################

# LSR_interventionchoice <- c("Low_High", "High_Low")[1]
# LSR_speedchoice <- c("None", "Fast", "Slow")[1]
ylimits <- c(-11, 11)


# Patch so we don't have to re-run
samplesByRun <- lapply(
  samplesByRun,
  function(x) {
    x %>% dplyr::mutate(
      Type = dplyr::case_when(
        Type == "Space-For-Time" ~ "Space\nFor Time",
        Type == "Time-For-Time" ~ "Time\nFor Time",
        TRUE ~ Type
      )
    )
  }
)

for(LSR_interventionchoice in c("Low_High", "High_Low")) {
  for(LSR_speedchoice in c("None", "Fast", "Slow")) {
    

###### Plot of Local Species Richness Gain Without True Values ################
# LSR_nointervent_notrue <- dplyr::bind_rows(samplesByRun)  %>% dplyr::filter(
#   Intervention == LSR_interventionchoice, Sampled,
#   InterventionSpeed == LSR_speedchoice
# ) %>% dplyr::group_by(
#   `Local Species Richness Gain (vs. Control)`, Type, `Species Subset`
# ) %>% dplyr::mutate(
#   `Local Species Richness Gain (vs. Control)` =
#     `Local Species Richness Gain (vs. Control)` +
#     if(length(`Local Species Richness Gain (vs. Control)`) > 1) {
#       rep(c(-1/6, 1/6), c(
#         floor(length(`Local Species Richness Gain (vs. Control)`)/2),
#         ceiling(length(`Local Species Richness Gain (vs. Control)`)/2)
#       ))
#     } else {
#       0
#     },
#   groupsize = dplyr::n()
# ) %>% ggplot2::ggplot(
#   ggplot2::aes(x = Type,
#                y = `Local Species Richness Gain (vs. Control)`,
#                fill = Type)
#   # ) + ggplot2::geom_violin(
#   #   draw_quantiles = c(
#   #     #0.01, 0.05, # Only 20 samples!
#   #     0.1, 0.25, 0.5, 0.75, 0.9#,
#   #     #0.95, 0.99
#   #   )
# ) + ggplot2::geom_dotplot(
#   binaxis = "y", stackdir = "center", position = "dodge", binwidth = 1/5,
#   dotsize = 1.,
#   show.legend = FALSE
# ) + ggplot2::facet_wrap(
#   . ~ factor(`Species Subset (Base)`, levels = c(
#     "Overall", "Detected in Control", "Not Detected in Control"),
#     ordered = TRUE)
# ) + ggplot2::coord_cartesian(
#   ylim = ylimits, expand = FALSE
# ) + ggplot2::scale_y_continuous(
#   breaks = ylimits[1]:ylimits[2]
# ) + ggplot2::theme(
#   panel.grid.minor.y = ggplot2::element_blank()
# )
# ggplot2::ggsave(LSR_nointervent_notrue,
#                 filename = paste0("Image-PrConstrain-LSR-",
#                                   LSR_interventionchoice, "-",
#                                   LSR_speedchoice,
#                                   "-NoInt-WOTrue.png"),
#                 width = 12, height = 10, units = "cm"
# )
#
# ###### Plot of Local Species Richness Gain With True Values ###################
# LSR_nointervent_true <- LSR_nointervent_notrue + ggplot2::geom_point(
#   data = dplyr::bind_rows(samplesByRun) %>% dplyr::filter(
#    Intervention == LSR_interventionchoice,
#     Sampled == FALSE,
#    InterventionSpeed == LSR_speedchoice
#   ),
#   show.legend = FALSE
# )
# ggplot2::ggsave(LSR_nointervent_true,
#                 filename = paste0("Image-PrConstrain-LSR-",
#                                   LSR_interventionchoice, "-",
#                                   LSR_speedchoice,
#                                   "-NoInt-WTrue.png"),
#                 width = 12, height = 10, units = "cm"
# )

###### Plot of Sampling for True Intervention. ################################

###### Plot of Local Species Richness Gain Without True Values ################
LSR_intervent_notrue <- dplyr::bind_rows(samplesByRun) %>% dplyr::filter(
  Intervention == LSR_interventionchoice, Sampled,
  InterventionSpeed == LSR_speedchoice
) %>% dplyr::group_by(
  `Local Species Richness Gain (vs. Control)`, Type, `Species Subset`
) %>% dplyr::mutate(
  `Local Species Richness Gain (vs. Control)` =
    `Local Species Richness Gain (vs. Control)` +
    if(length(`Local Species Richness Gain (vs. Control)`) > 1) {
      rep(c(-1/6, 1/6), c(
        floor(length(`Local Species Richness Gain (vs. Control)`)/2),
        ceiling(length(`Local Species Richness Gain (vs. Control)`)/2)
      ))
    } else {
      0
    },
  groupsize = dplyr::n()
) %>% ggplot2::ggplot(
  ggplot2::aes(x = Type,
               y = `Local Species Richness Gain (vs. Control)`,
               fill = Type)
  # ) + ggplot2::geom_violin(
  #   draw_quantiles = c(
  #     #0.01, 0.05, # Only 20 samples!
  #     0.1, 0.25, 0.5, 0.75, 0.9#,
  #     #0.95, 0.99
  #   )
) + ggplot2::geom_dotplot(
  binaxis = "y", stackdir = "center", position = "dodge", binwidth = 1/5,
  dotsize = 1.,
  show.legend = FALSE
) + ggplot2::facet_wrap(
  . ~ factor(`Species Subset (Base)`, levels = c(
    "Overall", "Detected in Control", "Not Detected in Control"),
    ordered = TRUE)
) + ggplot2::coord_cartesian(
  ylim = ylimits, expand = FALSE
) + ggplot2::scale_y_continuous(
  breaks = ylimits[1]:ylimits[2]
) + ggplot2::theme(
  panel.grid.minor.y = ggplot2::element_blank()
) + ggplot2::scale_fill_manual(
  values = c("#7f6000", "#b6d7a8")
)

ggplot2::ggsave(LSR_intervent_notrue,
                filename = paste0("Image-PrConstrain-LSR-",
                                  LSR_interventionchoice, "-",
                                  LSR_speedchoice,
                                  "-Int-WOTrue.png"),
                width = 12, height = 10, units = "cm"
)

###### Plot of Local Species Richness Gain With True Values ###################
LSR_intervent_true <- LSR_intervent_notrue + ggplot2::geom_point(
  data = dplyr::bind_rows(samplesByRun) %>% dplyr::filter(
    Intervention == LSR_interventionchoice, Sampled == FALSE,
    InterventionSpeed == LSR_speedchoice
  ),
  show.legend = FALSE
)
ggplot2::ggsave(LSR_intervent_true,
                filename = paste0("Image-PrConstrain-LSR-",
                                  LSR_interventionchoice, "-",
                                  LSR_speedchoice,
                                  "-Int-WTrue.png"),
                width = 12, height = 10, units = "cm"
)

##### Plot of Local Species Richness Gain With Counterfactual #################
LSR_intervent_cf <- LSR_intervent_true + ggplot2::geom_point(
  data = samplesTRUTH %>% dplyr::filter(
    Intervention == LSR_interventionchoice,
    InterventionSpeed == LSR_speedchoice
  ) %>% tidyr::pivot_longer(
    cols = `Overall`:`Not Detected in Control`,
    names_to = "Species Subset (Base)",
    values_to = "Local Species Richness Gain (vs. Control)"
  ),
  ggplot2::aes(x = 1.5),
  show.legend = FALSE, color = "purple", fill = "black", shape = 7, size = 2
)
ggplot2::ggsave(LSR_intervent_cf,
                filename = paste0("Image-PrConstrain-LSR-",
                                  LSR_interventionchoice, "-",
                                  LSR_speedchoice,
                                  "-Int-CF.png"),
                width = 12, height = 10, units = "cm"
)
}}

