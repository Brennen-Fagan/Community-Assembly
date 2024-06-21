# Plots used for ABC Presentation.
options(bitmapType = "cairo")

datfolders <- c(
  "TSTS_Simulations_18-1_6-6_2024-05-23"
)
# Libraries: ##################################################################
library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer) # Shading: stackoverflow.com/a/24436825
# library(patchwork) # Consistently getting invalid layout.pos.row/col.
library(ggpubr) # replace patchwork with ggarrange.

library(RMTRCode2)

source("TimeSpaceAndTimeSeries-0-Functions.R")

# Load Data: ##################################################################

results <- lapply(
  dir(datfolders, full.names = TRUE, pattern = "Result"), function(x) {
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
  dir(datfolders, full.names = TRUE, pattern = "PoolMats"), function(x) {
    names <- load(x)
    return(c(mget(names), "Dir" = x))
  })

# Convert formats for plotting: ###############################################

#

# Plotting: ###################################################################
# Parameters
targetpatches <- c(1, 2) # Both experimental

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
                                             fixed = TRUE)[[1]][4]
      filepropertiesIntervention <- ifelse(
        filepropertiesIntervention == "2", "2 ^ (1 - 2 euclid(m, n))",
        ifelse(filepropertiesIntervention == "3", "10 ^ (1 - 2 euclid(m, n))",
               "NA")
      )
    } else {
      filepropertiesIntervention <- "NA"
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
      # d$Ellipsis$Intervention$Time <-
      #   d$Ellipsis$Intervention$Time / d$Ellipsis$ReactionTime
    }

    d$Diversities %>% dplyr::mutate(
      Time = round(Time / Div_rounding) * Div_rounding
    ) %>% dplyr::group_by(
      Time, Environment, Measurement
    ) %>% dplyr::summarise(
      Value = median(Value),
      Count = dplyr::n(),
      .groups = "drop"
    ) %>% dplyr::mutate(
      Space = filepropertiesSpace,
      ExtirpationProportion = filepropertiesExt,
      Intervention = ifelse(
        length(fileproperties) == 6,
        "(0.5, 0.5) -> (0, 1)", "No Intervention"),
      InterventionIntensity = filepropertiesIntervention,
      ID = d$Dir
    )
  })
)

Div_spacechoice <- "5"

# Richness (local) and Space Specific plot in the first patch.
# Assumes that all data are drawn from the same time/history/system.
Div_experiment <- ggplot2::ggplot(
  diversitiesPlottable %>% dplyr::filter(
    Space == Div_spacechoice,
    Environment %in% as.character(targetpatches[1]),
    Measurement == "Richness"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = Intervention,
    linetype = factor(
      interaction(Intervention, InterventionIntensity),
      levels = rev(unique(interaction(Intervention, InterventionIntensity)))),
    group = ID
  )
) + ggplot2::geom_vline(
  xintercept = # SAFE ONLY IN OUR SPECIFIC CASE
    (diversities[[1]]$Ellipsis$Affinity$TimeIntervention /
       diversities[[1]]$Ellipsis$ReactionTime)
) + ggplot2::geom_line(
  # show.legend = FALSE
) + ggplot2::scale_color_manual(
  values = c(
    "No Intervention" = "#44546a",
    "(0.5, 0.5) -> (0, 1)" = "#ed7d31"
  )
) + ggplot2::ylab(
  paste0("Richness\n(Patch", targetpatches[1], ")")
) + ggplot2::coord_cartesian(
  ylim = range(diversitiesPlottable %>% dplyr::filter(
    Space == Div_spacechoice, Environment != "Gamma",
    Measurement == "Richness"
  ) %>% dplyr::pull(Value))
) + ggplot2::facet_wrap(
  . ~ ExtirpationProportion
) + ggplot2::scale_linetype(
  name = "Intervention+\nIntensity"
)

Div_control <- ggplot2::ggplot(
  diversitiesPlottable %>% dplyr::filter(
    Space == Div_spacechoice,
    Environment %in% as.character(targetpatches[2]),
    Measurement == "Richness"
  ),
  ggplot2::aes(
    x = Time,
    y = Value,
    color = Intervention,
    linetype = factor(
      interaction(Intervention, InterventionIntensity),
      levels = rev(unique(interaction(Intervention, InterventionIntensity)))),
    group = ID
  )
) + ggplot2::geom_vline(
  xintercept = # SAFE ONLY IN OUR SPECIFIC CASE
    (diversities[[1]]$Ellipsis$Affinity$TimeIntervention /
       diversities[[1]]$Ellipsis$ReactionTime)
) + ggplot2::geom_line(
  # show.legend = FALSE
) + ggplot2::scale_color_manual(
  values = c(
    "No Intervention" = "#44546a",
    "(0.5, 0.5) -> (0, 1)" = "#ed7d31"
  )
) + ggplot2::ylab(
  paste0("Richness\n(Patch", targetpatches[2], ")")
) + ggplot2::coord_cartesian(
  ylim = range(diversitiesPlottable %>% dplyr::filter(
    Space == Div_spacechoice, Environment != "Gamma",
    Measurement == "Richness"
  ) %>% dplyr::pull(Value))
) + ggplot2::facet_wrap(
  . ~ ExtirpationProportion
) + ggplot2::scale_linetype(
  name = "Intervention+\nIntensity"
)

ggplot2::ggsave(
  ggpubr::ggarrange(Div_experiment, Div_control,
                    common.legend = TRUE, nrow = 2, legend = "right"),
  filename = paste0("Image-ABC-Div-",
                    Div_spacechoice,
                    ".png"),
  width = 12, height = 8, units = "in"
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
targettimes <- c(min(dplyr::bind_rows(samples)$TimeBase), # Pre-intervent
                 5) # 0 is intervention, 10 is timespan, 50% is symmetry.
timeType <- "Time\nFor Time"
spaceType <- "Space\nFor Time"

samplesByRun <- dplyr::bind_rows(samples) %>% dplyr::ungroup(
) %>% dplyr::group_by(
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
               NA) # TECHDEBT (Being lazy...)
      )

    interventionIntensity <-
      ifelse(length(.yProperties) == 2, "No Intervention",
             ifelse(.yProperties[[3]][4] == "2", "2 ^ (1 - 2 euclid(m, n))",
               ifelse(.yProperties[[3]][4] == "3",  "10 ^ (1 - 2 euclid(m, n))",
                      "NA"
               ))
      )

    extirpation <-
      ifelse(.yProperties[[1]][5] == "2", 1,
             ifelse(.yProperties[[1]][5] == "4", 0.9,
                    ifelse(.yProperties[[1]][5] == "6", 0,
                           NA)))


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
    # ) %>% dplyr::mutate( # Shouldn't be necessary, 6e lines 185 -- 187
    #   SamplingNonZeroSpecies =
    #     paste0(
    #       strsplit(SamplingNonZeroSpecies, ", ", fixed = TRUE)[[1]][
    #         as.numeric(strsplit(SamplingNonZeroAbundances, ", ",
    #                             fixed = TRUE)[[1]]) > 1e-4 # TODO TECHDEBT
    #         ], collapse = ", "
    #     )
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

#TODO I've stopped here for the weekend.

# We want the true niche1 species still present, the true niche2 species that
# have appeared, and the overall change in the number of species in comparison
# to the same patch at the same type but with the alternative treatment.
# (These are what detected in control, not detected in control, and overall
#  are attempting to approximate.)

samplesTRUTH <- dplyr::bind_rows(samples[[1]]) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  ParentRun
  # samples[[1]][[1]] %>% dplyr::ungroup(
) %>% dplyr::group_modify(
  .f = function(.x, .y) {

    .yProperties <- strsplit(basename(as.character(.y$ParentRun)), "-")[[1]]

    truevalues <- .x %>% dplyr::filter(
      Patch == targetpatches[1] & # same time and place.
        TimeBase == targettimes[2]
    ) %>% dplyr::distinct(
      # Stuff to Keep
      Time, Patch, PatchType, TimeBase, Control, TimeActualRow, TimeActual,
      # The bit we're going to operate on.
      SamplingNonZeroSpecies, SamplingNonZeroAbundances
    )

    truevalues$Space <- .yProperties[6]
    truevalues$Intervention <- ifelse(
      .yProperties[length(.yProperties)] == "Intervention.RData",
      "Pool Swap", "No Intervention")

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

    pool <- which(dirname(.y$ParentRun) ==
                    unlist(lapply(poolmats, function(p) dirname(p$Dir))))
    pool <- poolmats[[pool]]$Pool

    # IN OUR SPECIFIC CASE: (The column especially, but also single row)
    speciesIDs <- speciesIDs[[1]]
    niches <- pool$Niche_Cat[match(speciesIDs, pool$ID)]
    nichesPool <- sort(unique(pool$Niche_Cat))

    # Split the speciesIDs up into their types.
    niches_split <- lapply(nichesPool, function(niche) {
      toString(speciesIDs[niches == niche])
    })
    niches_split <- data.frame(niches_split)
    names(niches_split) <- nichesPool

    .x$SamplingNonZeroSpecies <- toString(speciesIDs)

    return(cbind(.x, niches_split))

  }
)

samplesTRUTH <- samplesTRUTH %>% dplyr::ungroup(
) %>% dplyr::mutate(
  dplyr::across(
    .cols = c(SamplingNonZeroSpecies, Niche1, Niche2),
    .fns = list(Richness = function(.x) {
      unlist(lapply(strsplit(.x, split = ", ", fixed = TRUE), length))
    })
  ),
  dplyr::across(
    .cols = dplyr::ends_with("Richness"),
    .fns = ~ ifelse(Intervention == "No Intervention", -.x, .x)
  )
) %>% dplyr::group_by(
  Space, Patch, Time
) %>% dplyr::summarise(
  `Overall` = sum(SamplingNonZeroSpecies_Richness),
  `Detected in Control` = sum(Niche1_Richness),
  `Not Detected in Control` = sum(Niche2_Richness),
  .groups = "drop"
)

##### Plotting: ###############################################################

LSR_spacechoice <- "Inf"
ylimits <- c(-8, 8)
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

###### Plot of Local Species Richness Gain Without True Values ################
LSR_nointervent_notrue <- dplyr::bind_rows(samplesByRun)  %>% dplyr::filter(
  Space == LSR_spacechoice, Intervention == "No Intervention", Sampled
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
ggplot2::ggsave(LSR_nointervent_notrue,
                filename = paste0("Image-ABC-LSR-",
                                  LSR_spacechoice,
                                  "-NoInt-WOTrue.png"),
                width = 12, height = 10, units = "cm"
)

###### Plot of Local Species Richness Gain With True Values ###################
LSR_nointervent_true <- LSR_nointervent_notrue + ggplot2::geom_point(
  data = dplyr::bind_rows(samplesByRun) %>% dplyr::filter(
    Space == LSR_spacechoice, Intervention == "No Intervention",
    Sampled == FALSE
  ),
  show.legend = FALSE
)
ggplot2::ggsave(LSR_nointervent_true,
                filename = paste0("Image-ABC-LSR-",
                                  LSR_spacechoice,
                                  "-NoInt-WTrue.png"),
                width = 12, height = 10, units = "cm"
)

###### Plot of Sampling for True Intervention. ################################

###### Plot of Local Species Richness Gain Without True Values ################
LSR_intervent_notrue <- dplyr::bind_rows(samplesByRun) %>% dplyr::filter(
  Space == LSR_spacechoice, Intervention == "Pool Swap", Sampled
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
                filename = paste0("Image-ABC-LSR-",
                                  LSR_spacechoice,
                                  "-Int-WOTrue.png"),
                width = 12, height = 10, units = "cm"
)

###### Plot of Local Species Richness Gain With True Values ###################
LSR_intervent_true <- LSR_intervent_notrue + ggplot2::geom_point(
  data = dplyr::bind_rows(samplesByRun) %>% dplyr::filter(
    Space == LSR_spacechoice, Intervention == "Pool Swap", Sampled == FALSE
  ),
  show.legend = FALSE
)
ggplot2::ggsave(LSR_intervent_true,
                filename = paste0("Image-ABC-LSR-",
                                  LSR_spacechoice,
                                  "-Int-WTrue.png"),
                width = 12, height = 10, units = "cm"
)

##### Plot of Local Species Richness Gain With Counterfactual #################
LSR_intervent_cf <- LSR_intervent_true + ggplot2::geom_point(
  data = samplesTRUTH %>% dplyr::filter(
    Space == LSR_spacechoice
  ) %>% tidyr::pivot_longer(
    cols = `Overall`:`Not Detected in Control`,
    names_to = "Species Subset (Base)",
    values_to = "Local Species Richness Gain (vs. Control)"
  ),
  ggplot2::aes(x = 1.5),
  show.legend = FALSE, color = "purple", fill = "black", shape = 7, size = 2
)
ggplot2::ggsave(LSR_intervent_cf,
                filename = paste0("Image-ABC-LSR-",
                                  LSR_spacechoice,
                                  "-Int-CF.png"),
                width = 12, height = 10, units = "cm"
)

