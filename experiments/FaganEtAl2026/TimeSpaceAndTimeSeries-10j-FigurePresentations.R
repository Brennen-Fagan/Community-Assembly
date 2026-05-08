# Presentation variations:

variant <- c("Networks")[1]
figures <- c(7)
# 2:8 # Available if variant == "Networks"
pref <- c("50% 0, 50% 1", "Uniform(0, 1)")[1]
# Break up the preferences to save memory for Networks 5.
prefstring <- switch(
  pref,
  "100% 0" = "1000",
  "50% 0, 50% 1" = "5050",
  "Uniform(0, 1)" = "Unif"
)

if (variant == "Networks") {
  ### Networks oriented: ######################################################
  #### Setup: #################################################################
  ##### Resources: ############################################################
  # Data:
  source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
  source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
  source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")
  source("TimeSpaceAndTimeSeries-10i-PreparationsTimeBC.R")

  load("TSTS_Interventions_10a1.RData")
  InterventionTimes <- InterventionTimes |> tidytable::select(
    -tidytable::starts_with("Intervention")
  ) # Only for the base cases, so Intervention information all NA.

  # Functions:
  source(file.path("R", "flattenDiversity.R")) # Req'd by below
  source(file.path("R", "generateNetworks.R")) # To create inset graphs.
  library(patchwork) # Plot assembly
  library(ggExtra) # N6 marginal distributions

  # Directories:
  dirImages <- file.path(".", "TSTS_Images_Networks")
  if (!dir.exists(dirImages)) {
    dir.create(dirImages, showWarnings = FALSE)
  }

  ##### Data Management: ######################################################
  figureNetworks <- list(
    graph = list(
      # Note that, unusually, we need all max(time), but only seed for others.
      seed = "2", # "11", "17", "2"!,                 # for examples
      time = c(100, 2000, 25000),                     # for examples
      timeInterventions = c(0, 10, 100, 1000, 10000), # for examples
      pref = c("100% 0", "Uniform(0, 1)"),            # for KDEs
      interventions = c("(0)", "(0.5)", "(1)",        # for KDEs
                        "(0)->(0.5)", "(0.5)->(0)")   # for examples
    ),
    interventions = c("(0)", "(0)->(0.5)", "(0.5)->(0)",
                      "(0.5)", "(0.5)->(1)", "(1)"),
    ci = 0.75
  )

  # Apply initial specification set-up, but then identify practical times
  # post-intervention rather than the approximates supplied above for the
  # examples.
  figureNetworks$graph$specification <-
    diversitiesRichness |> tidytable::select(c(
      # Which network:
      "Time", "Environment1",
      # Which File (Base):
      "PoolPatch", "PoolPatchSeed", "Interactions", "InteractionsSeed",
      "Events", "EventsSeed", "InitialConditions", "InitialConditionsSeed",
      "Dispersal", "NicheDistance",
      # Oops, there was a collision causing human readable to replace machine.
      # Will be replaced SpeciesAffinity#2 with -> SpeciesPreferences.
      "SpeciesAffinity",
      "SpeciesAffinitySeed", "PatchAffinity", "PatchAffinitySeed",
      # Which File (Intervention):
      "InterventionPatchType", "InterventionPatchSeed", "InterventionTimeType",
      "InterventionTimeSeed", "InterventionDispersal",
      "InterventionNicheDistance",
      # Ease of Use
      "SpeciesPreferences", "Intervention"
    )) |> tidytable::left_join(
      InterventionTimes |> tidytable::select(
        TimeIntervention, PoolPatch:PatchAffinitySeed
      ),
      by = c("PoolPatch", "PoolPatchSeed", "Interactions",
             "InteractionsSeed", "Events",
             "EventsSeed", "InitialConditions",
             "InitialConditionsSeed", "Dispersal",
             "NicheDistance", "SpeciesAffinity",
             "SpeciesAffinitySeed", "PatchAffinity",
             "PatchAffinitySeed")
    ) |> tidytable::mutate(
      TimeSinceIntervention = Time - TimeIntervention
    ) |> tidytable::filter(
      SpeciesPreferences %in% figureNetworks$graph$pref,
      NicheDistance == defaultNicheDistance,
      Intervention %in% figureNetworks$graph$interventions | # for KDEs
        (PoolPatchSeed %in% figureNetworks$graph$seed &      # for examples
           Intervention %in% figureNetworks$interventions),  # for examples
      Time == max(figureNetworks$graph$time) |               # for KDEs
        PoolPatchSeed %in% figureNetworks$graph$seed         # for examples
    ) |> tidytable::distinct(
    )

  figureNetworks$graph$timeInterventions <-
    figureNetworks$graph$specification |> tidytable::filter(
      PoolPatchSeed %in% figureNetworks$graph$seed
    ) |> with(
      TimeSinceIntervention[
        outer(
          TimeSinceIntervention,
          figureNetworks$graph$timeInterventions,
          function(x, y) abs(x - y)
        ) |> apply(2, which.min)
        ]
    )

  figureNetworks$graph$specification <-
    figureNetworks$graph$specification |> tidytable::filter(
      round(Time, 6) %in% round(figureNetworks$graph$time, 6) | # for examples
        (round(TimeSinceIntervention, 6) %in%
           round(figureNetworks$graph$timeInterventions, 6))
    )

  figureNetworks$dataRich <- tidytable::bind_rows(
    diversitiesRichness |> tidytable::filter(
      NicheDistance == defaultNicheDistance,
      Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
      PoolPatchSeed %in% basePoolPatchSeeds,
      Metric == "Alpha Hill:0"
    ),
    diversitiesRichness |> tidytable::filter(
      SpeciesPreferences %in% figureNetworks$graph$pref,
      NicheDistance == defaultNicheDistance,
      Intervention %in% figureNetworks$graph$interventions,
      PoolPatchSeed == figureNetworks$graph$seed,
      Metric == "Alpha Hill:0"
    )
  ) |> tidytable::distinct()

  figureNetworks$dataAbund <- diversitiesAbund |> tidytable::filter(
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance"
  )

  figureNetworks$dataTurnover <- diversitiesTimeBC |> tidytable::filter(
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "TimeBrayCurtis: 10"
  )

  # Why to the level of summary? Because the PlotMeanAndInner function
  # isn't built to handle the multiple resolutions that we have in the
  # actual data, which makes it harder to portray the data accurately.
  figureNetworks$dataSummary <- tidytable::bind_rows(
    diversitiesRichness,
    if (6 %in% figures) {# Only case where we need time-progression
      tidytable::bind_rows(
        diversitiesAbund,
        diversitiesTimeBC
      )
    }
  ) |> tidytable::filter(
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0)->(0.5)", # (0)->(0.5) used for N6
                        "(0.5)", "(0.5)->(0)", # all others used in N3.
                        "(0.5)->(1)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric %in% c("Alpha Hill:0", "Alpha Abundance", "TimeBrayCurtis: 10"),
    is.na(Subset) # Overall values
  ) |> tidytable::left_join(
    InterventionTimes
  ) |> tidytable::mutate(
    Metric = factor(
      Metric,
      levels = c("Alpha Hill:0", "Alpha Abundance", "TimeBrayCurtis: 10"),
      labels = c("Richness", "Abundance", "Bray-Curtis Dissimilarity"),
      ordered = TRUE
    ),
    Time = round(Time - TimeIntervention, 6) # remove false differences
  ) |> tidytable::filter(
    Time < 16050, # Need the start for the inset.
    # Avoid singletons.
    abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0,
    Time > -1000
  ) |> tidytable::mutate(
    Time = tidytable::case_when( # Create groupings for times.
      Time < -5 ~ round(Time, -1),
      # In between sampling regimes vs Initial Conditions
      Time < 0 & (InterventionInitial == InterventionFinal) ~ -5,
      Time < 0 & (InterventionInitial != InterventionFinal) ~ -0.1,
      Time <= 50 ~ round(Time, 0), # After 50, jumps to 10s
      Time < 1115 ~ round(Time, -1), # After 1115, slowly expands until ~1863
      Time < 1117 ~ 1110, # Marginal
      Time < 1298 ~ (round(Time/11)*11), # ~1300
      Time < 1395 ~ (round(Time/15)*15), # ~1400
      Time < 1500 ~ (round(Time/30)*30),
      Time < 1600 ~ (round(Time/50)*50),
      Time < 1831 ~ (round(Time/60)*60),
      Time < 16350 ~ round(Time, -2), # 100 gap hols until the end.
      TRUE ~ Time
    )
  ) |> tidytable::group_by(
    # Average Over the now grouped times to make each sim equally weighted.
    # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
    # simple and match each other with the seeds. More complicated set-ups
    # will want to adjust the groupings here.
    Intervention, InterventionInitial, InterventionFinal, Metric,
    PoolPatchSeed, SpeciesPreferences, Time
  ) |> tidytable::summarise(
    Value = median(Value), .groups = "drop"
  ) |> tidytable::group_by(
    # Average across simulations
    Intervention, InterventionInitial, InterventionFinal, Metric,
    SpeciesPreferences, Time
  ) |> tidytable::summarise(
    Lower = quantile(
      Value,
      p = (1 - figureNetworks$ci) - (1 - figureNetworks$ci)/2,
      na.rm = TRUE
    ),
    Average = mean(Value),
    Upper = quantile(
      Value,
      p = figureNetworks$ci + (1 - figureNetworks$ci)/2,
      na.rm = TRUE
    )
  )

  ##### Functions: ############################################################
  figureNetworks$makeViolins <- function(dat) {
    ggplot2::ggplot(
      dat |> tidytable::filter(
        Time > Start, Time < Stop
      ) |> tidytable::group_by(
        PoolPatchSeed, Intervention, SpeciesPreferences, Subset, Metric
      ) |> tidytable::summarise(
        Value = mean(Value)
      ),
      ggplot2::aes(
        x = Intervention,
        y = Value,
        color = Intervention,
        group = paste(Intervention, SpeciesPreferences, Subset)
      )
      # OVERALL Violins
    ) + ggplot2::geom_violin(
      data = function(x) x |> tidytable::filter(is.na(Subset)),
      position = ggplot2::position_dodge(0.9), scale = "count"
    ) + ggplot2::geom_boxplot(
      data = function(x) x |> tidytable::filter(is.na(Subset)),
      notch = TRUE, outlier.size = 0,
      position = ggplot2::position_dodge(0.9),
      width = 0.13
    ) + ggplot2::scale_color_manual(
      values = colorPalette, aesthetics = c("color", "fill"),
      name = "Habitat Type"
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      plot.tag.position = c(0.01, 1),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(1, "lines")
    ) + ggplot2::guides(
      color = "none",
      fill = "none"
    )
  }

  figureNetworks$renameSpeciesPreferences <- function(dat) {
    dat |> tidytable::mutate(
      SpeciesPreferences = tidytable::case_when(
        SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
        SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
        SpeciesPreferences == "Uniform(0, 1)" ~ "Many Adaptation Types",
        TRUE ~ SpeciesPreferences
      )
    )
  }

  #### Cluster of Single Adaptation Type Figures: #############################
  if (2 %in% figures) {
    ##### Figure 2: Richness, Networks through Time ###########################
    ###### Main Plot: #########################################################
    # This is essentially figure 2, but with some more reference points, which
    # we will anchor additional network insets to.
    figureNetworks$plot2A <- plotMeanAndInner(
      figureNetworks$dataRich |> tidytable::filter(
        Intervention %in% c("(0)", "(0.5)", "(1)"),
        SpeciesPreferences == "100% 0",
        is.na(Subset)
      ) |> figureNetworks$renameSpeciesPreferences(
      ),
      CIs = figureNetworks$ci
    ) + ggplot2::geom_point(
      data = function(x) {x |> tidytable::filter(
        PoolPatchSeed == figureNetworks$graph$seed,
        Intervention %in% c("(0)", "(0.5)", "(1)"),
        Time %in% figureNetworks$graph$time
      )},
      mapping = ggplot2::aes(fill = Intervention),
      shape = 21,
      color = "black"
    ) + ggplot2::labs(
      y = "Richness"
    ) + ggplot2::guides(
      linetype = "none",
      color = ggplot2::guide_legend(ncol = 5, title = "Habitat Type"),
      fill = ggplot2::guide_legend(ncol = 5, title = "Habitat Type")
    ) + ggplot2::theme(
      legend.position = c(0.8, 0.09),
      plot.tag.position = c(0.025, 0.95),
      panel.spacing = ggplot2::unit(1, "lines"),
      strip.text = ggplot2::element_text(size = 12)
    ) + ggplot2::coord_cartesian(
      xlim = c(0, 31000), ylim = c(0, richnessYMax), expand = FALSE
    ) + ggplot2::scale_x_continuous(
      breaks = c(0, 1, 10, 100, 1000, (0:3)*10000),
      transform = scales::transform_pseudo_log(sigma = 10)
    ) + ggplot2::facet_grid(
      cols = ggplot2::vars(SpeciesPreferences)
    )

    ###### Inset Plots: #######################################################
    # Network plots at various times, and fuel for the KDEs.
    figureNetworks$graph$networks <- generateNetworks(
      figureNetworks$graph$specification,
      Date = "2025-07-30", split = TRUE
    ); gc() # Tend to have lots of leftover memory usage.

    figureNetworks$indices <-
      figureNetworks$graph$networks$Index |> tidytable::filter(
        # SpeciesPreferences == figureNetworks$pref,
        NicheDistance == defaultNicheDistance,
        Intervention %in% figureNetworks$interventions,
        PoolPatchSeed %in% basePoolPatchSeeds
      ) |> tidytable::arrange(
        ID
      )

    figureNetworks$plot2B <- figureNetworks$indices |> tidytable::filter(
      Time %in% figureNetworks$graph$time,
      PoolPatchSeed %in% figureNetworks$graph$seed,
      SpeciesPreferences == "100% 0",
      Intervention %in% c("(0)", "(0.5)", "(1)")
    ) |> tidytable::pull(ID) |> lapply(
      function(id) {
        list(
          ID = figureNetworks$indices[id, ],
          plt =
            figureNetworks$graph$networks$Envs[[id]]$singletonGraphs[[1]] +
            ggplot2::theme(
              axis.title.x = ggplot2::element_blank(),
              panel.border = ggplot2::element_rect(
                color = "black", fill = NA
              ),
              legend.position = "none",
              axis.text.y = ggplot2::element_text(
                margin = ggplot2::margin(r = -23),
                size = 9, vjust = 0.2),
              axis.text.x = ggplot2::element_blank()
            ) + ggplot2::coord_cartesian(
              xlim = c(-0.5, 1), ylim = c(0.01, 3)
            ) + ggplot2::ylab(
              ""
            ) + ggplot2::scale_size(
              range = c(0.1, 2)
            )
        )
      }
    )

    ###### Summary Plot: ######################################################
    # KDE plots where we are looking at
    figureNetworks$kdes <- lapply(
      figureNetworks$graph$networks$Envs, function(e) {
        e$trophics$EdgeVertexLists[[1]][[1]]$Vertices |> tidytable::select(
          node, Type, Size, N
        ) |> cbind(
          e$Row |> tidytable::select(Time, PoolPatch:Intervention)
        ) |> tidytable::mutate(
          AffinityVals = e$result$Ellipsis$Affinity$SpeciesAffinities[
            as.numeric(substring(node, 2))
            ]
        )
      }
    ) |> tidytable::bind_rows(
    )

    figureNetworks$plot2C <- ggplot2::ggplot(
      figureNetworks$kdes |> tidytable::filter(
        Time == max(figureNetworks$graph$time),
        Intervention %in% c("(0)", "(0.5)", "(1)"),
        SpeciesPreferences == "100% 0"
      )
    ) + ggplot2::geom_density(
      mapping = ggplot2::aes(
        y = Size
      ),
      trim = TRUE, fill = "grey"
    ) + ggplot2::facet_grid(
      factor(Intervention, # Set the order precisely
             levels = c("(1)", "(0.5)", "(0)"), ordered = TRUE) ~ .
    ) + ggplot2::scale_y_log10(
    ) + ggplot2::theme_minimal(
    ) + ggplot2::coord_cartesian(
      xlim = c(0, 4)
    )

    ###### Save: ##############################################################
    ggplot2::ggsave(
      # Use Patchwork to Combine
      figureNetworks$plot2A + figureNetworks$plot2C + patchwork::plot_layout(
        ncol = 2, widths = c(18, 7)
      ),
      path = dirImages,
      filename = "FigureN2_NoIntervention_Base.png",
      units = "cm", width = 25, height = 11
    )

    ggplot2::ggsave(
      # Use Patchwork to Combine
      figureNetworks$plot2A + (
        figureNetworks$plot2C + ggplot2::geom_density(
          mapping = ggplot2::aes(
            y = Size, fill = Type#, color = Intervention
          ),
          alpha = 0.25,
          trim = TRUE, show.legend = FALSE
        ) + ggplot2::scale_fill_manual(
          values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
        )
      ) + patchwork::plot_layout(
        ncol = 2, widths = c(18, 7)
      ),
      path = dirImages,
      filename = "FigureN2_NoIntervention_Layered.png",
      units = "cm", width = 25, height = 11
    )

    # Save insets separately in order to animate them on the presentation.
    lapply(figureNetworks$plot2B, function(lst) {
      ggplot2::ggsave(
        lst$plt,
        path = dirImages,
        filename = paste0(
          "FigureN2_", lst$ID$ID, "_",
          # anything in (, ), ., -, or > needs to be eliminated for filename.
          gsub(lst$ID$Intervention, pattern = "[()(.)>-]", replacement = ""),
          "_", lst$ID$Time, ".png"
        ),
        units = "cm", width = 4, height = 3
      )
    })
  }

  ##### Figure 3: Intervention, Richness, Networks through Time ###############
  # Figure 3 looks a lot like figure 2, but we're going to focus on the
  # post-intervention period, and the comparison with Figure 2 instead of on
  # the KDEs.
  if (3 %in% figures) {
    ###### Main Plot: #########################################################
    figureNetworks$plot3A <- ggplot2::ggplot(
      figureNetworks$dataSummary |> tidytable::filter(
        Metric == "Richness", Time >= -100,
        Intervention %in% c("(0.5)->(0)", "(0.5)", "(0.5)->(1)"),
        SpeciesPreferences == "100% 0"
      ) |> figureNetworks$renameSpeciesPreferences(
      ),
      aes(x = Time, y = Average,
          color = Intervention,
          fill = Intervention
      )
    ) + ggplot2::geom_vline(
      xintercept = 0, color = "black", linetype = "dashed"
    ) + ggplot2::geom_line(
    ) + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      alpha = 0.25, linewidth = 0.25
    ) + ggplot2::geom_point(
      data = figureNetworks$dataRich |> tidytable::filter(
        # Two step filter to reduce computation as much as possible.
        PoolPatchSeed == figureNetworks$graph$seed,
        Intervention %in% c("(0.5)->(0)", "(0.5)", "(0.5)->(1)"),
        is.na(Subset), SpeciesPreferences == "100% 0"
      ) |> tidytable::left_join(
        InterventionTimes
      ) |> tidytable::mutate(
        Time = Time - TimeIntervention
      ) |> tidytable::filter(
        Time %in% figureNetworks$graph$timeInterventions |
          round(Time, 6) %in% figureNetworks$graph$timeInterventions,
        Time < 2000 # Don't want the 10000 here.
      ) |> figureNetworks$renameSpeciesPreferences(
      ),
      mapping = ggplot2::aes(fill = Intervention, y = Value),
      shape = 21,
      color = "black"
    ) + ggplot2::scale_color_manual(
      values = colorPalette, aesthetics = c("color", "fill"),
      name = ""
    ) + ggplot2::guides(
      linetype = "none",
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      legend.position = c(0.58, 0.88),
      panel.spacing = ggplot2::unit(1, "lines"),
      strip.text = ggplot2::element_text(size = 12)
    ) + ggplot2::labs(
      x = "Time Since Intervention",
      y = "Richness"
    ) + ggplot2::coord_cartesian(
      ylim = c(0, richnessYMax), xlim = c(-20, NA),
      expand = FALSE
    ) + ggplot2::scale_x_continuous(
      breaks = c(0, 1, 10, 100, 1000, 10000, 15000),
      transform = scales::transform_pseudo_log(sigma = 10)
    ) + ggplot2::facet_grid(
      cols = ggplot2::vars(SpeciesPreferences)
    )

    ###### Inset Plots: #######################################################
    # Network plots at various times.
    # If we've already had the KDEs from F2, then we don't need to redo all of
    # this. If we don't then we want to only do the bits we'll use. KDEs mean
    # doing pretty much all (so minor addt'l cost), but here we use only some.
    figureNetworks$graph$networksSubset <-
      if("networks" %in% names(figureNetworks$graph)) {
        # NOT TESTED:
        targets <- figureNetworks$graph$networks$Index |> tidytable::filter(
          PoolPatchSeed == figureNetworks$graph$seed,
          Intervention %in% c("(0.5)->(0)", "(0.5)", "(0.5)->(1)"),
          SpeciesPreferences == "100% 0",
          TimeSinceIntervention > -1
        ) |> tidytable::pull(ID)
        list(
          Envs = figureNetworks$graph$networks[targets],
          Index = figureNetworks$graph$networks$Index[
            targets,
            ] |> tidytable::mutate(
              ID = tidytable::row_number()
            )
        )
      } else {
        generateNetworks(
          figureNetworks$graph$specification |> tidytable::filter(
            PoolPatchSeed == figureNetworks$graph$seed,
            Intervention %in% c("(0.5)->(0)", "(0.5)", "(0.5)->(1)"),
            SpeciesPreferences == "100% 0",
            TimeSinceIntervention > -1
          ),
          Date = "2025-07-30", split = TRUE
        ); gc() # Tend to have lots of leftover memory usage.
      }

    figureNetworks$plot3B <-
      figureNetworks$graph$networksSubset$Index |> tidytable::pull(ID) |> lapply(
        function(id) {
          list(
            ID = figureNetworks$graph$networksSubset$Index[id, ],
            plt =
              figureNetworks$graph$networksSubset$Envs[[
                id
                ]]$singletonGraphs[[1]] +
              ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(
                  color = "black", fill = NA
                ),
                legend.position = "none",
                axis.text.y = ggplot2::element_text(
                  margin = ggplot2::margin(r = -23),
                  size = 9, vjust = 0.2),
                axis.text.x = ggplot2::element_blank()
              ) + ggplot2::coord_cartesian(
                xlim = c(-0.5, 1), ylim = c(0.01, 3)
              ) + ggplot2::ylab(
                ""
              ) + ggplot2::scale_size(
                range = c(0.1, 2)
              )
          )
        }
      )

    ###### Comparison Plot: ###################################################
    figureNetworks$plot3C <- ggplot2::ggplot(
      figureNetworks$dataSummary |> tidytable::filter(
        Metric == "Richness", Time >= 9000,
        Intervention %in% c("(0)", "(0.5)", "(1)"),
        SpeciesPreferences == "100% 0"
      ) |> figureNetworks$renameSpeciesPreferences(
      ),
      aes(x = Time, y = Average,
          color = Intervention,
          fill = Intervention
      )
    ) + ggplot2::geom_vline(
      xintercept = 0, color = "black", linetype = "dashed"
    ) + ggplot2::geom_line(
    ) + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      alpha = 0.25, linewidth = 0.25
    ) + ggplot2::scale_color_manual(
      values = colorPalette, aesthetics = c("color", "fill"),
      name = ""
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      legend.position = c(0.5, 0.15),
      panel.spacing = ggplot2::unit(1, "lines"),
      strip.text = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_blank()
    ) + ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) + ggplot2::labs(
      x = "Time Since Intervention"
    ) + ggplot2::coord_cartesian(
      ylim = c(0, richnessYMax), xlim = c(10000, 15000),
      expand = FALSE
    ) + ggplot2::scale_x_continuous(
      breaks = c(10000, 15000),
      transform = scales::transform_pseudo_log(sigma = 10)
    )

    ###### Save: ##############################################################
    ggplot2::ggsave(
      # Use Patchwork to Combine
      figureNetworks$plot3A + figureNetworks$plot3C + patchwork::plot_layout(
        ncol = 2, widths = c(23, 2)
      ),
      path = dirImages,
      filename = "FigureN3_Intervention_Base.png",
      units = "cm", width = 20, height = 11
    )

    # Save insets separately in order to animate them on the presentation.
    lapply(figureNetworks$plot3B, function(lst) {
      ggplot2::ggsave(
        lst$plt,
        path = dirImages,
        filename = paste0(
          "FigureN3_", lst$ID$ID, "_",
          # anything in (, ), ., -, or > needs to be eliminated for filename.
          gsub(lst$ID$Intervention, pattern = "[()(.)>-]", replacement = ""),
          "_", lst$ID$Time, ".png"
        ),
        units = "cm", width = 4, height = 3
      )
    })
  }

  ##### Figure 4: Richness, Abundance, Turnover, Complexity (RATC) ############
  if (4 %in% figures) {
    ###### Prep Connectance Data: #############################################
    # This is buried in the graphs of all of the relevant simulations. We'll
    # take specifically from the t = 25000 CTU data, but note that we need to
    # keep things comparable for RAT as well, so we need all non-int. data.
    #

    intermediate <- generateNetworks( # high memory intermediate to be rm'd.
      # 5 Interventions x 44 simulations = 220
      figureNetworks$dataRich |> tidytable::filter(
        Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
        SpeciesPreferences == "100% 0", is.na(Subset),
        Time == 25000
      ) |> tidytable::select(-Metric, -Value),
      Date = "2025-07-30", split = TRUE
    ); gc() # Tend to have lots of leftover memory usage.

    figureNetworks$dataConnectance <- lapply(
      intermediate$Envs,
      function(env) {
        val <- env$graphs[[1]] %E>% filter(# only 1!
          to != from # Remove self-edges/loops; bias connectance
        ) |> igraph::edge_density(loops = FALSE)
        return(env$Row |> tidytable::mutate(
          Value = val,
          Metric = "Connectance"
        ))
      }
    ) |> tidytable::bind_rows()

    rm(intermediate); gc()

    ###### Violin Plots: ######################################################
    figureNetworks$plot4A <- (
      figureNetworks$dataRich |> tidytable::filter(
        Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
        SpeciesPreferences %in% c("100% 0")
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::coord_cartesian(
      ylim = c(0, richnessYMax), expand = FALSE
    ) + ggplot2::labs(
      y = "Richness",
      x = "Habitat Type"
    )

    figureNetworks$plot4B <- (
      figureNetworks$dataAbund |> tidytable::filter(
        SpeciesPreferences %in% c("100% 0")
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::labs(
      y = "Abundance",
      x = "Habitat Type"
    ) + ggplot2::scale_y_continuous(
      transform = "pseudo_log", breaks = 10^(4+(-4:4)/4),
      label = scales::label_log(digits = 2), limits = 10^c(3.25, 4.75)
    )

    figureNetworks$plot4C <- (
      figureNetworks$dataTurnover |> tidytable::filter(
        SpeciesPreferences %in% c("100% 0")
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::labs(
      y = "Bray-Curtis Dissimilarity",
      x = "Habitat Type"
    )

    figureNetworks$plot4D <- (
      figureNetworks$dataConnectance |> tidytable::filter(
        SpeciesPreferences %in% c("100% 0")
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::labs(
      y = "Connectance",
      x = "Habitat Type"
    )

    ###### Sense-check: #######################################################
    if (require("GGally")) {
      figureNetworks$sensecheck4 <- tidytable::bind_rows(
        figureNetworks$dataRich,
        figureNetworks$dataAbund |> tidytable::mutate(
          Metric = gsub(pattern = "Alpha", replacement = "Log10",
                        x = Metric, fixed = TRUE),
          Value = log10(Value)
        ),
        figureNetworks$dataTurnover
      ) |> tidytable::filter(
        Time == 25000, is.na(Subset), SpeciesPreferences == "100% 0"
      ) |> tidytable::bind_rows(
        figureNetworks$dataConnectance
      ) |> tidytable::pivot_wider(
        names_from = "Metric", values_from = "Value"
      ) |> tidytable::mutate(
        # For inspection, but not otherwise used.
        PerSpeciesAbundance = 10^`Log10 Abundance` / `Alpha Hill:0`,
        PerSpeciesBC = `TimeBrayCurtis: 10` / `Alpha Hill:0`,
        PerIndividualBC =  `TimeBrayCurtis: 10` / (10^`Log10 Abundance`),
        PerNodeEdges = (`Alpha Hill:0` - 1) * Connectance,
        Edges = `Alpha Hill:0` * (`Alpha Hill:0` - 1) * Connectance
      )

      figureNetworks$plot4E <- GGally::ggpairs(
        figureNetworks$sensecheck4,
        columns = c( # Each should be length(unique(...)) == 1.
          figureNetworks$dataRich$Metric[1],
          gsub(pattern = "Alpha", replacement = "Log10",
               x = figureNetworks$dataAbund$Metric[1], fixed = TRUE),
          figureNetworks$dataTurnover$Metric[1],
          figureNetworks$dataConnectance$Metric[1]
        ),
        mapping = ggplot2::aes(
          color = Intervention,
          group = Intervention,
          alpha = 0.25
        )
      ) + ggplot2::scale_color_manual(
        values = colorPalette, aesthetics = c("color", "fill"),
        name = "Habitat Type"
      ) + ggplot2::theme_minimal(
      )
    }

    ###### Save: ##############################################################
    ggplot2::ggsave(
      # Use Patchwork to Combine
      figureNetworks$plot4A + figureNetworks$plot4B +
        figureNetworks$plot4C + figureNetworks$plot4D + patchwork::plot_layout(
          ncol = 2, nrow = 2
        ),
      path = dirImages,
      filename = "FigureN4_Solos.png",
      units = "cm", width = 20, height = 11
    )

    if (require("GGally")) {
      ggplot2::ggsave(
        figureNetworks$plot4E,
        path = dirImages,
        filename = "FigureN4_Combos.png",
        units = "cm", width = 25, height = 14
      )
    }

    ((figureNetworks$sensecheck4 |> tidytable::pivot_longer(
      names_to = "Metric", values_to = "Value", cols = c(
        figureNetworks$dataRich$Metric[1],
        gsub(pattern = "Alpha", replacement = "Log10",
             x = figureNetworks$dataAbund$Metric[1], fixed = TRUE),
        figureNetworks$dataTurnover$Metric[1],
        figureNetworks$dataConnectance$Metric[1],
        "PerSpeciesAbundance",
        "PerSpeciesBC",
        "PerIndividualBC" ,
        "PerNodeEdges",
        "Edges"
      )) |> figureNetworks$makeViolins()) + ggplot2::facet_wrap(
        ggplot2::vars(Metric), scales = "free"
      ) + ggplot2::theme(
        strip.text = ggplot2::element_text()
      )) |> ggplot2::ggsave(
        path = dirImages,
        filename = "FigureN4_Overview.png",
        units = "cm", width = 25, height = 14
      )
  }

  #### 2 and Many Adaptation Type Figures: ################################
  ##### Figure 5: Richness w/Time, Abundance, Turnover, Complexity ############
  if (5 %in% figures) {
    ###### Prep Connectance Data: #############################################
    # This is buried in the graphs of all of the relevant simulations. We'll
    # take specifically from the t = 25000 CTU data, but note that we need to
    # keep things comparable for RAT as well, so we need all non-int. data.
    #

    intermediate <- generateNetworks( # high memory intermediate to be rm'd.
      # 5 Interventions x 44 simulations x 1 Prefs = 220
      figureNetworks$dataRich |> tidytable::filter(
        Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
        SpeciesPreferences == pref, is.na(Subset),
        Time == 25000
      ) |> tidytable::select(-Metric, -Value),
      Date = "2025-07-30", split = TRUE
    ); gc() # Tend to have lots of leftover memory usage.

    figureNetworks$dataConnectance <- lapply(
      intermediate$Envs,
      function(env) {
        val <- env$graphs[[1]] %E>% filter(# only 1!
          to != from # Remove self-edges/loops; bias connectance
        ) |> igraph::edge_density(loops = FALSE)
        return(env$Row |> tidytable::mutate(
          Value = val,
          Metric = "Connectance"
        ))
      }
    ) |> tidytable::bind_rows()

    rm(intermediate); gc()

    ###### Time Series Plot: ##################################################
    figureNetworks$plot5A <- plotMeanAndInner(
      figureNetworks$dataRich |> tidytable::filter(
        Intervention %in% c("(0)", "(0.5)", "(1)"),
        SpeciesPreferences != "100% 0", # Both prefs here.
        is.na(Subset)
      ) |> figureNetworks$renameSpeciesPreferences(
      ),
      CIs = figureNetworks$ci
    ) + ggplot2::labs(
      y = "Richness"
    ) + ggplot2::guides(
      linetype = "none",
      color = ggplot2::guide_legend(ncol = 1, title = "Habitat Type"),
      fill = ggplot2::guide_legend(ncol = 1, title = "Habitat Type")
    ) + ggplot2::theme(
      legend.position = c(0.25, 0.8),
      plot.tag.position = c(0.025, 0.95),
      panel.spacing = ggplot2::unit(1, "lines"),
      strip.text = ggplot2::element_text(size = 12)
    ) + ggplot2::coord_cartesian(
      xlim = c(0, 31000), ylim = c(0, richnessYMax), expand = FALSE
    ) + ggplot2::scale_x_continuous(
      breaks = c(0, 10,
                 100, 1000,
                 (1:3)*10000),
      labels = c("0", "10",
                 expression("10"^"2"), expression("10"^"3"),
                 expression(paste("10"^"4", phantom(100))), "",
                 expression(paste("3"%*%"10"^"5", phantom(0)))),
      transform = scales::transform_pseudo_log(sigma = 10)
    ) + ggplot2::facet_wrap(
      ggplot2::vars(SpeciesPreferences), nrow = 2
    )

    ###### RATC: ##############################################################
    ###### Violin Plots: ######################################################
    figureNetworks$plot5B1 <- (
      figureNetworks$dataRich |> tidytable::filter(
        SpeciesPreferences %in% pref
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::coord_cartesian(
      ylim = c(0, richnessYMax), expand = FALSE
    ) + ggplot2::labs(
      y = "Richness",
      x = "Habitat Type"
    )

    figureNetworks$plot5C1 <- (
      figureNetworks$dataAbund |> tidytable::filter(
        SpeciesPreferences %in% pref
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::labs(
      y = "Abundance",
      x = "Habitat Type"
    ) + ggplot2::scale_y_continuous(
      transform = "pseudo_log", breaks = 10^(4+(-4:4)/4),
      label = scales::label_log(digits = 2), limits = 10^c(3.25, 4.75)
    )

    figureNetworks$plot5D1 <- (
      figureNetworks$dataTurnover |> tidytable::filter(
        SpeciesPreferences %in% pref
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::coord_cartesian(
      ylim = c(0, 0.125), expand = FALSE
    ) + ggplot2::labs(
      y = "Bray-Curtis Dissimilarity",
      x = "Habitat Type"
    )

    figureNetworks$plot5E1 <- (
      figureNetworks$dataConnectance |> tidytable::filter(
        SpeciesPreferences %in% pref
      ) |> figureNetworks$makeViolins()
    ) + ggplot2::coord_cartesian(
      ylim = c(0, 0.5), expand = TRUE
    ) + ggplot2::labs(
      y = "Connectance",
      x = "Habitat Type"
    )

    ###### Sense-check: #######################################################
    if (require("GGally")) {
      figureNetworks$sensecheck5 <- tidytable::bind_rows(
        figureNetworks$dataRich,
        figureNetworks$dataAbund |> tidytable::mutate(
          Metric = gsub(pattern = "Alpha", replacement = "Log10",
                        x = Metric, fixed = TRUE),
          Value = log10(Value)
        ),
        figureNetworks$dataTurnover
      ) |> tidytable::filter(
        Time == 25000, is.na(Subset), SpeciesPreferences == pref
      ) |> tidytable::bind_rows(
        figureNetworks$dataConnectance
      ) |> tidytable::pivot_wider(
        names_from = "Metric", values_from = "Value"
      ) |> tidytable::mutate(
        # For inspection, but not otherwise used.
        PerSpeciesAbundance = 10^`Log10 Abundance` / `Alpha Hill:0`,
        PerSpeciesBC = `TimeBrayCurtis: 10` / `Alpha Hill:0`,
        PerIndividualBC =  `TimeBrayCurtis: 10` / (10^`Log10 Abundance`),
        PerNodeEdges = (`Alpha Hill:0` - 1) * Connectance,
        Edges = `Alpha Hill:0` * (`Alpha Hill:0` - 1) * Connectance
      )

      figureNetworks$plot5F1 <- GGally::ggpairs(
        figureNetworks$sensecheck5 |> tidytable::filter(
          SpeciesPreferences == pref
        ),
        columns = c( # Each should be length(unique(...)) == 1.
          figureNetworks$dataRich$Metric[1],
          gsub(pattern = "Alpha", replacement = "Log10",
               x = figureNetworks$dataAbund$Metric[1], fixed = TRUE),
          figureNetworks$dataTurnover$Metric[1],
          figureNetworks$dataConnectance$Metric[1]
        ),
        mapping = ggplot2::aes(
          color = Intervention,
          group = Intervention,
          alpha = 0.25
        )
      ) + ggplot2::scale_color_manual(
        values = colorPalette, aesthetics = c("color", "fill"),
        name = "Habitat Type"
      ) + ggplot2::theme_minimal(
      ) + ggplot2::ggtitle(
        pref
      )
    }

    ###### Save: ##############################################################
    ggplot2::ggsave(
      # Use Patchwork to Combine
      figureNetworks$plot5A +
        figureNetworks$plot5B1 + figureNetworks$plot5C1 +
        figureNetworks$plot5D1 + figureNetworks$plot5E1 +
        patchwork::plot_layout(
          design = "AABBCC
                    AABBCC
                    AADDEE
                    AADDEE"
        ),
      path = dirImages,
      filename = paste0("FigureN5_", prefstring, ".png"),
      units = "cm", width = 20, height = 11
    )

    if (require("GGally")) {
      ggplot2::ggsave(
        figureNetworks$plot5F1,
        path = dirImages,
        filename = paste0("FigureN5_Combos_", prefstring, ".png"),
        units = "cm", width = 25, height = 14
      )
    }

    ((figureNetworks$sensecheck5 |> tidytable::filter(
      SpeciesPreferences == pref
    ) |> tidytable::pivot_longer(
      names_to = "Metric", values_to = "Value", cols = c(
        figureNetworks$dataRich$Metric[1],
        gsub(pattern = "Alpha", replacement = "Log10",
             x = figureNetworks$dataAbund$Metric[1], fixed = TRUE),
        figureNetworks$dataTurnover$Metric[1],
        figureNetworks$dataConnectance$Metric[1],
        "PerSpeciesAbundance",
        "PerSpeciesBC",
        "PerIndividualBC" ,
        "PerNodeEdges",
        "Edges"
      )) |> figureNetworks$makeViolins()) + ggplot2::facet_wrap(
        ggplot2::vars(Metric), scales = "free"
      ) + ggplot2::theme(
        strip.text = ggplot2::element_text()
      ) + ggplot2::ggtitle(
        pref
      )) |> ggplot2::ggsave(
        path = dirImages,
        filename = paste0("FigureN5_Overview_", prefstring, ".png"),
        units = "cm", width = 25, height = 14
      )
  }
  ##### Figure 6: Short Term Int. RAT #########################################
  if (6 %in% figures) {
    # Two variants:
    #   Time series of the RAT (no C, would need to massively refactor code).
    #     Unfortunately, due to holes in my sampling strategy, BrayCurtis
    #     doesn't have all the points necessary to look convincing.
    #       (1, 2, & 3 end early, 10 has holes, and these are substitutable.)
    #   Progression of the networks with R only.
    #
    ###### Prep KDE Data: #####################################################
    # Goal: x = "Affinity", y = "Size", scatter plot at a time, with marginal
    # KDE plots, so that we can show how the distribution of species shifts
    # after the intervention for a system with uniform habitat adaptations.

    ColExtInterventionStrings <- ColExt %>% dplyr::select(
      PatchAffinity, PoolPatch, InterventionPatchType
    ) %>% dplyr::distinct(
    ) %>% dplyr::mutate(
      Intervention = unlist(mapply(
        FUN = interventionNamingScheme,
        PatchAffinity, PoolPatch, InterventionPatchType
      ))
    )

    figureNetworks$ColExt <- ColExt |> tidytable::filter(
      EventType != "Present", Success, # Only Events that Did Something.
      NicheDistance == defaultNicheDistance # Base case.
    ) |> tidytable::mutate( # Make human readable formats.
      SpeciesPreferences =
        speciesAffinityDictionaryOrigin$SpeciesAffinities[
          as.numeric(SpeciesAffinity)
          ]
    ) |> tidytable::left_join(
      ColExtInterventionStrings,
      by = c("PatchAffinity", "PoolPatch", "InterventionPatchType"),
      multiple = "all"
    ) |> changePreferencesLevels(
    ) |> changeInterventionLevels(
    ) |> tidytable::filter( # Reduce to the levels we're interested in.
      SpeciesPreferences == "Uniform(0, 1)",
      Intervention %in% c("(0)->(0.5)", "(0.5)->(0)",  # To display.
                          "(0)", "(0.5)") # To sense check.
    ) |> tidytable::left_join(
      endTimes |> dplyr::select(-Times)
    ) |> tidytable::left_join( # So we can calculate our post-int. slices
      InterventionTimes
    ) |> tidytable::group_by(
      # Species Characteristics
      Species, Environment, SpeciesType, Size, ReproductionRate, Speed,
      Affinity,
      # Simulation Characteristics
      PoolPatch:TimeIntervention
    ) |> tidytable::mutate( # Create the groupings for when present in system.
      In = ifelse(EventType == "Arrival", "In", "Out"),
      InOutNumber = cumsum(In == "In"),
      Time = Times - TimeIntervention
    ) |> tidytable::pivot_wider( # Rearrange so we have time ranges
      names_from = In, values_from  = Time,
      id_cols = c(
        # Species Characteristics
        Species, Environment, SpeciesType, Size, ReproductionRate, Speed,
        Affinity,
        # Simulation Characteristics
        PoolPatch:TimeIntervention,
        InOutNumber
      )
    ) |> tidytable::filter( # Remove ranges that aren't included in any slice.
      outer( # Check all combinations intervention and range, then summarize.
        FUN = function(X, Y, low, high) {low[X] < Y & Y < high[X]},
        X = seq_along(In),
        Y = figureNetworks$graph$timeInterventions,
        low = In, high = Out
      ) |> apply(MARGIN = 1, FUN = any)
    )


    ###### Main Plot: #########################################################
    figureNetworks$plot6A <- lapply(
      unique(figureNetworks$dataSummary$Metric), function(met, dat) {
        plt <- ggplot2::ggplot(
          dat |> tidytable::filter(
            Time >= -100, Metric == met,
            Intervention %in% c("(0)", "(0.5)->(0)", "(0.5)", "(0)->(0.5)"),
            SpeciesPreferences == "Uniform(0, 1)"
          ) |> figureNetworks$renameSpeciesPreferences(
          ),
          aes(
            x = Time, y = Average,
            color = Intervention,
            fill = Intervention,
            group = interaction(Intervention, Metric),
            alpha = Intervention
          )
        ) + ggplot2::geom_vline(
          xintercept = 0, color = "black", linetype = "dashed"
        ) + ggplot2::geom_line(
        ) + ggplot2::scale_color_manual(
          values = colorPalette, aesthetics = c("color", "fill"),
          name = ""
        ) + ggplot2::guides(
          linetype = "none",
          fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
        ) + ggplot2::theme_minimal(
        ) + ggplot2::theme(
          legend.position = if(met == "Richness") {c(0.05, 0.35)} else {"none"},
          panel.spacing = ggplot2::unit(1, "lines"),
          strip.text = ggplot2::element_text(size = 12)
        ) + ggplot2::labs(
          x = "Time Since Intervention",
          y = met
        ) + ggplot2::coord_cartesian(
          xlim = c(-20, NA)
        ) + ggplot2::scale_x_continuous(
          breaks = c(0, 1, 10, 100, 1000, 10000, 15000),
          transform = scales::transform_pseudo_log(sigma = 10)
          # ) + ggplot2::facet_grid(
          #   cols = ggplot2::vars(SpeciesPreferences)
        ) + ggplot2::scale_alpha_manual(
          values = c("(0)" = 1, "(0.5)->(0)" = 1,
                     "(0.5)" = 1, "(0)->(0.5)" = 1),
          guide = "none"
        )
        
        return(plt)
      },
      dat = figureNetworks$dataSummary
    )

    ###### Network Plots: ####################################################
    figureNetworks$plot6B <- with(
      list(
        plotfun = function(inter, time) {
          list(
            (
              figureNetworks$ColExt |> tidytable::filter(
                In < time, time < Out, Intervention == inter
              ) |> ggplot2::ggplot(
                ggplot2::aes(
                  x = Affinity, y = Size, color = SpeciesType
                )
              ) + ggplot2::geom_point(
                show.legend = FALSE
              ) + ggplot2::scale_y_log10(
              ) + ggplot2::scale_x_continuous(
                breaks = c(0, 0.5, 1)
              ) + ggplot2::theme_minimal(
              ) + ggplot2::labs(
                x = if(time == 0) "Habitat Adaptation",
                y = if(time == 0) "Species Size"
              ) + scale_color_manual(
                values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
              ) + ggplot2::coord_cartesian(
                xlim = c(0, 1), ylim = c(10^-2, 10^0.5), expand = FALSE
              )
            ) |> ggExtra::ggMarginal(
            )
          )
        }
      ),
      expand_grid(
        Intervention = sort(unique(figureNetworks$ColExt$Intervention)),
        Time = figureNetworks$graph$timeInterventions
      ) |> tidytable::mutate_rowwise(
        Plot = plotfun(Intervention, Time)
      )
    )

    ###### Save: #############################################################
    for (i in 1+5*(0:3)) {
      ggplot2::ggsave(
        # Use Patchwork to Combine
        with(figureNetworks,
             plot6A[[1]] + ggplot2::geom_vline(
               xintercept = plot6B$Time[i+0:4], linetype = "dashed"
             )  + ggplot2::scale_alpha_manual(
               values = (table(plot6B$Intervention[i+0:4])+1)/6,
               guide = "none"
             ) +
               plot6B$Plot[[i]] + plot6B$Plot[[i+1]] + plot6B$Plot[[i+2]] +
               plot6B$Plot[[i+3]] + plot6B$Plot[[i+4]] +
               patchwork::plot_layout(
                 design =
                   "BBCCDDEEFF
                    BBCCDDEEFF
                    AAAAAAAAAA
                    AAAAAAAAAA"
               )
        ),
        path = dirImages,
        filename = paste0("FigureN6_RichnessAndKDEs_", i, ".png"),
        units = "cm", width = 20, height = 11
      )
    }

    ggplot2::ggsave(
      # Use Patchwork to Combine
      with(figureNetworks,
           plot6A[[1]] + ggplot2::geom_vline(
             xintercept = plot6B$Time[i+0:4], linetype = "dashed"
           ) +
             plot6B$Plot[[i]] + plot6B$Plot[[i+1]] + plot6B$Plot[[i+2]] +
             plot6B$Plot[[i+3]] + plot6B$Plot[[i+4]] +
             patchwork::plot_layout(
               design =
                 "BBCCDDEEFF
                  BBCCDDEEFF
                  AAAAAAAAAA
                  AAAAAAAAAA"
             )
      ),
      path = dirImages,
      filename = paste0("FigureN6_RichnessAndKDEs_", i, "_noalpha.png"),
      units = "cm", width = 20, height = 11
    )

    if (require("magick")) {
      magick::image_animate(
        c(
          magick::image_read(
            file.path(dirImages, "FigureN6_RichnessAndKDEs_16.png")
          ),
          magick::image_read(
            file.path(dirImages, "FigureN6_RichnessAndKDEs_11.png")
          )
        ),
        fps = 1, optimize = TRUE
      ) |> magick::image_write(
        file.path(dirImages, "FigureN6_RichnessAndKDEs_1611.gif")
      )
    }
  }

  #### Summary Images: ########################################################
  ##### Figure 7: Parameters Cause RATC #######################################
  if (7 %in% figures) {
    ((diversitiesRichness |> tidytable::filter(
      SpeciesPreferences %in% c("100% 0", "Uniform(0, 1)"),
      NicheDistance == defaultNicheDistance,
      Metric == "Alpha Hill:0", is.na(Subset),
      Start < Time, Time < Stop
    ) |> tidytable::group_by(
      Metric, Environment1, Environment2, PoolPatch:Stop
    ) |> tidytable::summarise(
      Average = mean(Value), # Average Richness within Simulation
      .groups = "drop"
    ) |> tidytable::group_by(
      Metric, PoolPatch, Interactions, Events, InitialConditions, Dispersal,
      NicheDistance, SpeciesAffinity, PatchAffinity, InterventionPatchType,
      InterventionTimeType, InterventionDispersal, InterventionNicheDistance,
      Intervention, SpeciesPreferences, InterventionInitial, InterventionFinal,
      DispersalParam
    ) |> tidytable::summarise(
      Average = mean(Average), # Average Richness across Simulations
      .groups = "drop" # (Simulations evenly weighted)
    ) |> tidytable::mutate(
      SpeciesPreferences = tidytable::case_when(
        SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
        SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
        SpeciesPreferences == "Uniform(0, 1)" ~ "Many Adaptation Types",
        TRUE ~ SpeciesPreferences
      )
    )) |> ggplot2::ggplot(
      ggplot2::aes(
        x = InterventionInitial,
        y = InterventionFinal,
        fill = Average
      )
    ) + ggplot2::geom_tile(
      width = 1, height = 1, color = NA
    ) + ggplot2::geom_text(
      ggplot2::aes(label = round(Average))
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      strip.text = ggplot2::element_text(size = 10.5)
    ) + ggplot2::labs(
      fill = "Avg.\nRichness",
      x = "Initial Habitat Type",
      y = "Final Habitat Type"
    ) + ggplot2::coord_fixed(
    ) + ggplot2::facet_wrap(
      ncol = 1, ggplot2::vars(SpeciesPreferences)
    ) + ggplot2::scale_fill_viridis_c(
      begin = 0.25
    )) |> ggplot2::ggsave(
      filename = "FigureN7_TileSummary.png",
      path = dirImages,
      units = "cm", width = 11, height = 11
    )
  }

  ##### Figure 8: Network reorganisation over short time scales ###############
  if (8 %in% figures) {
    ((diversitiesRichness |> tidytable::filter(
      SpeciesPreferences %in% c("100% 0", "Uniform(0, 1)"),
      Intervention %in% c("(0)", "(0)->(0.5)", "(0)->(1)",
                          "(0.5)", "(0.5)->(0)", "(0.5)->(1)",
                          "(1)", "(1)->(0.5)", "(1)->(0)"),
      NicheDistance == defaultNicheDistance,
      Metric == "Alpha Hill:0", is.na(Subset)
    ) |> tidytable::left_join(
      InterventionTimes
    ) |> tidytable::mutate(
      Time = Time - TimeIntervention
    ) |> tidytable::filter(
      abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0,
      Time > -25, Time <= 1000
    ) |> tidytable::mutate(
      # Follow what was done for the above figures.
      Time = tidytable::case_when( # Create groupings for times.
        Time < -5 ~ round(Time, -1),
        # In between sampling regimes vs Initial Conditions
        Time < 0 & (InterventionInitial == InterventionFinal) ~ -5,
        Time < 0 & (InterventionInitial != InterventionFinal) ~ -0.1,
        Time <= 50 ~ round(Time, 0), # After 50, jumps to 10s
        Time < 1115 ~ round(Time, -1), # After 1115, slowly expands until ~1863
        Time < 1117 ~ 1110, # Marginal
        Time < 1298 ~ (round(Time/11)*11), # ~1300
        Time < 1395 ~ (round(Time/15)*15), # ~1400
        Time < 1500 ~ (round(Time/30)*30),
        Time < 1600 ~ (round(Time/50)*50),
        Time < 1831 ~ (round(Time/60)*60),
        Time < 16350 ~ round(Time, -2), # 100 gap hols until the end.
        TRUE ~ Time
      )
    ) |> tidytable::group_by(
      # Average Over the now grouped times to make each sim equally weighted.
      # NOTE TO FUTURE USERS, I've been lazy here because the groupings are
      # simple and match each other with the seeds. More complicated set-ups
      # will want to adjust the groupings here.
      Intervention, InterventionInitial, InterventionFinal, Metric,
      PoolPatchSeed, SpeciesPreferences, Time
    ) |> tidytable::summarise(
      Value = median(Value), .groups = "drop"
    ) |> tidytable::group_by(
      # Average across simulations
      Intervention, InterventionInitial, InterventionFinal, Metric,
      SpeciesPreferences, Time
    ) |> tidytable::summarise(
      Average = mean(Value),.groups = "drop"
    ) |> tidytable::mutate(
      SpeciesPreferences = tidytable::case_when(
        SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
        SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
        SpeciesPreferences == "Uniform(0, 1)" ~ "Many Adaptation Types",
        TRUE ~ SpeciesPreferences
      ),
      Linetype = ifelse(InterventionInitial != InterventionFinal, "dotted", "solid")
    )) |> ggplot2::ggplot(
      ggplot2::aes(
        x = Time,
        y = Average,
        color = Intervention,
        linetype = Linetype
      )
    ) + ggplot2::geom_vline(
      xintercept = 0, color = "black", linetype = "dashed"
    ) + ggplot2::geom_line(
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      strip.text = ggplot2::element_text(size = 10.5),
      axis.text.x = ggplot2::element_text(hjust = 1)
    ) + ggplot2::scale_color_manual(
      values = colorPalette, aesthetics = c("color", "fill"),
      name = "", guide = "none"
    ) + ggplot2::labs(
      y = "Avg.\nRichness",
      x = "Time Since Intervention"
    ) + ggplot2::facet_wrap(
      ncol = 1, ggplot2::vars(SpeciesPreferences), scales = "free_y"
    ) + ggplot2::coord_cartesian(
      xlim = c(-20, NA),
      expand = FALSE
    ) + ggplot2::scale_x_continuous(
      breaks = c(0, 1, 10, 100, 1000, 10000, 15000),
      transform = scales::transform_pseudo_log(sigma = 10)
    ) + ggplot2::scale_linetype_manual(
      values = c("dashed" = "dashed", "solid" = "solid", "dotted" = "dotted"),
      guide = "none"
    )) |> ggplot2::ggsave(
      filename = "FigureN8_LineSummary.png",
      path = dirImages,
      units = "cm", width = 11, height = 11
    )
  }

}
