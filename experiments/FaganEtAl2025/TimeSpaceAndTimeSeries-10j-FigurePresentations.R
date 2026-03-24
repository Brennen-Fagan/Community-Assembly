# Presentation variations:

variant <- c("Networks")[1]
figures <- c(4)
# 2:3 # Networks

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
  library(patchwork)

  # Directories:
  dirImages <- file.path(".", "TSTS_Images_Networks")
  if (!dir.exists(dirImages)) {
    dir.create(dirImages, showWarnings = FALSE)
  }

  ##### Data Management: ######################################################
  figureNetworks <- list(
    graph = list(
      # Note that, unusually, we need all max(time), but only seed for others.
      seed = "2", # "11", "17", "2"!,          # for examples
      time = c(100, 2000, 25000),              # for examples
      timeInterventions = c(0, 10, 100, 1000), # for examples
      pref = c("100% 0", "Uniform(0, 1)"),     # for KDEs
      interventions = c("(0)", "(0.5)", "(1)") # for KDEs
    ),
    interventions = c("(0)", "(0.5)->(0)", "(0.5)", "(0.5)->(1)", "(1)"),
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

  figureNetworks$dataRich <- diversitiesRichness |> tidytable::filter(
    # SpeciesPreferences %in% figureNetworks$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0"
  )

  figureNetworks$dataAbund <- diversitiesAbund |> tidytable::filter(
    # SpeciesPreferences %in% figure2$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Abundance"
  )

  figureNetworks$dataTurnover <- diversitiesTimeBC |> tidytable::filter(
    # SpeciesPreferences %in% figure2$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "TimeBrayCurtis: 10"
  )

  # Why to the level of summary? Because the PlotMeanAndInner function
  # isn't built to handle the multiple resolutions that we have in the
  # actual data, which makes it harder to portray the data accurately.
  figureNetworks$dataSummary <- figureNetworks$dataRich |>
    tidytable::left_join(
      InterventionTimes
    ) |> tidytable::mutate(
      Metric = factor(Metric, levels = c("Alpha Hill:0", "Alpha Abundance"),
                      labels = c("Richness", "Abundance"), ordered = TRUE),
      Time = round(Time - TimeIntervention, 6) # remove false differences
    ) |> tidytable::filter(
      is.na(Subset), # Overall values
      Time <= 16000, # Need the start for the inset.
      # Avoid singletons.
      abs(Time - round(Time)) < 1e-6 | Time >= 55 | Time < 0,
      Time > -1000,
      Metric %in% c("Richness", "Abundance")
    ) |> tidytable::mutate(
      Time = tidytable::case_when( # Create groupings for times.
        Time < -50 ~ round(Time, -2),
        Time < 0 ~ -25, # In the last bin before regime change.
        Time <= 50 ~ round(Time, 0),
        Time < 1105 ~ round(Time, -1), # Skip breaks < 5, drop.
        Time < 16350 ~ round(Time, -2),
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

  #### Cluster of Single Adaptation Type Figures: #############################
  if (2 %in% figures) {
    ##### Figure 2: Richness, Networks through Time ###########################
    ###### Main Plot: #########################################################
    # This is essentially figure 2, but with some more reference points, which
    # we will anchor additional network insets to.
    figureNetworks$plot2A <- plotMeanAndInner(
      figureNetworks$dataRich |> tidytable::filter(
        InterventionFinal == InterventionInitial,
        SpeciesPreferences == "100% 0",
        is.na(Subset)
      ) |> tidytable::mutate(
        SpeciesPreferences = tidytable::case_when(
          SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
          SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
          SpeciesPreferences == "Uniform(0, 1)" ~ "Multiple Adaptation Types",
          TRUE ~ SpeciesPreferences
        )
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
      # ) + ggplot2::scale_color_manual(
      #   values = colorPalette,
      #   name = "Habitat Land-use"
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
      ) |> tidytable::mutate(
        SpeciesPreferences = tidytable::case_when(
          SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
          SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
          SpeciesPreferences == "Uniform(0, 1)" ~ "Multiple Adaptation Types",
          TRUE ~ SpeciesPreferences
        )
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
          round(Time, 6) %in% figureNetworks$graph$timeInterventions
      ) |> tidytable::mutate(
        SpeciesPreferences = tidytable::case_when(
          SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
          SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
          SpeciesPreferences == "Uniform(0, 1)" ~ "Multiple Adaptation Types",
          TRUE ~ SpeciesPreferences
        )
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
      ) |> tidytable::mutate(
        SpeciesPreferences = tidytable::case_when(
          SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
          SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
          SpeciesPreferences == "Uniform(0, 1)" ~ "Multiple Adaptation Types",
          TRUE ~ SpeciesPreferences
        )
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
      x = "Time Since Intervention"#,
      # y = "Richness"
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
          group = paste(Intervention, Subset)
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
        # ) + ggplot2::facet_grid(
        #   . ~ SpeciesPreferences
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

    figureNetworks$plot4A <- (
      figureNetworks$dataRich |> tidytable::filter(
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
      transform = "pseudo_log", breaks = 10^(0:4),
      label = scales::label_log(digits = 2)
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
      )) |> figureNetworks$makeViolins()) + ggplot2::facet_wrap(ggplot2::vars(Metric), scales = "free"
      ) + ggplot2::theme(
        strip.text = ggplot2::element_text()
      )) |> ggplot2::ggsave(
        path = dirImages,
        filename = "FigureN4_Overview.png",
        units = "cm", width = 25, height = 14
      )
  }

  #### 2 and Multiple Adaptation Type Figures: ################################
  ##### Figure 5: Richness w/Time, Abundance, Turnover, Complexity ############

  ##### Figure 6: Short Term Int. RATC ########################################

  #### Summary Images: ########################################################
  ##### Figure 7a: Parameters Cause RATC ######################################

  ##### Figure 7b: Network reorganisation over short time scales ##############

}
