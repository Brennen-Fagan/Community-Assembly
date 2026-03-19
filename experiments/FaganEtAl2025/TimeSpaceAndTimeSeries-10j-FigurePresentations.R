# Presentation variations:

variant <- c("Networks")[1]
figures <- c(3)
# 2:3 # Networks

if (variant == "Networks") {
  ### Networks oriented: ######################################################
  #### Setup: #################################################################
  ##### Resources: ############################################################
  # Data:
  source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
  source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
  source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

  load("TSTS_Interventions_10a1.RData")

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
      Time %in% figureNetworks$graph$time |                  # for examples
        (TimeSinceIntervention %in% figureNetworks$graph$timeInterventions)
    )

  figureNetworks$dataRich <- diversitiesRichness |> tidytable::filter(
    # SpeciesPreferences %in% figureNetworks$pref,
    NicheDistance == defaultNicheDistance,
    Intervention %in% figureNetworks$interventions,
    PoolPatchSeed %in% basePoolPatchSeeds,
    Metric == "Alpha Hill:0"
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
    )  + ggplot2::labs(
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
      # switch = "y",
      cols = ggplot2::vars(SpeciesPreferences)
    )

    ###### Inset Plots: #######################################################
    # Network plots at various times, and fuel for the KDEs.
    figureNetworks$graph$networks <- generateNetworks(
      figureNetworks$graph$specification,
      Date = "2025-07-30", split = TRUE
    )

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
        y = Size#, color = Intervention
      ),
      trim = TRUE, fill = "grey"
      # ) + ggplot2::geom_density(
      #   mapping = ggplot2::aes(
      #     y = Size, fill = Type, #color = Intervention
      #   ),
      #   alpha = 0.25,
      #   trim = TRUE
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


  ##### Figure 4: Richness, Abundance, Turnover, Complexity (RATC) ############

  #### 2 and Multiple Adaptation Type Figures: ################################
  ##### Figure 5: Richness w/Time, Abundance, Turnover, Complexity ############

  ##### Figure 6: Short Term Int. RATC ########################################

  #### Summary Images: ########################################################
  ##### Figure 7a: Parameters Cause RATC ######################################

  ##### Figure 7b: Network reorganisation over short time scales ##############

}
