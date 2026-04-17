# Setup: ######################################################################
# Plot of Richness as a function of species preferences and land-use,
# when species preferences are 100% 0.
# Also functinally an overview plot of network structure.

source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsPersistence.R")
source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")

library(cowplot) # plot arrangement with pagination (ggarrange).
library(tibble) # data.frames that can incorporate grobs.
library(ggpp) # for the geom_grob function -> factor placement of grobs.

# This is better as an environment, but that's more opaque.
figure5 <- list(
  pref = c("Uniform(0, 1)", "50% 0, 50% 1"),
  isSupplement = TRUE
)

# Main Plots: #################################################################
### Plot 5: ###################################################################
##### Data: ###################################################################
# Richness data: should be straightforward.
figure5$dataRich <- diversitiesRichness |> tidytable::filter(
  SpeciesPreferences %in% figure5$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Hill:0"
)
figure5$dataAbund <- diversitiesAbund |> tidytable::filter(
  SpeciesPreferences %in% figure5$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds,
  Metric == "Alpha Abundance"
)

# Persistence data: why? Because we're using persistence as a weight, followed
# with by species aggregation. That way approximately we are picking a random
# simulation, a random time, and then a random species, the plot shows the
# probability we would get a certain land-use preference out.
figure5$dataPers <- Pers |> tidytable::filter(
  SpeciesPreferences %in% figure5$pref,
  NicheDistance == defaultNicheDistance,
  Intervention %in% c("(0)", "(0.25)", "(0.5)", "(0.75)", "(1)"),
  PoolPatchSeed %in% basePoolPatchSeeds
) |> tidytable::filter(
  In < Stop, Out > Start # Not things outside of [Start, Stop]
) |> tidytable::mutate(
  # Shorten intervals for equivalent comparisons.
  InType = ifelse(In < Start, "Persistent", InType),
  OutType = ifelse(Out > Stop, "Persistent", OutType),
  In = ifelse(In < Start, Start, In),
  Out = ifelse(Out > Stop, Stop, Out),
  Persistence = Out - In
) |> tidytable::group_by(
  Species, Environment, SpeciesType, Size, ReproductionRate, Speed,
  Affinity, AffinityBins,
  PoolPatch:InterventionNicheDistance,
  Intervention, SpeciesPreferences, Start, Stop
) |> tidytable::summarise( # Sum over Appearances.
  Persistence = sum(Persistence),
  .groups = "drop"
)

##### a: ######################################################################
# Richness through time across simulations, showing stability and separation.
figure5$plotA <- lapply(
  unique(figure5$dataRich$SpeciesPreferences),
  function(sp) {
    plotMeanAndInner(
      rbind(
        figure5$dataRich |> tidytable::filter(
          Intervention %in% c("(0)", "(0.5)", "(1)"),
          is.na(Subset), SpeciesPreferences == sp
        ),
        # We want to appear in the legend but not on the plot!
        figure5$dataRich |> tidytable::filter(
          PoolPatchSeed == figure5$dataRich$PoolPatchSeed[1],
          Intervention %in% c("(0.25)", "(0.75)"),
          abs(Time - 10000) == min(abs(Time - 10000)),
          is.na(Subset), SpeciesPreferences == sp
        ) |> tidytable::mutate(
          Value = 10 # coord_cartesian will eliminate these points.
        )
      ) |> tidytable::mutate(
        SpeciesPreferences = tidytable::case_when(
          SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
          SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
          SpeciesPreferences == "Uniform(0, 1)" ~ "Multiple Adaptation Types",
          TRUE ~ SpeciesPreferences
        )
      ), CIs = 0.75
    ) + ggplot2::labs(
      y = "Richness",
      tag = "a)"
    ) + ggplot2::guides(
      linetype = "none",
      color = ggplot2::guide_legend(ncol = 5,
                                    title = "Habitat Type"),
      fill = ggplot2::guide_legend(ncol = 5,
                                   title = "Habitat Type")
    ) + ggplot2::theme(
      legend.position = c(0.5, 0.09),
      plot.tag.position = c(0.025, 0.95),
      panel.spacing = ggplot2::unit(1, "lines"),
      strip.text = ggplot2::element_text(size = 12)
    ) + ggplot2::coord_cartesian(
      ylim = c(0, 30), expand = FALSE
    ) + ggplot2::scale_x_continuous(
      breaks = (0:3)*10000
    ) + ggplot2::facet_wrap( # Instead of grid to have strips on top.
      # switch = "y",
      ggplot2::vars(SpeciesPreferences), ncol = 1
    )
  }
)

##### b: Violins ##############################################################
figure5$plotB <- lapply(
  unique(figure5$dataRich$SpeciesPreferences),
  function(sp) {
    ggplot2::ggplot(
      figure5$dataRich |> tidytable::filter(
        Time > Start, Time < Stop, SpeciesPreferences == sp
      ) |> tidytable::group_by(
        PoolPatchSeed, Intervention, SpeciesPreferences, Subset
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
      width = 0.28
    ) + ggplot2::scale_color_manual(
      values = colorPalette, aesthetics = c("color", "fill"),
      name = "Habitat Type"
    ) + ggplot2::scale_y_continuous(
      breaks = c(0, 10, 20, 30)
    ) + ggplot2::coord_cartesian(
      ylim = c(0, 30), expand = FALSE
    ) + ggplot2::facet_wrap( # Instead of grid to have strips on top.
      # switch = "y",
      ggplot2::vars(SpeciesPreferences), ncol = 1
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      plot.tag.position = c(0.01, 1),
      panel.grid.minor = ggplot2::element_blank(),
      # strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(color = NA),
      panel.spacing = ggplot2::unit(1, "lines")
    ) + ggplot2::labs(
      y = "Richness",
      x = "Habitat Type"
    ) + ggplot2::guides(
      color = "none",
      fill = "none"
    )
  }
)

figure5$plotC <- lapply(
  unique(figure5$dataAbund$SpeciesPreferences),
  function(sp) {
    ggplot2::ggplot(
      figure5$dataAbund |> tidytable::filter(
        Time > Start, Time < Stop, SpeciesPreferences == sp
      ) |> tidytable::group_by(
        PoolPatchSeed, Intervention, SpeciesPreferences, Subset
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
      width = 0.28
    ) + ggplot2::scale_color_manual(
      values = colorPalette, aesthetics = c("color", "fill"),
      name = "Habitat Type"
    ) + ggplot2::facet_wrap( # Instead of grid to have strips on top.
      # switch = "y",
      ggplot2::vars(SpeciesPreferences), ncol = 1
    ) + ggplot2::theme_minimal(
    ) + ggplot2::theme(
      plot.tag.position = c(0.01, 1),
      panel.grid.minor = ggplot2::element_blank(),
      # strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(color = NA),
      panel.spacing = ggplot2::unit(1, "lines")
    ) + ggplot2::labs(
      y = "Abundance",
      x = "Habitat Type"
    ) + ggplot2::guides(
      color = "none",
      fill = "none"
    ) + scale_y_log10(
      label = scales::label_log(),
      breaks = 10^(2+(0:20)/5)
    )
  }
)

if (figure5$isSupplement) {
  figure5$plotB <- lapply(
    figure5$plotB,
    function(plt) {
      plt + geom_flat_violin(
        data = function(x) x |> tidytable::filter(
          !is.na(Subset)
        ) |> tidytable::separate_wider_delim(
          cols = Subset, delim = "_", names = c("Subset", "Bins")
        ) |> tidytable::filter(
          Subset == "Basal"
        ) |> tidytable::group_by(
          PoolPatchSeed, Intervention, SpeciesPreferences, Subset
        ) |> tidytable::summarise(
          Value = sum(Value)
        ),
        flip = 1, color = "black",
        ggplot2::aes(fill = Intervention)
        # CONSUMER Violin
      ) + geom_flat_violin(
        data = function(x) x |> tidytable::filter(
          !is.na(Subset)
        ) |> tidytable::separate_wider_delim(
          cols = Subset, delim = "_", names = c("Subset", "Bins")
        ) |> tidytable::filter(
          Subset == "Consumer"
        ) |> tidytable::group_by(
          PoolPatchSeed, Intervention, SpeciesPreferences, Subset
        ) |> tidytable::summarise(
          Value = sum(Value)
        ),
        flip = 2, color = "black",
        ggplot2::aes(fill = Intervention)
        # LABELS
      ) + ggplot2::geom_text(
        data = function(x) {
          x |> tidytable::filter(
            !is.na(Subset)
          ) |> tidytable::separate_wider_delim(
            cols = Subset, delim = "_", names = c("Subset", "Bins")
          ) |> tidytable::group_by(
            PoolPatchSeed, Intervention, SpeciesPreferences, Subset
          ) |> tidytable::summarise(
            Value = sum(Value),
            .groups = "drop"
          ) |> tidytable::mutate(
            Offset = (max(Value) - min(Value)) * 0.08
          ) |> dplyr::group_by(
            Intervention, SpeciesPreferences, Subset
          ) |> dplyr::summarise(
            Value = min(Value - Offset),
            Label = substr(Subset[1], 0, 1),
            .groups = "drop"
          )
        },
        mapping = ggplot2::aes(label = Label),
        position = ggplot2::position_dodge(0.9)
      )
    }
  )

  figure5$plotC <- lapply(
    figure5$plotC,
    function(plt) {
      plt + geom_flat_violin(
        # BASAL VIOLINS
        data = function(x) x |> tidytable::filter(
          !is.na(Subset)
        ) |> tidytable::separate_wider_delim(
          cols = Subset, delim = "_", names = c("Subset", "Bins")
        ) |> tidytable::filter(
          Subset == "Basal"
        ) |> tidytable::group_by(
          PoolPatchSeed, Intervention, SpeciesPreferences, Subset
        ) |> tidytable::summarise(
          Value = sum(Value)
        ),
        flip = 1, color = "black",
        ggplot2::aes(fill = Intervention)
      ) + geom_flat_violin(
        # CONSUMER VIOLINS
        data = function(x) x |> tidytable::filter(
          !is.na(Subset)
        ) |> tidytable::separate_wider_delim(
          cols = Subset, delim = "_", names = c("Subset", "Bins")
        ) |> tidytable::filter(
          Subset == "Consumer"
        ) |> tidytable::group_by(
          PoolPatchSeed, Intervention, SpeciesPreferences, Subset
        ) |> tidytable::summarise(
          Value = sum(Value)
        ),
        flip = 2, color = "black",
        ggplot2::aes(fill = Intervention)
        # LABELS
      ) + ggplot2::geom_text(
        data = function(x) {
          x |> tidytable::filter(
            !is.na(Subset)
          ) |> tidytable::separate_wider_delim(
            cols = Subset, delim = "_", names = c("Subset", "Bins")
          ) |> tidytable::group_by(
            PoolPatchSeed, Intervention, SpeciesPreferences, Subset
          ) |> tidytable::summarise(
            Value = sum(Value),
            .groups = "drop"
          ) |> tidytable::mutate(
            Offset = (max(log10(Value)) - min(log10(Value))) * 0.08
          ) |> dplyr::group_by(
            Intervention, SpeciesPreferences, Subset
          ) |> dplyr::summarise(
            Value = 10^min(log10(Value) - Offset),
            Label = substr(Subset[1], 0, 1),
            .groups = "drop"
          )
        },
        mapping = ggplot2::aes(label = Label),
        position = ggplot2::position_dodge(0.9)
      )
    }
  )
}

##### c: Insets ###############################################################
figure5$insetGrobs <- figure5$dataPers |> tidytable::filter(
  # Match A
  Intervention %in% c("(0)", "(0.5)", "(1)")
) |> tidytable::mutate(
  SpeciesPreferences = tidytable::case_when(
    SpeciesPreferences == "100% 0" ~ "1 Adaptation Type",
    SpeciesPreferences == "50% 0, 50% 1" ~ "2 Adaptation Types",
    SpeciesPreferences == "Uniform(0, 1)" ~ "Multiple Adaptation Types",
    TRUE ~ SpeciesPreferences
  )
) |> dplyr::group_by(
  SpeciesPreferences, Intervention
) |> dplyr::group_modify(
  .f = function(data, key) {
    if (key$SpeciesPreferences == "2 Adaptation Types") {
      # Image:
      tibblegrob <- tibble::tibble(
        grob = list(
          (
            dplyr::bind_cols(data, key) |> tidytable::group_by(
              Intervention, AffinityBins, SpeciesType
            ) |> tidytable::summarise(
              Persistence = sum(Persistence)
            ) |> tidytable::group_by(
              Intervention
            ) |> tidytable::mutate(
              Persistence = Persistence / sum(Persistence)
            ) |> ggplot2::ggplot(
              mapping = ggplot2::aes(
                x = AffinityBins,
                y = Persistence,
                color = Intervention,
                fill = factor(SpeciesType, ordered = TRUE,
                              levels = c("Consumer", "Basal"))
              )
            ) + ggplot2::geom_col(
              show.legend = FALSE
            ) + ggplot2::facet_grid(
              . ~ Intervention
            ) + ggplot2::scale_color_manual(
              values = colorPalette,
              name = "Habitat Type"
            ) + ggplot2::scale_fill_manual(
              values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
            ) + ggplot2::theme_void(
            ) + ggplot2::theme(
              plot.background = ggplot2::element_rect(fill = NA,
                                                      color = NA),
              panel.background = ggplot2::element_rect(fill = NA),
              strip.text = ggplot2::element_text(size = 10)
            ) + ggplot2::coord_cartesian(
              expand = FALSE, ylim = c(0, 1)
            )
          ) |> ggplot2::ggplotGrob()
        )
      )

      # Coordinates:
      if (key$Intervention == "(0)") {
        tibbleplace <- tibble::tibble(
          Time = 7500, Value = 2.5
        )
      } else if (key$Intervention == "(0.25)") {
        tibbleplace <- tibble::tibble(
          Time = 12500, Value = 2.5
        )
      } else if (key$Intervention == "(0.5)") {
        tibbleplace <- tibble::tibble(
          Time = 17500, Value = 2.5
        )
      } else if (key$Intervention == "(0.75)") {
        tibbleplace <- tibble::tibble(
          Time = 22500, Value = 2.5
        )
      } else if (key$Intervention == "(1)") {
        tibbleplace <- tibble::tibble(
          Time = 27500, Value = 2.5
        )
      }
    } else if (key$SpeciesPreferences == "Multiple Adaptation Types") {
      # Image:
      tibblegrob <- tibble::tibble(
        grob = list(
          (
            dplyr::bind_cols(data, key) |> ggplot2::ggplot(
              mapping = ggplot2::aes(
                x = Affinity,
                weight = Persistence,
                fill = factor(SpeciesType, ordered = TRUE,
                              levels = c("Consumer", "Basal")),
                # group = interaction(Intervention, SpeciesType),
                color = Intervention
              )
            ) + ggplot2::geom_density(
              position = "stack",
              adjust = 1/2, linewidth = 0.5,
              show.legend = FALSE
            ) + ggplot2::facet_grid(
              . ~ Intervention
            ) + ggplot2::scale_color_manual(
              values = colorPalette,
              name = "Habitat Type"
            ) + ggplot2::scale_fill_manual(
              values = c("Basal" = "darkgreen", "Consumer" = "burlywood4")
            ) + ggplot2::theme_void(
            ) + ggplot2::theme(
              plot.background = ggplot2::element_rect(fill = NA,
                                                      color = NA),
              panel.background = ggplot2::element_rect(fill = NA),
              strip.text = ggplot2::element_text(size = 10)
            ) + ggplot2::coord_cartesian(
              expand = FALSE
            )
          ) |> ggplot2::ggplotGrob()
        )
      )

      # Coordinates:
      if (key$Intervention == "(0)") {
        tibbleplace <- tibble::tibble(
          Time = 7500, Value = 22.5
        )
      } else if (key$Intervention == "(0.25)") {
        tibbleplace <- tibble::tibble(
          Time = 12500, Value = 22.5
        )
      } else if (key$Intervention == "(0.5)") {
        tibbleplace <- tibble::tibble(
          Time = 17500, Value = 22.5
        )
      } else if (key$Intervention == "(0.75)") {
        tibbleplace <- tibble::tibble(
          Time = 22500, Value = 22.5
        )
      } else if (key$Intervention == "(1)") {
        tibbleplace <- tibble::tibble(
          Time = 27500, Value = 22.5
        )
      }
    }

    tibblegrob <- dplyr::bind_cols(tibblegrob, tibbleplace)

    return(tibblegrob)
  }
)

##### Combine: ################################################################

figure5$plotA[[2]] <- figure5$plotA[[2]] + ggplot2::labs(tag = "a)")
figure5$plotB[[2]] <- figure5$plotB[[2]] + ggplot2::labs(tag = "b)")
figure5$plotC[[2]] <- figure5$plotC[[2]] + ggplot2::labs(tag = "c)")
figure5$plotA[[1]] <- figure5$plotA[[1]] + ggplot2::labs(tag = "d)")
figure5$plotB[[1]] <- figure5$plotB[[1]] + ggplot2::labs(tag = "e)")
figure5$plotC[[1]] <- figure5$plotC[[1]] + ggplot2::labs(tag = "f)")

figure5$plot <- ggpubr::ggarrange(
  plotlist = list(
    figure5$plotA[[2]] + ggpp::geom_grob(
      data = figure5$insetGrobs |> tidytable::filter(
        SpeciesPreferences == "2 Adaptation Types"
      ),
      mapping = ggplot2::aes(x = Time, y = Value, label = grob)
    ),
    figure5$plotB[[2]],
    figure5$plotC[[2]],
    figure5$plotA[[1]] + ggpp::geom_grob(
      data = figure5$insetGrobs |> tidytable::filter(
        SpeciesPreferences == "Multiple Adaptation Types"
      ),
      mapping = ggplot2::aes(x = Time, y = Value, label = grob)
    ),
    figure5$plotB[[1]],
    figure5$plotC[[1]]
  ),
  nrow = 2, ncol = 3, widths = c(1/3, 1/3, 1/3), common.legend = TRUE
)

figure5$prefix <- "Figure5"
if (figure5$isSupplement) {
  figure5$affix <- "S_Prototype3"
} else {
  figure5$affix <- "_Prototype3"
}

with(
  figure5,
  mapply(
    list(plot, plot, plotA[[2]], plotB[[2]], plotC[[2]],
         plotA[[1]], plotB[[1]], plotC[[1]],
         ggplot2::ggplot(
         ) + ggplot2::coord_cartesian(
           xlim = c(3500, 30000), ylim = c(0, 27)
         ) + ggpp::geom_grob(
           data = figure5$insetGrobs,
           mapping = ggplot2::aes(x = Time, y = Value, label = grob)
         )),
    c("", "", "A", "B", "C", "D", "E", "F", "Inset"),
    c(".png", rep(".pdf", 8)),
    FUN = function(arg1, arg2, arg3) {
      ggplot2::ggsave(
        plot = arg1,
        filename = file.path(dirImages, paste0(prefix, arg2, affix, arg3)),
        units = "cm", width = 6.5*3, height = 6.5*3)
    }
  )
)
