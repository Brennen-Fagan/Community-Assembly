# Now with added Invadability.

# Problems with X11
options(bitmapType = "cairo")

ribbonQuantiles <- c(0.05, 0.95) # For plotting vs time ribbons.
centralEst <- function(x) mean(x, na.rm = TRUE)
# For plotting vs time central point estimate.

# https://personal.sron.nl/~pault/#sec:qualitative
colorscheme <- c(
  # "red", "yellow", "cyan", "purple", "green", "grey"
  "#EE6677", "#CCBB44", "#66CCEE", "#AA3377", "#228833", "#BBBBBB"
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Plotting parameters.
time_grouping_size <- 20 # Used for initial filtering down of data, median.
time_averaging_size <- 4 # Used for averaging after grouping has been done.
time_burnin <- 1 # Discard data before this time for summarisation.
time_burnout <- 6.5 # Discard data after this time for summarisation.
time_units <- 1E4 # Inherited from diversity creation. Controls time axis label.
plot_width = 9; plot_height = 6; plot_units = "in"; plot_dpi = 640 # ggsave.

# Locations.
directory <- '.'
outputLocation <- file.path(directory, paste0("SaveImages_", Sys.Date()))
if (!dir.exists(outputLocation)) {
  dir.create(outputLocation, showWarnings = FALSE)
}

# Going to get a bit crafty here. Load into environments, then extract again.
temp <- lapply(1:4,
               function(i) {
                 load(file.path(directory, paste0("Set",i,"-AllLoadedDataBC.RData")))
                 return(list(
                   a = DiversitiesAlpha,
                   b = DiversitiesBeta,
                   g = DiversitiesGamma,
                   attr = Attributes
                 ))
               }
)

DiversitiesAlpha <- dplyr::bind_rows(temp[[1]]$a, temp[[2]]$a, temp[[3]]$a, temp[[4]]$a)
DiversitiesBeta <- dplyr::bind_rows(temp[[1]]$b, temp[[2]]$b, temp[[3]]$b, temp[[4]]$b)
DiversitiesGamma <- dplyr::bind_rows(temp[[1]]$g, temp[[2]]$g, temp[[3]]$g, temp[[4]]$g)
Attributes <- dplyr::bind_rows(temp[[1]]$attr, temp[[2]]$attr, temp[[3]]$attr, temp[[4]]$attr)

# Going to get a bit crafty here. Load into environments, then extract again.
temp <- lapply(1:4,
               function(i) {
                 load(file.path(directory, paste0("Set",i,"-AllInvadabilityData.RData")))
                 return(list(
                   i = Invadabilities,
                   attr = Attributes
                 ))
               }
)

Invadabilities <- dplyr::bind_rows(temp[[1]]$i, temp[[2]]$i, temp[[3]]$i, temp[[4]]$i)
AttrInv <- dplyr::bind_rows(temp[[1]]$attr, temp[[2]]$attr, temp[[3]]$attr, temp[[4]]$attr)

temp <- lapply(1:4,
               function(i) {
                 load(file.path(directory, paste0("Set",i,"-AlTimeJacData.RData")))
                 return(list(
                   j = TimeJaccards,
                   attr = Attributes
                 ))
               }
)

TimeJaccards <- dplyr::bind_rows(temp[[1]]$j, temp[[2]]$j, temp[[3]]$j, temp[[4]]$j)
AttrTimeJacs <- dplyr::bind_rows(temp[[1]]$attr, temp[[2]]$attr, temp[[3]]$attr, temp[[4]]$attr)

library(RMTRCode2) # source("Viking_HandleDiversity_HelperFunctions.R")
## Create Plotting Data Frames: ############################################
DiversitiesAlphaGamma <- dplyr::full_join(
  DiversitiesAlpha %>% dplyr::group_by(
    Time, Environment, Set, Number, History, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(dplyr::across(.fns = mean)),
  DiversitiesGamma %>% dplyr::group_by(
    Time, Set, Number, History, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(dplyr::across(.fns = mean)),
  # Why summarise? You'd expect each entry to be essentially unique, and you'd
  # be right in that they should be. Because we split up the runs into halves
  # and then binned them, this created duplicates at the edge of the runs.
  # Unfortunately, the edge then gets created twice with no method of telling
  # the two halves apart. To circumvent this issue, we run means, which does
  # nothing if there is a single entry, but resolves the case where there are 2.
  by = c(
    "Time", "Set", "Number", "History",
    "Pool", "Noise", "Neutral", "Space"
  ),
  suffix = c(", Alpha", ", Gamma")
)

## Create Helper Labels: ###################################################
levelsPoolNoise <- paste(
  DiversitiesGamma$Pool,
  DiversitiesGamma$Noise,
  sep = ";;"
) %>% unique() %>% strsplit(
  split = ";;", fixed = TRUE
)

spaceToDispersal <- data.frame(
  Space = unique(AttrInv$Space)
) %>% dplyr::mutate(
  Dispersal = 1 - exp( -2 / as.numeric(Space) ),
  Dispersal = paste0(
    formatC(Dispersal))
) %>% dplyr::arrange(
  as.numeric(Dispersal)
) %>% dplyr::mutate(
  # Edit to change the order of the factor levels
  Dispersal = factor(Dispersal, ordered = TRUE,
                     levels = Dispersal),
  Space = factor(Space, ordered = TRUE, levels = rev(Space))
)

neutralToLambdas <- data.frame(
  Neutral = unique(Attributes$Neutral)
) %>% tidyr::separate(
  Neutral, into = c("Immigration", "Extinction"), remove = FALSE,
  sep = "[,][ ]"
) %>% dplyr::mutate(
  # Edit to change the order of the factor levels
  Neutral = factor(Neutral, ordered = TRUE, levels = c(
    "0.1, 10",                    # 1:100
    "0.1, 1", "1, 10",            # 1:10
    "0.1, 0.1", "1, 1", "10, 10", # 1:1
    "10, 1", "1, 0.1",            # 10:1
    "10, 0.1",                    # 100:1
    "1, 0"                        # Inf:0
  )),
  Immigration = factor(Immigration, ordered = TRUE, levels = c(
    "0.1", "1", "10"
  )),
  Extinction = factor(Extinction, ordered = TRUE, levels = c(
    "10", "1", "0.1", "0"
  )),
  Immigration = paste0("Imm.: ", Immigration),
  Extinction = paste0("Ext.: ", Extinction)
)

FacetPanels <- c("1: Richness", "2: Jaccard", "3: Invadability")

concludingplotdata <- dplyr::bind_rows(
  DiversitiesAlphaGamma %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = "")
    # ) %>% tidyr::unite( # we do not need, since we'll group by statistic.
    #   col = "group", Neutral, Space, remove = FALSE
  ) %>% dplyr::group_by(
    Set, Number, History, Pool, Noise, Neutral, Space, Dispersal
  ) %>% dplyr::filter(
    Time > time_burnin, Time < time_burnout
  ) %>% dplyr::summarise(
    # Summarising over Time and Environment so 1 value per simulation.
    `Local Richness` = median(`Richness, Alpha`, na.rm = TRUE),
    `Regional Richness` = median(`Richness, Gamma`, na.rm = TRUE),
    .groups = "drop"
  ) %>% tidyr::pivot_longer(
    names_to = "Area",
    values_to = "Value",
    cols = c("Local Richness", "Regional Richness")
  ),
  Invadabilities %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = ""),
    Area = "Median Invasibility", #"Average Invadability",
    Value = LocalMdn # Choose LocalMdn for Medians, LocalAvg for Means.
  ) %>% dplyr::select(
    -Occupancy, -LocalAvg, -LocalMdn
  ),
  TimeJaccards %>% dplyr::filter(
    Measurement == "JaccardTemporal"
  ) %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = ""),
    Area = "Temporal Jaccard",
    Value = Mdn # Choose Mdn for Medians, Avg for Means.
  ) %>% dplyr::select(
    -Measurement, -Avg, -Mdn
  ),
  DiversitiesBeta %>% dplyr::left_join(
    spaceToDispersal, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = "")
  ) %>% dplyr::group_by(
    Set, Number, History, Pool, Noise, Neutral, Space, Dispersal
  ) %>% dplyr::filter(
    Time > time_burnin, Time < time_burnout
  ) %>% dplyr::summarise(
    # Summarising over Time and Environment so 1 value per simulation.
    JaccardSpatial = median(`Jaccard`, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::mutate(
    Area = "Spatial Jaccard",
    Value = JaccardSpatial
  ) %>% dplyr::select(
    -JaccardSpatial
  )
) %>% dplyr::distinct(
) %>% dplyr::mutate(
  FacetPanel = dplyr::case_when(
    Area == "Local Richness" ~ FacetPanels[1],
    Area == "Regional Richness" ~ FacetPanels[1],
    Area == "Temporal Jaccard" ~ FacetPanels[2],
    Area == "Spatial Jaccard" ~ FacetPanels[2],
    Area == "Median Invasibility" ~ FacetPanels[3],
    Area == "Average Invasibility" ~ FacetPanels[3]
  )
)

labelsforplot <- concludingplotdata %>% dplyr::select(
  Space, Dispersal
) %>% dplyr::distinct() %>% dplyr::arrange(Dispersal) %>% dplyr::mutate(
  Space = as.character(Space),
  Label = Dispersal#paste(Space, Dispersal, sep = "\n")
)

# Main Figure 4: ###############################################################

concludingplotdata_target <- concludingplotdata %>% dplyr::filter(
  Pool == "0, 0, 0, 0", # Noise == "1", # Incl. Het (Noise 1) and Hom (Noise 0)
  Neutral == "(1, 1)"
) %>% dplyr::mutate(
  Noise = dplyr::case_when(
    Noise == "2" ~ "Highly Het. Env.",
    Noise == "1" ~ "Heterogeneous Env.",
    Noise == "0" ~ "Homogeneous Env."
  ),
  Noise = factor(Noise, ordered = TRUE, levels = c(
    "Homogeneous Env.", "Heterogeneous Env.", "Highly Het. Env."
  ))
)

objs <- lapply(FacetPanels, function(Fac) {
  d <- concludingplotdata_target %>% dplyr::filter(FacetPanel == Fac)
  fills <- unique(d$Area)

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Dispersal, y = Value,
                 fill = Area,
                 group = interaction(Area, Dispersal))
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.8),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
    #title = title,
    #subtitle = paste0("Pool Modifier = ", pbs, ", ",
    #                  "Noise Modifier = ", inm)#,
    #y = ylabel
  ) + ggplot2::scale_x_discrete(
    name = #"Inter-community resistance and unit proportion dispersing",
      #"Proportion Dispersing in 1 Time Unit",
      "Dispersal Rate",# (Abundance Proportion)",
    # name = " ",
    #breaks = labelsforplot$Dispersal,
    #labels = labelsforplot$Label
    # if Square (width / 2, height * 2) or singlecolumn (width / 2)
    breaks = labelsforplot$Dispersal[c(1,2,4,6,8,10,11)],
    labels = c(as.character(labelsforplot$Label[c(1,2,4,6,8)]),
               "0.18   " ,"   0.86")
  ) + ggplot2::geom_line(
    data = d %>% dplyr::group_by(
      Area, Dispersal, FacetPanel, Noise
    ) %>% dplyr::summarise(
      Value = mean(Value, na.rm = TRUE), .groups = "drop"
    ),
    mapping = ggplot2::aes(group = Area, color = Area)
    # ) + ggplot2::scale_fill_viridis_d(
    #   name = "", begin = 0.35
    # ) + ggplot2::scale_color_viridis_d(
    #   name = "", begin = 0.35
  ) + ggplot2::scale_fill_manual(
    name = "", aesthetics = c("color", "fill"),
    values = c(
      "Local Richness" = colorscheme[1],
      "Regional Richness" = colorscheme[2],
      "Temporal Jaccard" = colorscheme[3],
      "Spatial Jaccard" = colorscheme[4],
      "Median Invasibility" = colorscheme[5],
      "Average Invasibility" = colorscheme[5]
    ),
    breaks = fills
  ) + ggplot2::theme(
    legend.position = #c(0.87, 0.93), Wide
      c(.226, 0.9), # Square (half width, twice height) or single column
    legend.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14)
  ) + ggplot2::ylab(
    if (Fac == FacetPanels[1]) {
      "Richness"
    } else if (Fac == FacetPanels[2]) {
      "Turnover"
    } else if (Fac == FacetPanels[3]) {
      "Invasibility"
    } else {
      "Value"
    }
  ) + ggplot2::facet_wrap(nrow = 1, facets = "Noise")

  if (Fac != FacetPanels[3]) {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    ) + ggplot2::expand_limits(
      y = 0
    )#+ ggplot2::xlab(" ")
  } else {
    obj <- obj + ggplot2::geom_text(
      data = data.frame(
      x = rep(c(1.8, 10.2), 2),
      y = rep(c(-14, -14), 2),  # if single column (half width)
      size = 4,
      label = c("     No Dispersal","", "", "Full Dispersal     ")
      ), ggplot2::aes(
        x = x, y = y, label = label
      ),
      inherit.aes = FALSE
    ) + ggplot2::coord_cartesian(
      clip = "off", ylim = c(0, 40)
    ) + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)
      )
  }

  obj
})

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "RichnessComparison-Invadability-BaseOnly.pdf"
  ),
  plot = (
    objs[[1]]/objs[[2]]/objs[[3]] +
      patchwork::plot_layout(ncol = 1)#, heights = c(3, 3, 3, 0.5))#, 1, 1))
  ),
  width = plot_width, height = plot_height * 1.25, units = plot_units,
  dpi = plot_dpi
)

# Supplemental Variations: #####################################################
## Grouped by Neutral: #########################################################
concludingplotdata_target <- concludingplotdata %>% dplyr::filter(
  Pool == "0, 0, 0, 0", Noise == "1"#, Neutral == "(1, 1)"
)

objs <- lapply(FacetPanels, function(Fac) {
  d <- concludingplotdata_target %>% dplyr::filter(
    FacetPanel == Fac
  ) %>% tidyr::separate(
    Neutral, sep = ", ", into = c("Imm.", "Ext.")
  ) %>% dplyr::mutate(
    Imm. = paste("Imm.:", substr(Imm., 2, nchar(Imm.))),
    Ext. = paste("Ext.:",substr(Ext., 1, nchar(Ext.) - 1))
  )
  fills <- unique(d$Area)

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Dispersal, y = Value,
                 fill = Area,
                 group = interaction(Area, Dispersal))
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.8),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
  ) + ggplot2::scale_x_discrete(
    name = #"Inter-community resistance and unit proportion dispersing",
      #"Proportion Dispersing in 1 Time Unit",
      "Dispersal Rate",
    # name = " ",
    #breaks = labelsforplot$Dispersal,
    #labels = labelsforplot$Label
    # if Square (width / 2, height * 2) or singlecolumn (width / 2)
    breaks = labelsforplot$Dispersal[c(1,2,4,6,8,10,11)],
    labels = c(as.character(labelsforplot$Label[c(1,2,4,6,8)]),
               "0.18   " ,"   0.86")
  ) + ggplot2::geom_line(
    data = d %>% dplyr::group_by(
      Area, Dispersal, FacetPanel, Imm., Ext.
    ) %>% dplyr::summarise(
      Value = mean(Value, na.rm = TRUE), .groups = "drop"
    ),
    mapping = ggplot2::aes(group = interaction(Area, Imm., Ext.), color = Area)
  ) + ggplot2::scale_fill_manual(
    name = "", aesthetics = c("color", "fill"),
    values = c(
      "Local Richness" = colorscheme[1],
      "Regional Richness" = colorscheme[2],
      "Temporal Jaccard" = colorscheme[3],
      "Spatial Jaccard" = colorscheme[4],
      "Median Invasibility" = colorscheme[5],
      "Average Invasibility" = colorscheme[5]
    ),
    breaks = fills
  ) + ggplot2::theme(
    legend.position = #c(0.87, 0.93), Wide
      c(0.125, 0.225), # Square (half width, twice height) or single column
    legend.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14)
  ) + ggplot2::ylab(
    if (Fac == FacetPanels[1]) {
      "Richness"
    } else if (Fac == FacetPanels[2]) {
      "Turnover"
    } else if (Fac == FacetPanels[3]) {
      "Invasibility"
    } else {
      "Value"
    }
  ) + ggplot2::facet_grid(Imm. ~ Ext.)

  if (Fac != FacetPanels[3]) {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  } else {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
  }

  obj
})

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "RichnessComparison-Invadability-ByNeutral.pdf"
  ),
  plot = (
    objs[[1]]/objs[[2]]/objs[[3]] +
      patchwork::plot_layout(ncol = 1)#, heights = c(3, 3, 3, 0.5))#, 1, 1))
  ),
  width = plot_width*1.5, height = plot_height*2.5, units = plot_units,
  dpi = plot_dpi
)

axesLimits_neutral <- concludingplotdata_target %>% dplyr::group_by(
  FacetPanel
) %>% dplyr::summarise(
  Min = min(Value),
  Max = max(Value)
)

## Grouped by Pool & Noise: ####################################################
concludingplotdata_target <- concludingplotdata %>% dplyr::filter(
  #Pool == "0, 0, 0, 1", Noise == "1",
  Neutral == "(1, 1)"
) %>% tidyr::unite(
  "PoolNoise", Pool, Noise, sep = " & "
) %>% dplyr::mutate(
  PoolNoise = dplyr::case_when(
    PoolNoise == "0, 0, 0, 0 & 0" ~ "Homogeneous Env.",
    PoolNoise == "0, 0, 0, 0 & 1" ~ "Base Case",
    PoolNoise == "0, 0, 0, 0 & 2" ~ "Highly Het. Env.",
    PoolNoise == "-1, 0, 0, 0 & 1" ~ "Smaller Basal",
    PoolNoise == "0, 0, 0, 1 & 1" ~ "Larger Consumer",
    PoolNoise == "-1, 0, 0, 1 & 1" ~ "Smaller and Larger",
    TRUE ~ "Oops"
  )
)

objs <- lapply(FacetPanels, function(Fac, Axes) {
  d <- concludingplotdata_target %>% dplyr::filter(
    FacetPanel == Fac
  ) %>% tidyr::separate(
    Neutral, sep = ", ", into = c("Imm.", "Ext.")
  ) %>% dplyr::mutate(
    Imm. = paste("Imm.:", substr(Imm., 2, nchar(Imm.))),
    Ext. = paste("Ext.:",substr(Ext., 1, nchar(Ext.) - 1))
  )
  fills <- unique(d$Area)

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Dispersal, y = Value,
                 fill = Area,
                 group = interaction(Area, Dispersal))
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.8),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
  ) + ggplot2::scale_x_discrete(
    name = #"Inter-community resistance and unit proportion dispersing",
      #"Proportion Dispersing in 1 Time Unit",
      "Dispersal Rate",
    # name = " ",
    #breaks = labelsforplot$Dispersal,
    #labels = labelsforplot$Label
    # if Square (width / 2, height * 2) or singlecolumn (width / 2)
    breaks = labelsforplot$Dispersal[c(1,2,4,6,8,10,11)],
    labels = c(
      "0  ", "  2e-09",
      as.character(labelsforplot$Label[c(4, 6)]),
      "2e-03", "0.18   " ,"  0.86")
  ) + ggplot2::scale_y_continuous(
    limits = unlist(Axes[Axes$FacetPanel == Fac, 2:3])
  ) + ggplot2::geom_line(
    data = d %>% dplyr::group_by(
      Area, Dispersal, FacetPanel, PoolNoise
    ) %>% dplyr::summarise(
      Value = mean(Value, na.rm = TRUE), .groups = "drop"
    ),
    mapping = ggplot2::aes(group = interaction(Area, PoolNoise), color = Area)
  ) + ggplot2::scale_fill_manual(
    name = "", aesthetics = c("color", "fill"),
    values = c(
      "Local Richness" = colorscheme[1],
      "Regional Richness" = colorscheme[2],
      "Temporal Jaccard" = colorscheme[3],
      "Spatial Jaccard" = colorscheme[4],
      "Median Invasibility" = colorscheme[5],
      "Average Invasibility" = colorscheme[5]
    ),
    breaks = fills
  ) + ggplot2::theme(
    legend.position = #c(0.87, 0.93), Wide
      c(0.125, 0.86), # Square (half width, twice height) or single column
    legend.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14)
  ) + ggplot2::ylab(
    if (Fac == FacetPanels[1]) {
      "Richness"
    } else if (Fac == FacetPanels[2]) {
      "Turnover"
    } else if (Fac == FacetPanels[3]) {
      "Invasibility"
    } else {
      "Value"
    }
  ) + ggplot2::facet_grid(. ~ PoolNoise)

  if (Fac != FacetPanels[3]) {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  } else {
    obj <- obj + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
  }

  obj
}, Axes = axesLimits_neutral)

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    "RichnessComparison-Invadability-ByPoolNoise.pdf"
  ),
  plot = (
    objs[[1]]/objs[[2]]/objs[[3]] +
      patchwork::plot_layout(ncol = 1)#, heights = c(3, 3, 3, 0.5))#, 1, 1))
  ),
  width = plot_width*1.5, height = plot_height*1.75, units = plot_units,
  dpi = plot_dpi
)

## Alpha Gamma Space Plots: ####################################################
concludingplotdata_target <- DiversitiesAlphaGamma%>% dplyr::filter(
  Space %in% c("1", "1000", "1e+06", "1e+09"),
  Neutral %in% c("1, 0.1", "1, 1", "0.1, 1"),
  Pool == "0, 0, 0, 0", Noise == "1"
)

plotsGammaOverAlpha <- lapply(
  1,#seq_along(levelsPoolNoise),
  GammaOverAlpha,
  d = concludingplotdata_target,
  levels = levelsPoolNoise,
  spacelabs = spaceToDispersal,
  neutrallabs = neutralToLambdas
)

ggplot2::ggsave(
  filename = file.path(
    outputLocation,
    paste0("GammaOverAlpha-BaseOnly.pdf")
  ),
  plot = plotsGammaOverAlpha[[1]] + ggplot2::labs(
    x = "Local Richness", y = "Regional Richness",
    title = NULL, subtitle = NULL
  ) + ggplot2::facet_grid(
    Immigration + factor(Extinction, ordered = TRUE,
                         levels = c("Ext.: 1", "Ext.: 0.1"))
    ~ factor(
      Dispersal,
      levels = c(2e-09, 2e-06, 0.001998, 0.8647),
      labels = c("Little Dispersal", "< Medium Dispersal",
                 "> Medium Dispersal", "Full Dispersal"), ordered = TRUE
    )
  ) + ggplot2::geom_text(
    data = data.frame(
      x = c(45, 45),
      y = c(25, 25),
      label = c("Reduced\nImmigration\n(from pool)", "Reduced\nStochastic\nExtirpation"),
      Immigration = c("Imm.: 0.1", "Imm.: 1"),
      Extinction = c("Ext.: 1", "Ext.: 0.1"),
      Dispersal = c(0.8647, 0.8647)
    ),
    mapping = ggplot2::aes(x = x, y = y, label = label)
  ) + ggplot2::coord_cartesian(
    xlim = range(concludingplotdata_target$`Richness, Alpha`),
    ylim = range(concludingplotdata_target$`Richness, Gamma`),
    clip = FALSE
  ),
  width = plot_width, height = plot_height, units = plot_units,
  dpi = plot_dpi
)
