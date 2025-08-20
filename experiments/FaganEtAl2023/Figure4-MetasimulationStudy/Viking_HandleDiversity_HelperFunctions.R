# Helper functions for Viking_HandleDiversity_ParametersAndPlots3.R (v3+)


# Functions: ###################################################################
## Attributes: ################################################################
extractAttributes <- function(Diversity, idNums) {
  data.frame(
    Pool = toString(Diversity$PoolMod),
    Noise = toString(Diversity$NoiseMod),
    Neutral = toString(Diversity$NeutralMod),
    Space = toString(Diversity$SpaceMod),
    Set = idNums[1],
    CaseNumber = idNums[2],
    History = if(is.na(idNums[3])) {
      # All Histories
      1:10
    } else idNums[3],
    Part = if(is.na(idNums[4])) {
      # All Parts
      1
    } else idNums[4],
    stringsAsFactors = FALSE
  )
}

assignAttributes <- function(i, d, a) {
  d %>% dplyr::bind_rows(
  ) %>% dplyr::mutate(
    Set = a$Set[i],
    Number = a$CaseNumber[i],
    History = a$History[i],
    Pool = a$Pool[i],
    Noise = a$Noise[i],
    Neutral = a$Neutral[i],
    Space = a$Space[i]
  )
}

## Diversity Metrics: ##########################################################
# Primary distinctions are what metric is extracted in function(i, d, a),
# what columns might be dropped, what variables are used to group, and flooring.
extractAlphas <- function(Diversity, Attributes) {
  dplyr::bind_rows(
    lapply(
      seq_along(Diversity$Diversities),
      function(i, d, a) assignAttributes(i, d[[i]]$alpha, a),
      d = Diversity$Diversities,
      a = Attributes
    )
  ) %>% dplyr::select(-Species) %>% dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) %>% dplyr::group_by(
    Time, Environment, Set, Number, History, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    Richness = floor(median(Richness)),
    .groups = "drop"
  )
}

extractBetas <- function(Diversity, Attributes) {
  dplyr::bind_rows(
    lapply(
      seq_along(Diversity$Diversities),
      function(i, d, a) assignAttributes(i, d[[i]]$beta, a),
      d = Diversity$Diversities,
      a = Attributes
    )
  ) %>% dplyr::mutate(
    Time = floor(Time * time_grouping_size)/time_grouping_size
  ) %>% dplyr::group_by(
    Time, Env1, Env2, Set, Number, History, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    Jaccard = median(Jaccard),
    .groups = "drop"
  )
}

extractGammas <- function(Diversity, Attributes) {
  dplyr::bind_rows(lapply(
    seq_along(Diversity$Diversities),
    function(i, d, a) assignAttributes(i, d[[i]]$gamma, a),
    d = Diversity$Diversities,
    a = Attributes)) %>% dplyr::filter(
      Aggregation == "Gamma"
    ) %>% dplyr::mutate(
      Time = floor(Time * time_grouping_size)/time_grouping_size
    ) %>% dplyr::group_by(
      Time, Set, Number, History, Pool, Noise, Neutral, Space
    ) %>% summarise(
      Richness = floor(median(Richness)),
      .groups = "drop"
    )
}

## Load Function: ##############################################################
loadDiversity <- function(file, verbose = TRUE) {
  # Retrieve the trailing id numbers from before the file extension.
  # 1st: Set, 2nd: CaseNumber, 3rd: History, 4th: Part
  # Note, if last two are not present, all histories are bundled together.
  # If last two are present, each file has a single part of a single history.
  idNums <- suppressWarnings(na.omit(
    as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                             split = "-", fixed = TRUE)[[1]],
                    n = 4))),
    classes = "simpleWarning")

  if (verbose) {
    print(file)
    print(Sys.time())
  }

  tryCatch({
    load(file)

    if (!exists("Diversity")) {stop("Diversity does not exist.")}

    Attributes <- extractAttributes(Diversity, idNums)
    DiversitiesAlpha <- extractAlphas(Diversity, Attributes)
    DiversitiesBeta <- extractBetas(Diversity, Attributes)
    DiversitiesGamma <- extractGammas(Diversity, Attributes)

    return(list(
      Attr = Attributes,
      Alpha = DiversitiesAlpha,
      Beta = DiversitiesBeta,
      Gamma = DiversitiesGamma
    ))
  }, error = function(e) {
    print(paste(paste(idNums, collapse = " "), e, sep = ": "))
    return(NULL)
  })
}

## Average Function: ###########################################################
createAvg <- function(d, target, ps = ribbonQuantiles, cE = centralEst) {
  d %>% dplyr::mutate(
    Time = floor(Time * time_averaging_size)/time_averaging_size
  ) %>%  dplyr::group_by(
    Time, Pool, Noise, Neutral, Space
  ) %>% dplyr::summarise(
    Low = unlist(dplyr::across(dplyr::any_of(target),
                               .fns = ~ quantile(.x, p = ps[1], na.rm = TRUE))),
    High = unlist(dplyr::across(dplyr::any_of(target),
                                .fns = ~ quantile(.x, p = ps[2], na.rm = TRUE))),
    Central = unlist(dplyr::across(dplyr::any_of(target),
                                   .fns = ~ cE(.x))),
    .groups = "drop"
  )
}

## Plotting Functions: #########################################################
# Format is [Y-axis]Over[X-axis].
StatisticOverTime <- function(
  i, d, da, levels,
  spacelabs, neutrallabs,
  ytarget, ylabel, title,
  burnin = time_burnin,
  burnout = time_burnout,
  units = time_units
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]

  # Attach labels.
  d <- d %>% dplyr::filter(
    pbs == Pool,
    inm == Noise
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  ) %>% tidyr::unite(
    col = "group", !Time & !dplyr::all_of(ytarget), remove = FALSE
  )

  da <- da %>% dplyr::filter(
    pbs == Pool,
    inm == Noise
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  ) %>% tidyr::unite(
    col = "group", !Time & !Low & !High & !Central, remove = FALSE
  )

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = Time,
                 group = group)
  ) + ggplot2::geom_line(
    mapping = ggplot2::aes_string(
      y = ytarget
    ),
    alpha = 0.05
  ) + ggplot2::geom_line(
    data = da,
    mapping = ggplot2::aes(
      y = Central
    ),
    color = "red", size = 1.5
  ) + ggplot2::geom_ribbon(
    data = da,
    mapping = ggplot2::aes(
      ymin = Low,
      ymax = High,
      x = Time
    ),
    alpha = 0.4,
    fill = "blue",
    inherit.aes = FALSE
  ) + ggplot2::geom_vline(
    xintercept = c(burnin, burnout)
  ) + ggplot2::theme_bw(
  ) + ggplot2::facet_grid(
    Immigration + Extinction ~ Space + Dispersal,
    labeller = ggplot2::label_parsed
  ) + ggplot2::labs(
    title = title,
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm),
    ylab = ylabel
  ) + ggplot2::xlab(
    bquote("Time, " * 10^.(log10(time_units)) * " Units")
  )

  return(obj)
}

StatisticOverSpace <- function(
  i, d, levels,
  spacelabs,
  # ytarget, SET TO Richness
  ylabel, title,
  burnin = time_burnin,
  burnout = time_burnout
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]

  # Attach labels.
  d <- d %>% dplyr::filter(
    pbs == Pool,
    inm == Noise,
    Time > burnin, Time < burnout
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::mutate(
    Neutral = paste("(", Neutral, ")", sep = "")
  ) %>% tidyr::unite(
    col = "group", Neutral, Space, remove = FALSE
  )

  xlabs <- d %>% dplyr::select(
    Space, Dispersal
  ) %>% dplyr::distinct() %>% dplyr::arrange(Space) %>% dplyr::mutate(
    Space = as.character(Space),
    Label = paste(Space, Dispersal, sep = "\n")
  )

  obj <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = as.character(Space), y = Richness,
                 fill = Neutral,
                 group = group)
  ) + ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.5), notch = TRUE
  ) + ggplot2::stat_summary(
    fun = mean,
    position = ggplot2::position_dodge(0.5),
    shape = 4
  ) + ggplot2::theme_bw(
  ) + ggplot2::labs(
    title = title,
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm),
    y = ylabel
  ) + ggplot2::scale_x_discrete(
    name = "Inter-community distance and unit proportion dispersing",
    breaks = xlabs$Space,
    labels = xlabs$Label
  ) + ggplot2::scale_fill_discrete(
    name = "(Imm., Ext.)\nModifiers"
  )

  return(obj)
}

GammaOverAlpha <- function(
  i, d, levels,
  spacelabs, neutrallabs,
  burnin = time_burnin,
  burnout = time_burnout,
  # For each parameter set, there are 10 Pools, and in each...
  # Pool (for which there are 10 Histories),
  # History (for which there are 10 Environments), Environment
  NumberSelect = NULL ,#1,
  HistorySelect = NULL, #1,
  EnvirSelect = NULL #1
) {
  # Different levels need different plots; not directly comparable.
  pbs <- levels[[i]][1]
  inm <- levels[[i]][2]

  d <- d %>% dplyr::filter(
    Pool == pbs,
    Noise == inm
  ) %>% dplyr::left_join(
    spacelabs, by = "Space"
  ) %>% dplyr::left_join(
    neutrallabs, by = "Neutral"
  )

  if (!is.null(NumberSelect) &&
      !is.null(HistorySelect) &&
      !is.null(EnvirSelect)) {
    number <- d %>% dplyr::select(
      Space, Neutral, Number, History, Environment
    ) %>% dplyr::distinct() %>% dplyr::arrange(
      Number
    ) %>% dplyr::group_by(
      Space, Neutral
    ) %>% dplyr::filter(
      Number %in% (min(Number) - 1 + NumberSelect)
    ) %>% dplyr::group_by(
      Space, Neutral, Number
    ) %>% dplyr::filter(
      History %in% (min(History) - 1 + HistorySelect)
    ) %>% dplyr::group_by(
      Space, Neutral, Number, History
    ) %>% dplyr::filter(
      Environment %in% (min(Environment) - 1 + EnvirSelect)
    )

    dTraj <- d %>% dplyr::semi_join(number, by = c(
      "Space", "Neutral", "Number", "History", "Environment"
    ))

    dTrajStop <- dTraj %>% dplyr::group_by(
      Number, History, Environment
    ) %>% dplyr::slice_max(Time)

  } else if (!is.null(NumberSelect) ||
             !is.null(HistorySelect) ||
             !is.null(EnvirSelect)) {
    warning("Not all of Selects are null or non-null. Skipping trajectory.")
  }

  dga <- d %>% dplyr::filter(
    Time > burnin, Time < burnout
  )

  obj <- ggplot2::ggplot(
    dga,
    ggplot2::aes(x = `Richness, Alpha`,
                 y = `Richness, Gamma`)
  ) + ggplot2::geom_bin2d(
    breaks = seq(-0.5, max(dga$`Richness, Alpha`,
                           dga$`Richness, Gamma`) + 1, 1)
  ) + ggplot2::geom_point(
    data = dga %>% group_by(
      Immigration, Extinction, Space, Dispersal
    ) %>% summarise(
      `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
      `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
      .groups = "drop"
    ),
    shape = 4, size = 3, color = "red", alpha = 0.8
  ) + ggplot2::scale_fill_viridis_c(
    direction = -1, trans = "log10"
  ) + ggplot2::facet_grid(
    Immigration + Extinction ~ Space + Dispersal,
    labeller = ggplot2::label_parsed
  ) + ggplot2::theme_bw(
  ) + ggplot2::geom_abline(
    slope = 1, intercept = 0
  ) + ggplot2::labs(
    title = if (exists("dTraj")) {
      "Sample trajectories in diversity counts"
    } else {
      "Gamma by Alpha"
    },
    subtitle = paste0("Pool Modifier = ", pbs, ", ",
                      "Noise Modifier = ", inm)
  )

  if (exists("dTraj")) {
    obj <- obj + ggplot2::geom_path(
      data = dTraj,
      mapping = ggplot2::aes(
        group = interaction(Environment, Number, History),
        color = Time
      ), alpha = 0.3
    ) + ggplot2::geom_point(
      data = dTrajStop,
      shape = 1, size = 3, color = "yellow", alpha = 0.8
    ) + ggplot2::scale_color_viridis_c(
      option = "B"
    )
  }

  return(obj)
}

