# Subset Version 2 to only apply to
#  1. Increased Immigration, Xor Extinction, or Equal Rates and
#  2. Spaces of 10^9, 10^5, or 10^0
# In principle, one should load the MasterParameters.RData,
# but there was the problem in my code, see 31/01/2022 1400 Calendar Note.
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)

clust <- parallel::makeCluster(4, outfile = "")
doParallel::registerDoParallel(clust)

time_grouping_size <- 100
time_averaging_size <- 10

files <- dir(#"Viking_
             "SaveDiversity_2022-02-15",
             #_2021-12-28_2022-01-23",
             pattern = "[.]RData$", full.names = TRUE)

# Diversities <- list()

Diversities <- foreach(
  file = iterators::iter(files),
  .packages = c("dplyr", "ggplot2")
) %dopar% {
  # Retrieve the trailing id numbers from before the file extension.
  # 1st: Set, 2nd: CaseNumber, 3rd: History, 4th: Part
  # Note, if last two are not present, all histories are bundled together.
  # If last two are present, each file has a single part of a single history.
  idNums <- suppressWarnings(na.omit(
    as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                             split = "-", fixed = TRUE)[[1]],
                    n = 4))),
    classes = "simpleWarning")

  print(file)
  print(Sys.time())
  load(file)
  # Loads Diversity

  Attributes <-
    data.frame(
      Pool = toString(Diversity$PoolMod),
      Noise = toString(Diversity$NoiseMod),
      Neutral = toString(Diversity$NeutralMod),
      Space = toString(Diversity$SpaceMod),
      Set = idNums[1],
      CaseNumber = idNums[2],
      History = if(is.na(idNums[3])) {
        # All Histories and Parts
        # rep(1:10, each = 2)
        1:10
      } else idNums[3],
      Part = if(is.na(idNums[4])) {
        # rep(1:2, times = 10)
        1
      } else idNums[4],
      stringsAsFactors = FALSE
    )

  DiversitiesAlpha <-
    dplyr::bind_rows(
      lapply(
        seq_along(Diversity$Diversities),
        function(i, d, a) d[[i]]$alpha %>% dplyr::mutate(
          Set = a$Set[i],
          Number = a$CaseNumber[i],
          History = a$History[i],
          Pool = a$Pool[i],
          Noise = a$Noise[i],
          Neutral = a$Neutral[i],
          Space = a$Space[i]
        ),
        d = Diversity$Diversities,
        a = Attributes
      )
    ) %>% dplyr::select(-Species) %>% dplyr::mutate(
      Time = floor(Time * time_grouping_size)/time_grouping_size
    ) %>% dplyr::group_by(
      Time, Environment, Set, Number, History
    ) %>% dplyr::summarise(
      Richness = floor(median(Richness)),
      .groups = "drop"
    )

  DiversitiesBeta <-
    dplyr::bind_rows(
      lapply(
        seq_along(Diversity$Diversities),
        function(i, d, a) dplyr::bind_rows(
          d[[i]]$beta
        ) %>% dplyr::mutate(
          Set = a$Set[i],
          Number = a$CaseNumber[i],
          History = a$History[i],
          Pool = a$Pool[i],
          Noise = a$Noise[i],
          Neutral = a$Neutral[i],
          Space = a$Space[i]
        ),
        d = Diversity$Diversities,
        a = Attributes
      )
    ) %>% dplyr::mutate(
      Time = floor(Time * time_grouping_size)/time_grouping_size
    ) %>% dplyr::group_by(
      Time, Env1, Env2, Set, Number, History
    ) %>% dplyr::summarise(
      Jaccard = median(Jaccard),
      .groups = "drop"
    )

  DiversitiesGamma <-
    dplyr::bind_rows(lapply(
      1:10,
      function(i, d) d[[i]]$gamma %>% dplyr::mutate(
        History = i, Case = case, Number = caseNumber
      ),
      d = Diversity)) %>% dplyr::filter(
        Aggregation == "Gamma"
      ) %>% dplyr::mutate(
        Time = floor(Time * time_grouping_size)/time_grouping_size
      ) %>% dplyr::group_by(
        Time, History, Case
      ) %>% summarise(
        Richness = floor(median(Richness)),
        .groups = "drop"
      ) %>% dplyr::mutate(
        Pool = toString(Diversity$PoolMod),
        Noise = toString(Diversity$NoiseMod),
        Neutral = toString(Diversity$NeutralMod),
        Space = toString(Diversity$SpaceMod)
      )

  return(list(
    Attr = Attributes,
    Alpha = DiversitiesAlpha,
    Beta = DiversitiesBeta,
    Gamma = DiversitiesGamma
  ))
}

parallel::stopCluster(clust)

Attributes <- dplyr::bind_rows(lapply(
  Diversities, function(d) {d$Attr}
))
DiversitiesAlpha <- dplyr::bind_rows(lapply(
  Diversities, function(d) {d$Alpha}
))
DiversitiesBeta <- dplyr::bind_rows(lapply(
  Diversities, function(d) {d$Beta}
))
DiversitiesGamma <- dplyr::bind_rows(lapply(
  Diversities, function(d) {d$Gamma}
))

print("Pool")
print(table(Attributes$Pool))
# print(table(unlist(lapply(Diversities, function(x) toString(x$PoolMod)))))
print("Noise")
print(table(Attributes$Noise))
# print(table(unlist(lapply(Diversities, function(x) toString(x$NoiseMod)))))
print("Neutral")
print(table(Attributes$Neutral))
# print(table(unlist(lapply(Diversities, function(x) toString(x$NeutralMod)))))
print("Space")
print(table(Attributes$Space))
# print(table(unlist(lapply(Diversities, function(x) toString(x$SpaceMod)))))

print("Preprocessing")

rm(Diversities)

DiversitiesAlphaGamma <- dplyr::full_join(
  DiversitiesAlpha,
  DiversitiesGamma,
  by = c("Time", "History", "Case",
         "Pool", "Noise", "Neutral", "Space"),
  suffix = c(", Alpha", ", Gamma")
)

DiversitiesAlphaAvg <- DiversitiesAlpha %>% dplyr::mutate(
  Time = floor(Time * time_averaging_size)/time_averaging_size
) %>%  dplyr::group_by(
  Time, Pool, Noise, Neutral, Space
) %>% dplyr::summarise(
  RichnessLow = quantile(Richness, p = 0.1),
  RichnessHi = quantile(Richness, p = 0.9),
  Richness = mean(Richness),
  History = "Mean",
  Case = "Mean",
  Environment = "Mean",
  .groups = "drop"
)

DiversitiesBetaAvg <- DiversitiesBeta %>% dplyr::mutate(
  Time = floor(Time * time_averaging_size)/time_averaging_size
) %>%  dplyr::group_by(
  Time, Pool, Noise, Neutral, Space
) %>% dplyr::summarise(
  JaccardLow = quantile(Jaccard, p = 0.1, na.rm = TRUE),
  JaccardHi = quantile(Jaccard, p = 0.9, na.rm = TRUE),
  Jaccard = mean(Jaccard, na.rm = TRUE),
  History = "Mean",
  Case = "Mean",
  Env1 = "Mean",
  Env2 = "Mean",
  .groups = "drop"
)

DiversitiesGammaAvg <- DiversitiesGamma %>% dplyr::mutate(
  Time = floor(Time * time_averaging_size)/time_averaging_size
) %>%  dplyr::group_by(
  Time, Pool, Noise, Neutral, Space
) %>% dplyr::summarise(
  RichnessLow = quantile(Richness, p = 0.1),
  RichnessHi = quantile(Richness, p = 0.9),
  Richness = mean(Richness),
  History = "Mean",
  Case = "Mean",
  .groups = "drop"
)

levelsPoolNoise <- paste(
  DiversitiesGamma$Pool,
  DiversitiesGamma$Noise,
  sep = ";;"
) %>% unique() %>% strsplit(
  split = ";;", fixed = TRUE
)

print("Plotting Alpha")

plotsAlpha <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    da <- DiversitiesAlpha %>% dplyr::filter(
      pbs == Pool,
      inm == Noise
    )

    daa <- DiversitiesAlphaAvg %>% dplyr::filter(
      pbs == Pool,
      inm == Noise
    )

    obj <- ggplot2::ggplot(
      da,
      ggplot2::aes(x = Time, y = Richness,
                   group = interaction(History, Case, Environment))
    ) + ggplot2::geom_line(
      alpha = 0.2
    ) + ggplot2::geom_line(
      data = daa,
      color = "red"
    ) + ggplot2::geom_ribbon(
      data = daa,
      mapping = ggplot2::aes(
        ymin = RichnessLow,
        ymax = RichnessHi,
        x = Time
      ),
      alpha = 0.4,
      fill = "blue",
      inherit.aes = FALSE
    ) + ggplot2::theme_bw(
    ) + ggplot2::facet_grid(
      Neutral ~ Space
    ) + ggplot2::labs(
      title = "Alpha Diversities (Black), Mean (red), and Inner 80% (Blue)",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )

    return(obj)
  }
)

print("Plotting Beta")

plotsBeta <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    db <- DiversitiesBeta %>% dplyr::filter(
      pbs == Pool,
      inm == Noise, History == 1
    )

    dba <- DiversitiesBetaAvg %>% dplyr::filter(
      pbs == Pool,
      inm == Noise
    )

    obj <- ggplot2::ggplot(
      db,
      ggplot2::aes(x = Time, y = Jaccard,
                   group = interaction(History, Case, Env1, Env2))
    ) + ggplot2::geom_line(
      alpha = 0.2
    ) + ggplot2::geom_line(
      data = dba,
      color = "red"
    ) + ggplot2::geom_ribbon(
      data = dba,
      mapping = ggplot2::aes(
        ymin = JaccardLow,
        ymax = JaccardHi,
        x = Time
      ),
      alpha = 0.4,
      fill = "blue",
      inherit.aes = FALSE
    ) + ggplot2::theme_bw(
    ) + ggplot2::facet_grid(
      Neutral ~ Space
    ) + ggplot2::labs(
      title = "Beta Diversities (Black), Mean (red), and Inner 80% (Blue)",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm),
      y = "Jaccard Distance"
    )

    return(obj)
  }
)

print("Plotting Gamma")

plotsGamma <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    dg <- DiversitiesGamma %>% dplyr::filter(
      pbs == Pool,
      inm == Noise
    )

    dga <- DiversitiesGammaAvg %>% dplyr::filter(
      pbs == Pool,
      inm == Noise
    )

    obj <- ggplot2::ggplot(
      dg,
      ggplot2::aes(x = Time, y = Richness,
                   group = interaction(History, Case))
    ) + ggplot2::geom_line(
      alpha = 0.2
    ) + ggplot2::geom_line(
      data = dga,
      color = "red"
    ) + ggplot2::geom_ribbon(
      data = dga,
      mapping = ggplot2::aes(
        ymin = RichnessLow,
        ymax = RichnessHi,
        x = Time
      ),
      alpha = 0.4,
      fill = "blue",
      inherit.aes = FALSE
    ) + ggplot2::theme_bw(
    ) + ggplot2::facet_grid(
      Neutral ~ Space
    ) + ggplot2::labs(
      title = "Gamma Diversities (Black), Mean (red), and Inner 80% (Blue)",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )

    return(obj)
  }
)

print("Plotting Alpha + Gamma")

plotsAlphaGamma <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    dag <- DiversitiesAlphaGamma %>% dplyr::filter(
      Pool == pbs,
      Noise == inm,
      Time > 3
    )

    ggplot2::ggplot(
      dag,
      ggplot2::aes(x = `Richness, Alpha`,
                   y = `Richness, Gamma`)
    ) + ggplot2::geom_bin2d(
      breaks = seq(-0.5, max(dag$`Richness, Alpha`,
                             dag$`Richness, Gamma`) + 1, 1)
    ) + ggplot2::geom_point(
      data = dag %>% group_by(
        Neutral, Space
      ) %>% summarise(
        `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
        `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
        .groups = "drop"
      ),
      shape = 4, size = 3, color = "red", alpha = 0.8
    ) + ggplot2::scale_fill_viridis_c(
      direction = -1, trans = "log10"
    ) + ggplot2::facet_grid(
      Neutral ~ Space
    ) + ggplot2::theme_bw(
    ) + ggplot2::geom_abline(
      slope = 1, intercept = 0
    ) + ggplot2::labs(
      title = "Gamma by Alpha",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )
  }
)

plotsAlphaGamma3x3 <- lapply(
  1, function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    dag <- DiversitiesAlphaGamma %>% dplyr::filter(
      Pool == pbs,
      Noise == inm,
      Time > 3,
      Space %in% c(1, 1E5, 1E6, 1E9),
      Neutral %in% c(#"1, 1", Currently going to borrow from No Noise.
        "1, 10", "10, 1")
    )

    dag <- dplyr::bind_rows(
      dag,
      DiversitiesAlphaGamma %>% dplyr::filter(
        Pool == pbs,
        Noise == "0",
        Time > 3,
        Space %in% c(1, 1E5, 1E6, 1E9),
        Neutral %in% c("1, 1")
      )
    )

    dag <- dag %>% dplyr::mutate(
      Neutral = dplyr::case_when(
        Neutral == "1, 10" ~ "High Ext.",
        Neutral == "1, 1" ~ "Equal Im. & Ext.",
        Neutral == "10, 1" ~ "High Im.",
        TRUE ~ "Oops"
      ),
      Neutral = factor(
        Neutral, levels = c("High Im.", "Equal Im. & Ext.", "High Ext.")
      )
    )

    ggplot2::ggplot(
      dag,
      ggplot2::aes(x = `Richness, Alpha`,
                   y = `Richness, Gamma`)
    ) + ggplot2::geom_bin2d(
      breaks = seq(-0.5, max(dag$`Richness, Alpha`,
                             dag$`Richness, Gamma`) + 1, 1)
    ) + ggplot2::geom_point(
      data = dag %>% group_by(
        Neutral, Space
      ) %>% summarise(
        `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
        `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
        .groups = "drop"
      ),
      shape = 4, size = 3, color = "red", alpha = 0.8
    ) + ggplot2::scale_fill_viridis_c(
      direction = -1, trans = "log10"
    ) + ggplot2::facet_grid(
      Neutral ~ Space
    ) + ggplot2::theme_bw(
    ) + ggplot2::geom_abline(
      slope = 1, intercept = 0
    ) + ggplot2::labs(
      title = "Gamma by Alpha"#,
      # subtitle = paste0("Pool Modifier = ", pbs, ", ",
      #                   "Noise Modifier = ", inm)
    )
  }
)

print("Plotting Alpha + Gamma + Space")

plotsAlphaGammaSpace <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    dag <- DiversitiesAlphaGamma %>% dplyr::filter(
      Pool == pbs,
      Noise == inm,
      Time > 3
    )

    ggplot2::ggplot(
      dag,
      ggplot2::aes(x = `Richness, Alpha`,
                   y = `Richness, Gamma`)
    ) + ggplot2::geom_bin2d(
      breaks = seq(-0.5, max(dag$`Richness, Alpha`,
                             dag$`Richness, Gamma`) + 1, 1),
      alpha = 0.2
    ) + ggplot2::geom_point(
      data = dag %>% dplyr::group_by(
        Neutral, Space
      ) %>% dplyr::summarise(
        `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
        `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
        .groups = "drop"
      ),
      mapping = ggplot2::aes(
        color = Neutral,
        shape = Space
      ),
      size = 4, alpha = 1
      # ) + ggplot2::geom_rect(
      #   data = dag %>% dplyr::group_by(
      #     Neutral, Space
      #   ) %>% dplyr::summarise(
      #     xmin = min(`Richness, Alpha`, na.rm = TRUE),
      #     xmax = max(`Richness, Alpha`, na.rm = TRUE),
      #     ymin = min(`Richness, Gamma`, na.rm = TRUE),
      #     ymax = max(`Richness, Gamma`, na.rm = TRUE),
      #     .groups = "drop"
      #   ),
      #   mapping = ggplot2::aes(
      #     xmin = xmin - 0.5, xmax = xmax + 0.5,
      #     ymin = ymin - 0.5, ymax = ymax + 0.5
      #   ), inherit.aes = FALSE, color = "black", fill = "transparent"
    ) + ggplot2::geom_path(
      data = dag %>% dplyr::group_by(
        Neutral, Space
      ) %>% dplyr::summarise(
        `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
        `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
        .groups = "drop"
      ) %>% dplyr::arrange(Space),
      mapping = ggplot2::aes(color = Neutral)
    ) + ggplot2::scale_fill_viridis_c(
      direction = -1, trans = "log10"
    ) + ggplot2::theme_bw(
    ) + ggplot2::geom_abline(
      slope = 1, intercept = 0
    ) + ggplot2::labs(
      title = "Alpha and Gamma By Space, Means",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )
  }
)

print("Plotting Alpha + Space")

plotsAlphaSpace <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    da <- DiversitiesAlpha %>% dplyr::filter(
      pbs == Pool,
      inm == Noise
    )

    obj <- ggplot2::ggplot(
      da,
      ggplot2::aes(x = Space, y = Richness,
                   fill = Neutral)
    ) + ggplot2::geom_violin(
    ) + ggplot2::theme_bw(
    ) + ggplot2::facet_grid(
      Neutral ~ .
    ) + ggplot2::labs(
      title = "Alpha Diversities by Island Distance and Pool Relationship",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )

    return(obj)
  }
)

plotsAlphaSpace2 <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    da <- DiversitiesAlpha %>% dplyr::filter(
      pbs == Pool,
      inm == Noise,
      Time > 3
    )

    obj <- ggplot2::ggplot(
      da,
      ggplot2::aes(x = Space, y = Richness,
                   fill = Neutral, group = interaction(Neutral, Space))
    ) + ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(0.5)
    ) + ggplot2::stat_summary(
      fun = mean,
      position = ggplot2::position_dodge(0.5),
      shape = 4
    ) + ggplot2::labs(
      title = "Alpha Diversities by Island Distances",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )

    return(obj)
  }
)

print("Plotting Gamma + Space")


plotsGammaSpace <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    dg <- DiversitiesAlpha %>% dplyr::filter(
      pbs == Pool,
      inm == Noise
    )

    obj <- ggplot2::ggplot(
      dg,
      ggplot2::aes(x = Space, y = Richness,
                   fill = Neutral)
    ) + ggplot2::geom_violin(
    ) + ggplot2::theme_bw(
    ) + ggplot2::facet_grid(
      Neutral ~ Space
    ) + ggplot2::labs(
      title = "Alpha Diversities by Island Distance and Pool Relationship",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )

    return(obj)
  }
)

plotsGammaSpace2 <- lapply(
  1:length(levelsPoolNoise), function(i) {

    pbs <- levelsPoolNoise[[i]][1]
    inm <- levelsPoolNoise[[i]][2]

    dg <- DiversitiesGamma %>% dplyr::filter(
      pbs == Pool,
      inm == Noise,
      Time > 3
    )

    obj <- ggplot2::ggplot(
      dg,
      ggplot2::aes(x = Space, y = Richness,
                   fill = Neutral, group = interaction(Neutral, Space))
    ) + ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(0.5)
    ) + ggplot2::stat_summary(
      fun = mean,
      position = ggplot2::position_dodge(0.5),
      shape = 4
    ) + ggplot2::labs(
      title = "Gamma Diversities by Island Distances",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )

    return(obj)
  }
)

plotAlphaGammaTrajectory <- function(
  dag, TimeFilter = 3,
  EnvirSelect = 1, CaseSelect = 1, HistorySelect = 1
) {

  cases <- sort(unique(dag$Case))
  case <- cases[CaseSelect]

  ggplot2::ggplot(
    dag %>% dplyr::filter(
      Time > TimeFilter
    ),
    ggplot2::aes(x = `Richness, Alpha`,
                 y = `Richness, Gamma`)
  ) + ggplot2::geom_bin2d(
    breaks = seq(-0.5, max(dag$`Richness, Alpha`,
                           dag$`Richness, Gamma`) + 1, 1)
  ) + ggplot2::geom_path(
    data = dag %>% dplyr::filter(
      Environment == EnvirSelect,
      Case == case,
      History == HistorySelect
    ),
    mapping = ggplot2::aes(
      color = Time
    ), alpha = 0.3
  ) + ggplot2::geom_point(
    data = dag %>% dplyr::filter(
      Time > TimeFilter
    ) %>% group_by(
      Neutral, Space
    ) %>% summarise(
      `Richness, Alpha` = mean(`Richness, Alpha`, na.rm = TRUE),
      `Richness, Gamma` = mean(`Richness, Gamma`, na.rm = TRUE),
      .groups = "drop"
    ),
    shape = 4, size = 3, color = "red", alpha = 0.8
  ) + ggplot2::scale_color_viridis_c(
    option = "B"
  ) + ggplot2::scale_fill_viridis_c(
    direction = -1, trans = "log10"
  ) + ggplot2::facet_grid(
    Neutral + Noise ~ Space + Pool
  ) + ggplot2::theme_bw(
  ) + ggplot2::geom_abline(
    slope = 1, intercept = 0
  )
}

plotAlphaGammaTrajectoryEG <- plotAlphaGammaTrajectory(
  dag = DiversitiesAlphaGamma %>% dplyr::filter(
    Pool == "0, 0, 0, 0",
    Noise == "1",
    Space == "1e+07",
    Neutral == "1, 1"
  ))

plotAlphaGammaTrajectoryEGBIMODAL <- plotAlphaGammaTrajectory(
  dag = DiversitiesAlphaGamma %>% dplyr::filter(
    Pool == "0, 0, 0, 0",
    Noise == "0",
    Space == "1e+09",
    Neutral == "1, 0"
  ))

plotLCABTrust <- DiversitiesAlphaGamma %>% dplyr::filter(
  Time > 3, Space != "NA"
) %>% tidyr::pivot_longer(
  cols = c("Richness, Gamma", "Richness, Alpha"),
  names_to = "Diversity", values_to = "Richness"
) %>% dplyr::mutate(
  Diversity = substring(Diversity, first = 11),
  Diversity = case_when(
    Diversity == "Alpha" ~ "Local",
    Diversity == "Gamma" ~ "Regional",
    TRUE ~ "Oops"
  )
) %>% ggplot2::ggplot(
  ggplot2::aes(x = Space, y = Richness, fill = Diversity)
) + ggplot2::geom_boxplot(
  position = ggplot2::position_dodge(0.8)
) + ggplot2::stat_summary(
  fun = mean,
  position = ggplot2::position_dodge(0.8),
  shape = 4
) + ggplot2::labs(
  x = "Island-Island Distances",
  y = "Number of Species"
)
