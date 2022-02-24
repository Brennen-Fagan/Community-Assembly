# Copied from Server_HandleDiversity_Paremeters_AndPlots2Subset.R
# The idea is to create N workers (N cpus-per-task) to parallelise loading the
# diversity data, then analyse and generate the plots.
# We'll see how well it works, but on the server it seemed like parallelising
# the load of the data was dealing with th bottleneck.
# Total memory at this point is North of 150 GB, so we should try 250.
# Note that Time has already been divided by 1E4.

# OLD ABSTRACT: ###############################################################
# # Subset Version 2 to only apply to
# #  1. Increased Immigration, Xor Extinction, or Equal Rates and
# #  2. Spaces of 10^9, 10^5, or 10^0
# # In principle, one should load the MasterParameters.RData,
# # but there was the problem in my code, see 31/01/2022 1400 Calendar Note.

# Setup: ######################################################################
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)

directory <- '.'
outputLocation <- file.path(directory, paste0("SaveImages_", Sys.Date()))
if (!dir.exists(outputLocation)) {
  dir.create(outputLocation, showWarnings = FALSE)
}

cargs <- as.numeric(commandArgs(TRUE))

clust <- parallel::makeCluster(cargs[1], outfile = "")
doParallel::registerDoParallel(clust)

time_grouping_size <- 100
time_averaging_size <- 10
time_burnin <- 3
time_units <- 1E4

# csvs <- dir(path = ".", pattern = "Cases[.]csv$",
#             full.names = TRUE)
# 
# contentsCsvs <- lapply(csvs, utils::read.csv)
# 
# # From old MasterParameters.RData
# # 10^9, 10^5, 10^0
# spaces <- c(1, 5, 8)
# # (1,1), (1,10), and (10,1)
# neutrals <- c(1, 3, 6)
# 
# IDs <- lapply(seq_along(contentsCsvs), function(i) {
#   csv <- contentsCsvs[[i]]
#   temp <- data.frame(Numbers = csv %>% dplyr::mutate(
#     Numbers = as.numeric(rownames(csv))
#   ) %>% dplyr::filter(
#     Space %in% spaces,
#     Neutral %in% neutrals
#   ) %>% dplyr::pull(Numbers)
#   )
#   
#   if (nrow(temp) > 0) {
#     temp$Case <- i
#     return(temp)
#   }
#   else {
#     return(NULL)
#   }
# }) %>% dplyr::bind_rows()

files <- dir("SaveDiversity_2022-01-25",
             pattern = "[.]RData$", full.names = TRUE)

# Diversities <- list()

Diversities <- foreach(
  file = iterators::iter(files),
  .packages = c("dplyr", "ggplot2")
) %dopar% {
  #for (file in files) {
  # Retrieve File Number == Case Number
  caseNumber <- strsplit(tools::file_path_sans_ext(file),
                         split = "-", fixed = TRUE)[[1]]
  case <- caseNumber[length(caseNumber) - 1]
  caseNumber <- caseNumber[length(caseNumber)]
  
  # if (IDs %>% dplyr::filter(
  #   Case == as.numeric(case),
  #   (as.numeric(caseNumber) %/% 10) + 1 == Numbers
  # ) %>% nrow() == 0) {
  #   return(NULL)
  # }
  
  print(file)
  print(Sys.time())
  load(file)
  # Loads Diversity
  
  # Diversities[[toString(idNums)]] <- Diversity
  Attributes <- # dplyr::bind_rows(
    # if (exists("Attributes")) {
    #   Attributes
    # },
    data.frame(
      Pool = toString(Diversity$PoolMod),
      Noise = toString(Diversity$NoiseMod),
      Neutral = toString(Diversity$NeutralMod),
      Space = toString(Diversity$SpaceMod),
      Case = caseNumber,
      stringsAsFactors = FALSE
    )
  # )
  
  DiversitiesAlpha <- #dplyr::bind_rows(
    # if (exists("DiversitiesAlpha")) {
    #   DiversitiesAlpha
    # },
    dplyr::bind_rows(lapply(
      1:10,
      function(i, d) d[[i]]$alpha %>% dplyr::mutate(
        History = i, Case = caseNumber
      ),
      d = Diversity)) %>% dplyr::select(-Species) %>% dplyr::mutate(
        Time = floor(Time * time_grouping_size)/time_grouping_size
      ) %>% dplyr::group_by(
        Time, History, Environment, Case
      ) %>% summarise(
        Richness = floor(median(Richness)),
        .groups = "drop"
      ) %>% dplyr::mutate(
        Pool = toString(Diversity$PoolMod),
        Noise = toString(Diversity$NoiseMod),
        Neutral = toString(Diversity$NeutralMod),
        Space = toString(Diversity$SpaceMod)
      )
  # )
  
  DiversitiesBeta <- # dplyr::bind_rows(
    # if (exists("DiversitiesBeta")) {
    #   DiversitiesBeta
    # },
    dplyr::bind_rows(lapply(
      1:10,
      function(i, d) dplyr::bind_rows(
        d[[i]]$beta
      ) %>% dplyr::mutate(History = i, Case = caseNumber),
      d = Diversity)) %>% dplyr::mutate(
        Time = floor(Time * time_grouping_size)/time_grouping_size
      ) %>% dplyr::group_by(
        Time, History, Env1, Env2, Case
      ) %>% summarise(
        Jaccard = median(Jaccard),
        .groups = "drop"
      ) %>% dplyr::mutate(
        Pool = toString(Diversity$PoolMod),
        Noise = toString(Diversity$NoiseMod),
        Neutral = toString(Diversity$NeutralMod),
        Space = toString(Diversity$SpaceMod)
      )
  # )
  
  DiversitiesGamma <- # dplyr::bind_rows(
    # if (exists("DiversitiesGamma")) {
    #   DiversitiesGamma
    # },
    dplyr::bind_rows(lapply(
      1:10,
      function(i, d) d[[i]]$gamma %>% dplyr::mutate(
        History = i, Case = caseNumber
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
  # )
  
  return(list(
    Attr = Attributes,
    Alpha = DiversitiesAlpha,
    Beta = DiversitiesBeta,
    Gamma = DiversitiesGamma
  ))
}

parallel::stopCluster(clust)

Attributes <- dplyr::bind_rows(parLapply(
  Diversities, function(d) {d$Attr}
))
DiversitiesAlpha <- dplyr::bind_rows(parLapply(
  Diversities, function(d) {d$Alpha}
))
DiversitiesBeta <- dplyr::bind_rows(parLapply(
  Diversities, function(d) {d$Beta}
))
DiversitiesGamma <- dplyr::bind_rows(parLapply(
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

levelsPoolNoise <- interaction(
  DiversitiesGamma$Pool,
  DiversitiesGamma$Noise,
  sep = ";;"
) %>% levels() %>% strsplit(
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
                        "Noise Modifier = ", inm),
      x = bquote(paste(Time, ",") ~ 10^.(log10(time_units)))
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
      y = "Jaccard Distance",
      x = bquote(paste(Time, ",") ~ 10^.(log10(time_units)))
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
                        "Noise Modifier = ", inm),
      x = bquote(paste(Time, ",") ~ 10^.(log10(time_units)))
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
      Time > time_burnin
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
      Time > time_burnin
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
    ) + ggplot2::geom_rect(
      data = dag %>% dplyr::group_by(
        Neutral, Space
      ) %>% dplyr::summarise(
        xmin = min(`Richness, Alpha`, na.rm = TRUE),
        xmax = max(`Richness, Alpha`, na.rm = TRUE),
        ymin = min(`Richness, Gamma`, na.rm = TRUE),
        ymax = max(`Richness, Gamma`, na.rm = TRUE),
        .groups = "drop"
      ),
      mapping = ggplot2::aes(
        xmin = xmin - 0.5, xmax = xmax + 0.5,
        ymin = ymin - 0.5, ymax = ymax + 0.5
      ), inherit.aes = FALSE, color = "black", fill = "transparent"
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
      Time > time_burnin
    )
    
    obj <- ggplot2::ggplot(
      da,
      ggplot2::aes(x = Space, y = Richness,
                   fill = Neutral, group = interaction(Neutral, Space))
    ) + ggplot2::geom_boxplot(
      notch = TRUE
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
      Time > time_burnin
    )
    
    obj <- ggplot2::ggplot(
      dg,
      ggplot2::aes(x = Space, y = Richness,
                   fill = Neutral, group = interaction(Neutral, Space))
    ) + ggplot2::geom_boxplot(
      notch = TRUE
    ) + ggplot2::labs(
      title = "Gamma Diversities by Island Distances",
      subtitle = paste0("Pool Modifier = ", pbs, ", ",
                        "Noise Modifier = ", inm)
    )
    
    return(obj)
  }
)

plotsList <- list(
  plotsAlpha,
  plotsBeta,
  plotsGamma,
  plotsAlphaGamma,
  plotsAlphaGammaSpace,
  plotsAlphaSpace,
  plotsAlphaSpace2,
  plotsGammaSpace,
  plotsGammaSpace2
)

for (p in seq_along(plotsList)) {
  for (pp in seq_along(plotsList[[p]])) {
    plt <- plotsList[[p]][[pp]]
    ggplot2::ggsave(
      plot = plt,
      filename = file.path(
        outputLocation,
        paste0("DiversityPlot-", p, "-", pp, "-.png")
      ),
      units = "cm",
      width = 16, height = 12
    )
  }
}

# Possibility that this might work:
save(
  time_grouping_size,
  time_averaging_size,
  time_burnin,
  time_units,
  Attributes,
  DiversitiesAlpha,
  DiversitiesBeta,
  DiversitiesGamma,
  file = "PlottingData_2022-02-01.RData"
)
