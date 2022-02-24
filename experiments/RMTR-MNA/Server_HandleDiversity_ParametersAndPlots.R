# Due to a typo, the parameters are failing to be saved with each run.
# To compensate, we need to load the required files and translate them again.
library(dplyr)
library(ggplot2)

time_grouping_size <- 100
time_averaging_size <- 1000

# Note alphabetical, using 1:3 <=> A:C, and one letter is the only difference.
filesParameters <- dir(path = ".", pattern = "^MNA[-].+Cases[.]csv$",
                       full.names = TRUE, recursive = TRUE,
                       include.dirs = TRUE, no.. = TRUE)
parametersMaster <- load("MNA-MasterParameters.RData")
parametersCSVs <- lapply(filesParameters, read.csv)

files <- dir("Viking_SaveDiversity_2021-12-28_2022-01-23",
             pattern = "[.]RData$", full.names = TRUE)

Diversities <- list()

for (file in files) {
  load(file)
  # Loads Diversity

  # Retrieve the trailing two id numbers from before the file extension.
  idNums <- as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                                     split = "-", fixed = TRUE)[[1]],
                            n = 2))

  theParameters <- parametersCSVs[[idNums[1]]][(idNums[2] %/% 10) + 1,]
  Diversity$PoolMod <- systemMods$PoolBodySizes[theParameters$Framework,]
  Diversity$NoiseMod <- systemMods$InteractionsNoiseMultiplier[theParameters$Framework]
  Diversity$NeutralMod<- systemMods$NeutralRateMultipliers[theParameters$Neutral,]
  Diversity$SpaceMod <- systemMods$SpaceDistanceMultiplier[theParameters$Space]

  Diversities[[toString(idNums)]] <- Diversity
}

print("Pool")
print(table(unlist(lapply(Diversities, function(x) toString(x$PoolMod)))))
print("Noise")
print(table(unlist(lapply(Diversities, function(x) toString(x$NoiseMod)))))
print("Neutral")
print(table(unlist(lapply(Diversities, function(x) toString(x$NeutralMod)))))
print("Space")
print(table(unlist(lapply(Diversities, function(x) toString(x$SpaceMod)))))


DiversitiesAlpha <- dplyr::bind_rows(lapply(
  Diversities, function(x) {
    temp <- dplyr::bind_rows(lapply(
      1:10,
      function(i, d) d[[i]]$alpha %>% dplyr::mutate(Simulation = i),
      d = x)) %>% dplyr::select(
        -Species
        ) %>% dplyr::mutate(
        Time = floor(Time * time_grouping_size)/time_grouping_size
      ) %>% dplyr::group_by(
        Time, Simulation, Environment
      ) %>% summarise(
        Richness = floor(median(Richness))
      )

    temp <- temp %>% dplyr::mutate(
      Pool = toString(x$PoolMod),
      Noise = toString(x$NoiseMod),
      Neutral = toString(x$NeutralMod),
      Space = toString(x$SpaceMod)
    )

    return(temp)
  }
))

# Note DiversitiesBeta is formed from lists of lists, hence the extra bind_rows.
DiversitiesBeta <- dplyr::bind_rows(lapply(
  Diversities, function(x) {
    temp <- dplyr::bind_rows(lapply(
      1:10,
      function(i, d) dplyr::bind_rows(
        d[[i]]$beta
      ) %>% dplyr::mutate(Simulation = i),
      d = x)) %>% dplyr::mutate(
        Time = floor(Time * time_grouping_size)/time_grouping_size
      ) %>% dplyr::group_by(
        Time, Simulation, Env1, Env2
      ) %>% summarise(
        Jaccard = floor(median(Jaccard))
      )

    temp <- temp %>% dplyr::mutate(
      Pool = toString(x$PoolMod),
      Noise = toString(x$NoiseMod),
      Neutral = toString(x$NeutralMod),
      Space = toString(x$SpaceMod)
    )

    return(temp)
  }
))

DiversitiesGamma <- dplyr::bind_rows(lapply(
  Diversities, function(x) {
    temp <- dplyr::bind_rows(lapply(
      1:10,
      function(i, d) d[[i]]$gamma %>% dplyr::mutate(Simulation = i),
      d = x)) %>% dplyr::mutate(
        Time = floor(Time * time_grouping_size)/time_grouping_size
      ) %>% dplyr::group_by(
        Time, Simulation
      ) %>% summarise(
        Richness = floor(median(Richness))
      )

    temp <- temp %>% dplyr::mutate(
      Pool = toString(x$PoolMod),
      Noise = toString(x$NoiseMod),
      Neutral = toString(x$NeutralMod),
      Space = toString(x$SpaceMod)
    )

    return(temp)
  }
))

print("Preprocessing")

rm(Diversities)
rm(Diversity)

DiversitiesGamma <- DiversitiesGamma %>% dplyr::filter(
  Aggregation == "Gamma"
)

DiversitiesGammaAvg <- DiversitiesGamma %>% dplyr::mutate(
  Time = floor(Time * time_averaging_size)/time_averaging_size
) %>%  dplyr::group_by(
  Time, Pool, Noise, Neutral, Space
) %>% dplyr::summarise(
  RichnessLow = quantile(Richness, p = 0.1),
  RichnessHi = quantile(Richness, p = 0.9),
  Richness = mean(Richness),
  Simulation = "Mean"
)

print("Plotting 1")

plotsGamma <- lapply(
  1:nrow(systemMods$PoolBodySizes), function(i) {

    pbs <- toString(systemMods$PoolBodySizes[i, ])
    inm <- toString(systemMods$InteractionsNoiseMultiplier[i])

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
      ggplot2::aes(x = Time, y = Richness, group = Simulation)
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

