load("diversitiesFlattened9a8_plottable.RData")
height <- 3600; width <- 4800; units <- "px"

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("rep_0"),
    NicheDistance == 7, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: 0; Time > 30,000; Effect Magnitude: 10"
)

ggplot2::ggsave("image0_magnitude10_2024-12-16.png",
                height = height, width = width, units = units)

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("rep_1"),
    NicheDistance == 7, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: 1; Time > 30,000; Effect Magnitude: 10"
)

ggplot2::ggsave("image1_magnitude10_2024-12-16.png",
                height = height, width = width, units = units)

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("evensplit_01"),
    NicheDistance == 7, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: Even Split 0+1; Time > 30,000; Effect Magnitude: 10"
)

ggplot2::ggsave("image01_magnitude10_2024-12-16.png",
                height = height, width = width, units = units)

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("runif"),
    NicheDistance == 7, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: Unif.@Rand.; Time > 30,000; Effect Magnitude: 10"
)

ggplot2::ggsave("imagerunif_magnitude10_2024-12-16.png",
                height = height, width = width, units = units)

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("rep_0"),
    NicheDistance == 3, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: 0; Time > 30,000; Effect Magnitude: 2"
)

ggplot2::ggsave("image0_magnitude2_2024-12-16.png",
                height = height, width = width, units = units)

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("rep_1"),
    NicheDistance == 3, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: 1; Time > 30,000; Effect Magnitude: 2"
)

ggplot2::ggsave("image1_magnitude2_2024-12-16.png",
                height = height, width = width, units = units)

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("evensplit_01"),
    NicheDistance == 3, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: Even Split 0+1; Time > 30,000; Effect Magnitude: 2"
)

ggplot2::ggsave("image01_magnitude2_2024-12-16.png",
                height = height, width = width, units = units)

plotDiversityOverview(
  diversitiesFlattened %>% dplyr::filter(
    SpeciesAffinity %in% c("runif"),
    NicheDistance == 3, Time > 30000,
    is.na(Subset)
  ) %>% dplyr::mutate(
    Intervention = factor(
      Intervention,
      levels = c(
        "(0)", "(0)->(0)", "(0)->(0.25)", "(0)->(0.5)", "(0)->(0.75)", "(0)->(1)",
        "(0.25)", "(0.25)->(0)", "(0.25)->(0.25)", "(0.25)->(0.5)", "(0.25)->(0.75)", "(0.25)->(1)",
        "(0.5)", "(0.5)->(0)", "(0.5)->(0.25)", "(0.5)->(0.5)", "(0.5)->(0.75)", "(0.5)->(1)",
        "(0.75)", "(0.75)->(0)", "(0.75)->(0.25)", "(0.75)->(0.5)", "(0.75)->(0.75)", "(0.75)->(1)",
        "(1)", "(1)->(0)", "(1)->(0.25)", "(1)->(0.5)", "(1)->(0.75)", "(1)->(1)"
      )
    )),
  "Alpha Hill:0"
) + ggplot2::facet_wrap(
  .~Intervention,
  nrow = 5
) + ggplot2::theme(
  legend.position = "none"
) + ggplot2::geom_boxplot(
  ggplot2::aes(x = 34000),
  width = 1000
) + ggplot2::labs(
  title = "Pool: Unif.@Rand.; Time > 30,000; Effect Magnitude: 2"
)

ggplot2::ggsave("imagerunif_magnitude2_2024-12-16.png",
                height = height, width = width, units = units)