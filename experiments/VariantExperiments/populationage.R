ages <- lapply(seq_along(environs), function(i, e) {
  agesByColumn <- apply(
    e[[i]]$Abundance, MARGIN = 2, function(column) {
      positives <- which(column > result$Parameters$EliminationThreshold)
      Populations <- cumsum(c(FALSE, diff(positives) > 1))
      numPopulationsMinus1 <- Populations[length(Populations)]
      if (length(numPopulationsMinus1) == 0) {
        NULL
      } else if (numPopulationsMinus1 == 0) {
        list(data.frame(
          time = times[positives],
          age = times[positives] - times[positives[1]],
          occurrence = 0,
          patch = i
        ))
      } else {
        lapply(0:(numPopulationsMinus1), function(j) {
          subset <- positives[Populations == j]
          data.frame(
            time = times[subset],
            age = times[subset] - times[subset[1]],
            occurrence = j,
            patch = i
          )
        })
      }
    })
  
  lapply(seq_along(agesByColumn), function(i, aBC, name) {
    lapply(aBC[[i]], function(a) dplyr::bind_cols(a, species = name[i]))
  }, aBC = agesByColumn, name = names(agesByColumn)) %>% dplyr::bind_rows()
},
e = environs) %>% dplyr::bind_rows()

ggplot2::ggplot(
  ages,
  ggplot2::aes(x = time, y = age, 
               group = interaction(occurrence, species, patch),
               color = factor(patch), linetype = factor(patch))
) + ggplot2::geom_line(
  # ) + ggplot2::geom_point(
  #   data = ages %>% dplyr::group_by(patch, time) %>% dplyr::summarise(
  #     age = mean(age),
  #     .groups = "drop"
  #   ), 
  #   mapping = ggplot2::aes(x = time, y = age, color = factor(patch)),
  #   shape = '.',
  #   inherit.aes = FALSE
) + ggplot2::geom_ribbon(
  data = ages %>% dplyr::group_by(patch, time) %>% dplyr::summarise(
    lowage = quantile(age, p = 0.25),
    highage = quantile(age, p = 1 - 0.25),
    .groups = "drop"
  ),
  mapping = ggplot2::aes(
    x = time, color = factor(patch),
    ymin = lowage, ymax = highage
  ),
  inherit.aes = FALSE, alpha = 0.1
)
