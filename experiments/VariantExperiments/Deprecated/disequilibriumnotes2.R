# Version 1 didn't work. Consider a 0 pred, a 0 prey, a 1 pred, a 1 prey in a
# 0 environment. Then the 0 pred and 0 prey take up most of the interactions.
# The 1-1 interaction is nonexistent. The 0-1 and 1-0 interactions are in
# between, and, because of the abundance scaling, are precisely in between
# such that (0,0) * (1,1) ~ (0,1) * (1,0) : Disequilibrium ~ 0.


# Sketch of idea:
# interactions <- matrix(c(-1, -1, 0, 1, -1, -1, 0, 1, -1), byrow = T, nrow = 3)
# entriesPredPrey <- (interactions < 0) * t((interactions > 0))
# currentAbundance <- 1:3
# poolAffinities <- c(1,0,1)
# consumption <- entriesPredPrey * abs(interactions) *
#   currentAbundance[row(interactions)] * currentAbundance[col(interactions)]
#
# factory <- function(consumptionMat, InteractionsMat) {
#   norm <- sum(consumptionMat)
#   ro <- row(InteractionsMat)
#   co <- col(InteractionsMat)
#   function(i, j) {
#     sum(consumptionMat *
#           ifelse(poolAffinities == i, 1, 0)[ro] *
#           ifelse(poolAffinities == j, 1, 0)[co]) / norm
#   }
# }
#
# product <- factory(consumption, interactions)
# product <- Vectorize(product)
#
# det(outer(c(0,1), c(0, 1), FUN = product))

affinities <- unique(unlist(lapply(environs, function(e) {e$Affinity})))
affinities <- sort(affinities)
diseq <- lapply(1:nrow(environs[[1]]$Abundance), function(timestep) {
  disequilibriums <- lapply(environs, function(e, ts) {
    interactions <- e$Matrix
    currentAbundance <- e$Abundance[ts, ]
    poolAffinities <- e$Affinity

    entriesPredPrey <- (interactions < 0) * t((interactions > 0))
    consumption <- entriesPredPrey * abs(interactions) *
      currentAbundance[row(interactions)] * currentAbundance[col(interactions)]


    factory <- function(consumptionMat, InteractionsMat) {
      norm <- sum(consumptionMat)
      ro <- row(InteractionsMat)
      co <- col(InteractionsMat)
      function(i, j) {
        sum(consumptionMat *
              ifelse(poolAffinities == i, 1, 0)[ro] *
              ifelse(poolAffinities == j, 1, 0)[co]) / norm
      }
    }

    product <- factory(consumption, interactions)
    product <- Vectorize(product)

    # det(outer(affinities, affinities, FUN = product))
    # sum(unlist(lapply(affinities, function(i) product(i,i))))
    values <- outer(affinities, affinities, FUN = product)
    sum(diag(values)) / (sum(values) - sum(diag(values)))

  }, ts = timestep)

  data.frame(
    time = times[timestep],
    environment = 1:length(environs),
    disequilibrium = unlist(disequilibriums)
  )
  }) %>% dplyr::bind_rows()

ggplot2::ggplot(
  diseq,
  ggplot2::aes(x = time, y = disequilibrium, color = factor(environment))
) + ggplot2::geom_point(
  shape = '.'
  ) + ggplot2::scale_y_log10(
  ) + ggplot2::geom_hline(yintercept = 1)
