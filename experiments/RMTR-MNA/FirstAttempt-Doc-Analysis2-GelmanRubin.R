Richness_Binned <- Diversity_AG_Binned %>% dplyr::filter(
  Environment != "Mean",
  TimeFloor > 2
) %>% dplyr::mutate(
  Time = TimeFloor * divide_time_by
) %>% dplyr::group_by(
  Time, Environment, 
  Modifier, ModIntensity, Space, Distance
) %>% dplyr::summarise(
  # Why the complicated expression and not mean or median?
  # I want it to be a state that we know we can reach,
  # and so I am taking the closest state to the median.
  # The tail makes sure I am only looking at distances 
  # to the median, as opposed to all distances.
  Value = Value[which.min(
    tail(dist(matrix(c(Value, median(Value)))), length(Value))
    )]
) %>% dplyr::ungroup(
) %>% dplyr::group_split(
  Modifier, ModIntensity, Space, Distance
)

RichnessMC <- Richness_Binned %>% lapply(
  function (rich) {
    coda::mcmc.list(
      lapply(rich %>% dplyr::group_split(Environment),
             function(r) {
               coda::mcmc(
                 r$Value,
                 start = min(r$Time),
                 end = max(r$Time),
                 thin = unique(diff((
                   r %>% dplyr::arrange(Time)
                 )$Time))
                 )
             })
    )
  }
)

names(RichnessMC) <- unlist(lapply(
  Richness_Binned, function(rb) {
  paste(rb[1, 3:6], collapse = " ")
}))

lapply(RichnessMC, autoburnin=FALSE, coda::gelman.diag)