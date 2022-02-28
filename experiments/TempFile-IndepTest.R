# Independence testing examples
# Without Binning
Diversity %>% dplyr::filter(
  Measurement == "Richness"
) %>% tidyr::pivot_wider(
  id_cols = c("Time", "Modifier", "ModIntensity", "Space", "Distance"), 
  names_from = "Environment",
  values_from = "Value"
) %>% dplyr::group_by(
  Modifier, ModIntensity, Space, Distance
) %>% dplyr::group_split(
)  -> temp

testcorr::rcorr.test(temp[[25]] %>% dplyr::select(`1`:`10`))

# With Binning
Diversity %>% dplyr::filter(
  Measurement == "Richness",
  Time > 3
) %>% dplyr::mutate(
  TimeFloor = floor(Time * 10) / 10
) %>% dplyr::group_by(
  TimeFloor, Environment, Modifier, ModIntensity, Space, Distance
) %>% dplyr::summarise(
  Value = median(Value)
) %>% dplyr::ungroup(
) %>% tidyr::pivot_wider(
  id_cols = c("TimeFloor", "Modifier", "ModIntensity", "Space", "Distance"), 
  names_from = "Environment",
  values_from = "Value"
) %>% dplyr::group_by(
  Modifier, ModIntensity, Space, Distance
) %>% dplyr::group_split(
)  -> temp
testcorr::rcorr.test((temp %>% lapply(FUN = function(x) {x %>% dplyr::select(`1`, `2`:`9`, `10`)}))[[25]])