# Okay, next, one common topic is essentially native vs invasive.
# In truth all species are more or less the same, originating from the true
# regional pool. Our researchers wouldn't know this.
# Instead, they would both use their control to estimate the natives.
computeSpeciesInControl <- function(sampling,
                                    Time = "Time",
                                    IDColumn = "SamplingIDs",
                                    OutPrefix = "Sampling") {
  IDColNum <- which(colnames(sampling) == IDColumn)

  controlSpecies <- sampling %>% dplyr::filter(
    Control == "Control"
  ) %>% dplyr::pull(
    IDColNum
  ) %>% strsplit(
    ", ", fixed = TRUE
  ) %>% unlist(
  ) %>% unique()

  sampling %>% dplyr::group_by(
    dplyr::across(dplyr::all_of(Time)), Patch
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      splitIDs <- strsplit(x[[IDColumn]], ", ", fixed = T)[[1]]
      iDsInControl <- splitIDs %in% controlSpecies

      # Note: group_by (Patch) enforces Alpha / Local scale
      x <- x %>% dplyr::mutate(
        IDsNative = toString(splitIDs[iDsInControl]),
        AbundanceNative = sum(iDsInControl),
        AlphaNative = length(unique(splitIDs[iDsInControl])),
        IDsInvasive = toString(splitIDs[!iDsInControl]),
        AbundanceInvasive = sum(!iDsInControl),
        AlphaInvasive = length(unique(splitIDs[!iDsInControl]))
      )

      colnames(x)[(ncol(x)-5):ncol(x)] <-
        paste0(OutPrefix, colnames(x)[(ncol(x)-5):ncol(x)])

      return(x)
    }
  )
}
