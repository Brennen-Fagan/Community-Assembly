`%>%`<- dplyr::`%>%` # Old piping was used in this project
ConvertPreparedToBeta <- function(
  prepared, columns, method, presenceabsence, postfixes, indicator
) {
  # Lazily didn't implement. Requires editing combine step.
  stopifnot(length(columns) == 3)
  # Lazily, also didn't remove the dependence on Bootstrap and Type.
  # This should probably be a complement.

  # # Prepared looks like, with indicator as the first grouping variable.
  # bootstrapSamplesPairedBeta_source %>% dplyr::filter(
  #   TimeGapNumber == preferredTimeGap
  # ) %>% dplyr::group_by(
  #   DistanceFromCenterExpRev, Bootstrap, Type
  # )
  # # or like
  # bootstrapSamplesPairedBeta_source %>% dplyr::filter(
  #   DistanceFromCenterExpRev == 0
  # ) %>% dplyr::group_by(
  #   TimeGapNumber, Bootstrap, Type
  # )

  # Generate the Beta Diversities Across the Columns
  retval <- lapply(
    columns,
    function(Column, dat, method) {
      temp <- dat %>% dplyr::group_modify(
        .f = SpeciesStringsToBeta,
        SpeciesColumn = Column,
        AbundanceColumn = Column,
        PresenceAbsence = presenceabsence,
        Method = method
      )
      colnames(temp)[colnames(temp) == "Beta"] <- method
      return(temp)
    },
    dat = prepared,
    method = method
  )

  # Rename in preparation for combining.
  retval <- lapply(
    seq_along(retval),
    function(i, dat, postfixes) {
      indices <- which(colnames(dat[[i]]) %in% c(method, "Meaningless"))
      colnames(dat[[i]])[indices] <-
        paste0(colnames(dat[[i]])[indices], postfixes[i])
      dat[[i]]
    },
    dat = retval, postfixes = postfixes
  )

  # Combine
  retval <- dplyr::full_join(
    retval[[1]],
    retval[[2]],
    by = c(dplyr::all_of(indicator), "Bootstrap", "Type", "Patch1", "Patch2")
  ) %>% dplyr::full_join(
    retval[[3]],
    by = c(dplyr::all_of(indicator), "Bootstrap", "Type", "Patch1", "Patch2")
  ) %>% tidyr::pivot_longer(
    cols = -c(dplyr::all_of(indicator), "Bootstrap", "Type", "Patch1", "Patch2"),
    names_to = c("Measure", "Subset"),
    names_sep = ", "
  ) %>% tidyr::pivot_wider(
    names_from = "Measure", values_from = "value"
  )

  return(retval)
}
