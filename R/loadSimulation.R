loadSimulation <- function(filepath1, filepath2 = NULL) {
  file1 <- load(filepath1)
  if (length(file1) > 1) {
    warning("filepath1 has more than one object. Defaulting to first.")
  }
  file1 <- get(file1[1])

  EntriesRequired <- c("Events", "Abundance")
  if (any(! EntriesRequired %in% names(file1))) {
    error("Events or Abundance not found in filepath1 object 1.")
  }

  if (!is.null(filepath2)) {
    file2 <- load(filepath2)
    if (length(file2) > 1) {
      warning("filepath2 has more than one object. Defaulting to first.")
    }
    file2 <- get(file2[1])
    if (any(! EntriesRequired %in% names(file2))) {
      error("Events or Abundance not found in filepath2 object 1.")
    }

    EntriesToCheck <- !names(file2) %in% c("Events", "Abundance")

    stopifnot(isTRUE(all.equal(file1[EntriesToCheck],
                               file2[EntriesToCheck])))

    file1$Abundance <- file1$Abundance[
      file1$Abundance[, 1] < min(file2$Abundance[, 1]),
      ]

    file1$Events <- rbind(file1$Events, file2$Events)
    file1$Events <- file1$Events %>% dplyr::distinct()
    file1$Abundance <- rbind(file1$Abundance, file2$Abundance)
  }

  return(file1)
}
