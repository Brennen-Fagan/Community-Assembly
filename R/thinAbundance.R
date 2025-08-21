thinAbundance <- function(abundance, events, threshold,
                          preferred_rows_per_event) {
  bythin <- floor((nrow(abundance)
                   / nrow(events))
                  / preferred_rows_per_event)

  abundance <- abundance[seq(from = 1,
                             to = nrow(abundance),
                             by = bythin), ]

  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- abundance[, -1] < threshold # & abundance[, -1] > 0
  abundance[, -1][toEliminate] <- 0

  return(abundance)
}
