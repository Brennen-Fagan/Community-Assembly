abundanceToOccupancy64 <- function(abundance, epsilon) { 
  # Takes abundance <- results[[run]]$Abundance
  # Returns integer64s from bit64 as character string.
  # Warning: slow, and negatives possible because first bit is sign bit.
  # Note: rightmost bit is 1's bit.
  rows <- 1:nrow(abundance)
  col <- ncol(abundance[, -1])
  states <- vector(mode = "character", length = length(rows))
  for(row in rows) {
    states[row] <- 
      paste0(lapply(
        seq(1, col, by = 64), 
        function(i) {
          as.character(bit64::as.integer64.bitstring(
            paste0(as.numeric(
              abundance[row, -1][i:min(i+64-1, col)] > epsilon
            ), collapse = "")
          ))
        }
      ), collapse = " ")
  }
  return(states)
}
