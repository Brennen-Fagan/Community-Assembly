# Function that automatically increments after use. Note that we need to store
# the result if we want to use the same index twice.
indexFactory <- function() {
  index <- 1
  function() {
    on.exit(index <<- index + 1)
    return(index)
  }
}
