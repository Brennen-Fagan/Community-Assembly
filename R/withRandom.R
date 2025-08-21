# We need a function to generate control random values.
# https://stackoverflow.com/a/59875367 Gwang-Jin Kim and
# https://stackoverflow.com/a/14324316 Romain Francois (same question).
withRandom <- function(expr, seed) {
  if (exists(".Random.seed")) {
    oldSeed <- .Random.seed
    on.exit({.Random.seed <<- oldSeed})
  }
  set.seed(seed)
  expr
}
