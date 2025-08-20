load_safe <- function(fname) {
  loaded <- tryCatch({load(fname)},
                     error = function(e) {
                       print(fname)
                       print(e)
                       return(NA)
                     })
  if (all(is.na(loaded))) {
    return(NA)
  } else {
    return(sapply(loaded, get,
                  envir = sys.frame(sys.parent(0)),
                  simplify = FALSE, USE.NAMES = TRUE))
  }
}