retrieveFunction <- function(funcstring) {
  funcstring <- strsplit(funcstring, split = "::")
  if (length(funcstring) > 1) {
    stop(paste0("Too many functions provided in string: ", length(funcstring)))
  } else {
    funcstring <- funcstring[[1]]
  }
  if (length(funcstring) > 2) {
    stop(paste0("Too many parts to function provided: ",
                length(funcstring)))
  } else if (length(funcstring) == 2) {
    funcstring <- get(funcstring[2], envir = loadNamespace(funcstring[1]))
  } else if (length(funcstring) == 1) {
    funcstring <- get(funcstring[1])
  } else {
    stop("No parts found for function.")
  }
  return(funcstring)
}
