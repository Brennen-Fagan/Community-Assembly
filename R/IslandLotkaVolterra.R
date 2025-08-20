IslandLotkaVolterra <- function(t, y, parms) {
  with(as.list(parms), {
    if (exists("Verbose")) {
      if (Verbose) print(paste(t, "Present:", toString(which(y > 0))))
    }
    retval <- list(as.numeric(y * (r + a %*% y) + d %*% y))
    # as.numeric since the solver doesn't know what Matrix::Matrices are.
    if (exists("Verbose")) {
      if (Verbose) print(paste(t, "Pos-Derivatives:",
                               toString(which(retval[[1]] > 0))))
    }
    return(retval)
  })
}