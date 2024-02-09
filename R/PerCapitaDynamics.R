PerCapitaDynamics_Type1 <- function(
  ReproductionRate, InteractionMatrix, NumEnvironments
) {
  isvalid <- function(func) {
    forms <- names(formals(func))
    if ("..." %in% forms || all(c("t", "y", "parms") %in% forms)) {TRUE}
    else {FALSE}
  }
  
  # Parms contains r = ReproductionRate, a = InteractionMatrix
  if (is.function(ReproductionRate)) {
    stopifnot(isvalid(ReproductionRate))
    
    if (is.function(InteractionMatrix)) {
      stopifnot(isvalid(InteractionMatrix))
      function(t, y, parms = NULL) {
        unlist(lapply(1:NumEnvironments, function(i) {
          ReproductionRate(t = t, y = y, parms = c(parms, Patch = i))
        })) +
          InteractionMatrix(t = t, y = y, parms = parms) %*% y
      }
      
    } else {
      function(t, y, parms = NULL) {
        unlist(lapply(1:NumEnvironments, function(i) {
          ReproductionRate(t = t, y = y, parms = c(parms, Patch = i))
        })) +
          InteractionMatrix %*% y
      }
    }
    
  } else {
    if (is.function(InteractionMatrix)) {
      stopifnot(isvalid(InteractionMatrix))
      function(t, y, parms = NULL) {
        rep(ReproductionRate, NumEnvironments) +
          InteractionMatrix(t = t, y = y, parms = parms) %*% y
      }
      
    } else {
      function(t, y, parms = NULL) {
        rep(ReproductionRate, NumEnvironments) +
          InteractionMatrix %*% y
      }
    }
  }
}