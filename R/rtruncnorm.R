# .rtruncnormDirectory <- getSrcDirectory(function(x){x})
rtruncnorm <- function(lo, hi, mean, std) {
  # if(exists("rtnorm")) {
    rtnorm(lo, hi, mean, std)
  # } else {
  #   tryCatch(
  #     {
  #       message("Loading, please wait.")
  #       source(file.path(.rtruncnormDirectory, "run_rtnorm.R"))
  #       rtnorm(lo, hi, mean, std)
  #     },
  #     error = function(e) {e}
  #   )
  # }
}
