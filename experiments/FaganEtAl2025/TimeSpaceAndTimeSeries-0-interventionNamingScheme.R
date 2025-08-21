# Not a finished function!
interventionNamingScheme <- function(aff, ppa, ipt) {
  aDO <- patchAffinityDictionaryOrigin[aff, ]
  ppDO <- poolpatchDictionaryOrigin[ppa, ]

  if (explicit <- grepl(pattern = "rep", aDO$PatchAffinities)) {
    initState <-
      paste0("(", paste( # NOT PRETTY FOR 10, MAY WANT TO JUST REPORT FUNC CALL
        vals <- retrieveFunction(aDO$PatchAffinities)(ppDO$NumberEnvironments),
        collapse = ", "), ")")
  } else {
    initState <-
      paste0(aDO$PatchAffinities, "(", ppDO$NumberEnvironments, ")")
  }

  if(is.na(ipt)) {return(initState)}

  ipDO <- interventionPatchDictionaryOrigin[ipt, ]

  if (ppDO$NumberEnvironments == 1) {
    finState <- paste0(
      "(", retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments), ")"
    )
  } else if (is.na(ipDO$InterventionLocation) ||
             !explicit ||
             !grepl(pattern = "rep", ipDO$PatchAffinities)) {

    finState <- # InterventionPercentage is a bit of a misnomer!
      paste0(ipDO$InterventionPercentage * 100, "%", ipDO$PatchAffinities)

  } else if (ipDO$InterventionLocation == 0) {# Left

    valsnew <- retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments)

    finState <- paste0("(", paste(
      valsnew[1:(ppDO$NumberEnvironments*ipDO$InterventionPercentage)],
      vals[
        (ppDO$NumberEnvironments*ipDO$InterventionPercentage + 1):
          ppDO$NumberEnvironments],
      collapse = ", ", sep = ", "), ")")

  } else if (ipDO$InterventionLocation == 1) {# Right

    valsnew <- retrieveFunction(ipDO$PatchAffinities)(ppDO$NumberEnvironments)

    finState <- paste0("(", paste(
      vals[1:(ppDO$NumberEnvironments*(1 - ipDO$InterventionPercentage))],
      valsnew[
        (ppDO$NumberEnvironments*(1 - ipDO$InterventionPercentage) + 1):
          ppDO$NumberEnvironments],
      collapse = ", ", sep = ", "), ")")
  }

  return(paste0(initState, "->", finState))
}
