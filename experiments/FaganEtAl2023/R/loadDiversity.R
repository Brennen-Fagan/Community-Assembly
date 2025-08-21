loadDiversity <- function(file, verbose = TRUE) {
  # Retrieve the trailing id numbers from before the file extension.
  # 1st: Set, 2nd: CaseNumber, 3rd: History, 4th: Part
  # Note, if last two are not present, all histories are bundled together.
  # If last two are present, each file has a single part of a single history.
  idNums <- suppressWarnings(na.omit(
    as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                             split = "-", fixed = TRUE)[[1]],
                    n = 4))),
    classes = "simpleWarning")
  
  if (verbose) {
    print(file)
    print(Sys.time())
  }
  
  tryCatch({
    load(file)
    
    if (!exists("Diversity")) {stop("Diversity does not exist.")}
    
    Attributes <- extractAttributes(Diversity, idNums)
    DiversitiesAlpha <- extractAlphas(Diversity, Attributes)
    DiversitiesBeta <- extractBetas(Diversity, Attributes)
    DiversitiesGamma <- extractGammas(Diversity, Attributes)
    
    return(list(
      Attr = Attributes,
      Alpha = DiversitiesAlpha,
      Beta = DiversitiesBeta,
      Gamma = DiversitiesGamma
    ))
  }, error = function(e) {
    print(paste(paste(idNums, collapse = " "), e, sep = ": "))
    return(NULL)
  })
}