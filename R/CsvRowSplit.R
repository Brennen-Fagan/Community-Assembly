CsvRowSplit <- function(csv) {
  return(as.numeric(unlist(strsplit(csv, split = ", "))))
}