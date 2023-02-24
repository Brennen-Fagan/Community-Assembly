# Check all functions used, sort them by whether they are internal to this pkg.

library(NCmisc)

# Assumed this dir is root package directory.
thisdir <- '.'
targdirs <- file.path(thisdir, "experiments", c(
  "Figure3-ExampleOutcomes", "Figure4-MetasimulationStudy", "Tests"
))

files <- dir(
  path = targdirs, pattern = "[.]R$",
  full.names = TRUE, include.dirs = TRUE,
  no.. = TRUE, ignore.case = TRUE,
  recursive = TRUE # should be excessive
)

functionsused <- lapply(files, NCmisc::list.functions.in.file)
