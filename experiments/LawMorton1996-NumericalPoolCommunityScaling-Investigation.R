dir <- getSrcDirectory(function(x) {})

# Call the relevant settings file.
versionNum <- 3
start <- i <- 7
seedNum <- 13
# 3 -  7 - 13 was "successful", but it looks more like it was buggy(!).
# 3 - 25 - 41 was successful

# load packages

packages <- c(
  "foreach",
  "parallel",
  "doParallel",
  "dplyr",
  "tidyr",
  "Rcpp"
)
for (package in packages) {
  library(package, character.only = TRUE)
}
library(RMTRCode2)

source(file.path(dir,
  paste0("LCAB_LawMorton1996-NumericalPoolCommunityScaling", versionNum),
  paste0(
  "LawMorton1996-NumericalPoolCommunityScaling-Settings", versionNum, ".R"
)))

load(
  file = file.path(dir,
    paste0("LCAB_LawMorton1996-NumericalPoolCommunityScaling", versionNum),
    paste0(
    "LawMorton1996-NumericalPoolCommunityScaling-PoolMats", versionNum, ".RDS"
  ))
)

ourPool <- pools[[i]]
ourMat <- communityMats[[i]]
seed <- seedsRun[((i - 1) * runs + 1) : (i * runs)][seedNum]

aRun <- LawMorton1996_NumericalAssembly(
  Pool = ourPool,
  CommunityMat = ourMat,
  IntegratorTimeStep = 10000,
  InnerTimeStepSize = 5000,
  ArrivalEvents = events,
  ReturnValues = c("Abundance", "Sequence"),
  seed = seed
)
