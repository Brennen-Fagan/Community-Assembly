
#datfolders <- dir(pattern = "TSTS_Simulations_.+2025-07-30")
#datfolders <- dir(pattern = "TSTS_Simulations_.+2025-09-04")
datfolders <- dir(pattern = "TSTS_Simulations_.+2026-03-11")

overwrite <- TRUE # FALSE # TRUE

prefix <- "diversitiesFlattened10"
#suffix <- "2025-07-30"
#suffix <- "2025-09-04"
suffix <- "2026-03-11"

# Libraries: ##################################################################
directory <- '.'
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

#update.packages(ask = FALSE, repos = c(CRAN = "https://cloud.r-project.org"), 
#                instlib = librarypath, oldPkgs = "dplyr")

source("TimeSpaceAndTimeSeries-10-Dictionaries.R")
#source('TimeSpaceAndTimeSeries-0-Functions.R')
library(tidytable)
library(RMTRCode2)

# functions: ##################################################################

source(file.path("R", "flattenDiversity.R"))
source(file.path("R", "interventionNamingScheme.R"))

# Load Data: ##################################################################

for (datfolder in datfolders) {
  datfolderID <- paste0(
    strsplit(datfolder, split = "_")[[1]][-c(1:2)],
    collapse = "_")

  filestring <- paste0(prefix, "_", datfolderID,".RData")

  if (file.exists(filestring) && !overwrite) {next()}

  diversities <- lapply(
    dir(datfolder, full.names = TRUE, pattern = "Diversity"), function(x) {
      names <- load(x)
      stopifnot(length(names) == 1)
      obj <- get(names)
      return(c(obj, "Dir" = dirname(x), "File" = basename(x)))
    })

  diversitiesFlattened <- vector(mode = "list", length = length(diversities))
  for(i in 1:length(diversitiesFlattened)) {
    # Pop front of diversities, process, put in flattened, remove.
    # Hence, len(diversities) changes, but pre-alloc => len(flattened) is not.
    diversitiesFlattened[[i]] <- flattenDiversity(diversities[[1]])
    diversities[[1]] <- NULL
    gc()
  }

  rm(diversities)

  diversitiesFlattened <- do.call(rbind, diversitiesFlattened)

  # Human readable patch affinities
  diversitiesInterventionStrings <- diversitiesFlattened %>% dplyr::select(
    PatchAffinity, PoolPatch, InterventionPatchType
  ) %>% dplyr::distinct(
  ) %>% dplyr::mutate(
    Intervention = unlist(mapply(
      FUN = interventionNamingScheme,
      PatchAffinity, PoolPatch, InterventionPatchType
    ))
  )

  # Col-wise append
  diversitiesFlattened <- diversitiesFlattened %>% dplyr::left_join(
    diversitiesInterventionStrings,
    by = c("PatchAffinity", "PoolPatch", "InterventionPatchType"),
    multiple = "all"
  )

  # Human readable species affinities
  diversitiesFlattened <- diversitiesFlattened %>% dplyr::mutate(
    SpeciesPreferences =
      speciesAffinityDictionaryOrigin$SpeciesAffinities[as.numeric(SpeciesAffinity)]
  )

  # Correct the NA for richness values
  diversitiesFlattened <- diversitiesFlattened %>% dplyr::mutate(
    Value = dplyr::case_when(
      Metric == "Alpha Hill:0" & is.na(Value) ~ 0,
      TRUE ~ Value
    )
  )

  save(diversitiesFlattened,
       file = filestring)

  rm(diversitiesFlattened)
  gc()
}

print("files checked")

diversitiesAll <- NULL
for (datfolder in datfolders) {
  datfolderID <- paste0(
    strsplit(datfolder, split = "_")[[1]][-c(1:2)],
    collapse = "_")

  filestring <- paste0(prefix, "_", datfolderID, ".RData")

  stopifnot(file.exists(filestring))

  print(filestring)

  objnames <- load(filestring)
  obj <- get(objnames)

  # We just can't get past the memory barrier, so we'll need to reduce the
  # amount of data we are looking at.
  obj <- obj %>% tidytable::filter(
    # Multiples of 100 for long term behaviour, everything for short term.
    (round(Time, digits = -1) %% 100) == 0 | (Time > 15000 & Time < 18000)
  )

  if(!is.null(diversitiesAll)) {
    diversitiesAll <- tidytable::bind_rows(diversitiesAll, obj)
  } else {
    diversitiesAll <- obj
  }

  rm(obj)
  rm(list = objnames)
  gc()
}

print("divAll computed")

# The all-in-one file.
save(diversitiesAll, file = paste0(prefix, "a1_", suffix, ".RData"))

print("divAll saved")

diversitiesRichness <- diversitiesAll %>% tidytable::filter(
  Metric == "Alpha Hill:0"
)

print("divR computed")

save(diversitiesRichness, file = paste0(prefix, "a1_", suffix, "_Richness.RData"))
rm(diversitiesRichness); gc()

print("divR saved")

diversitiesTimeBC <- diversitiesAll %>% tidytable::filter(
  grepl(x = Metric, pattern = "TimeBrayCurtis")
)

print("divBC computed")

save(diversitiesTimeBC, file = paste0(prefix, "a1_", suffix, "_TimeBC.RData"))

print("divBC saved")
