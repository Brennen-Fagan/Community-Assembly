# Viking Adaptation of Server_HandleDiversity_ParametersAndPlots3.R
# Note that we are assuming (fingers crossed!) that all of the parameters were
# collected correctly this time.
# Due to memory constraints, we'll adopt a by-Set approach.
# For a given set, in parallel
#   load the data,
#   process the data,
#   and combine the data,
# and then serially plot the data and save the plots
# before proceeding to the next set.

# Note: when we say case number, we mean parameter id number.

# Setup: #######################################################################
# Control parameters. Neglecting a set can be used to skip that calculation.
sets <- c(1) # Choices: 1:3.
setsLetters <- LETTERS[sets]

# Problems with X11
options(bitmapType = "cairo")

ribbonQuantiles <- c(0.05, 0.95) # For plotting vs time ribbons.
centralEst <- function(x) mean(x, na.rm = TRUE)
# For plotting vs time central point estimate.

if (length(sets) == 0) stop("No sets chosen.")

librarypath <- file.path(".", "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

allLibraryPaths <- .libPaths()

# Manually verified the following packages are installed.
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)
library(RMTRCode2)

source("HelperFunctions.R") # local override with the patch so don't rebuild.

# Plotting parameters.
plot_width = 8; plot_height = 6; plot_units = "in"; plot_dpi = 320 # ggsave.

# Locations.
directory <- '.'

inputLocation <- c(
  #file.path(directory, "SaveInvadability_2022-05-27"),
  file.path(directory, "SaveTimeJaccard_2023-03-06")
)
if (!dir.exists(inputLocation)) stop(paste(inputLocation, "does not exist."))

# Cluster.
cargs <- as.numeric(commandArgs(TRUE))
# cargs <- 4

setsFiles <- lapply(
  setsLetters,
  function(sL) {
    c(#dir(path = inputLocation,
      #    pattern = paste0(".+Master", sL, ".+[.]RData$"),
      #    full.names = TRUE),
      #dir(path = inputLocation,
      #    pattern = paste0(".+HiDisp", sL, ".+[.]RData$"),
      #    full.names = TRUE),
      dir(path = inputLocation,
          pattern = paste0(".+Hik6", sL, ".+[.]RData$"),
          full.names = TRUE))
  }
)

# Main: ########################################################################
toExport <- ls()

success <- lapply(
  setsFiles,
  # List whose entries correspond to sets and contain associated files.
  # Each File's name contains its set, case number, history, and part.
  function(files) {
    ## Parallel Load: ##########################################################

    clust <- parallel::makeCluster(cargs[1], outfile = "")
    doParallel::registerDoParallel(clust)
    TimeJaccards <- foreach::foreach(
      file = iterators::iter(files),
      .packages = c("dplyr", "Matrix"), #, "RMTRCode2"),
      .inorder = FALSE,
      .export = toExport # Technically different environments, so be explicit.
    ) %do% {#par%  {
      gc()
      
      .libPaths(allLibraryPaths)
      library(RMTRCode2)


      idNums <- suppressWarnings(na.omit(
        as.numeric(tail(strsplit(tools::file_path_sans_ext(basename(file)),
                                 split = "-", fixed = TRUE)[[1]],
                        n = 4))),
        classes = "simpleWarning")

      temp <- tryCatch(load(file),
                       error = function(e) {
                         paste(file, e)
                       })

      if ("error" %in% class(temp)) {
        return(temp)
      } else if(!exists("TimeJaccard")) {
        return(paste(paste(idNums, collapse = " "), "TimeJaccard Not Found"))
      } else {
        Attributes <- extractAttributes(TimeJaccard, idNums)

        frame <- dplyr::bind_rows(
          lapply(
            seq_along(TimeJaccard$Betas),
            function(i, d, a) {
              assignAttributes(
                i = i, a = a,
                d[[i]] %>% dplyr::filter("Environment" != "Mean")
              )
            },
            d = TimeJaccard$Betas,
            a = Attributes
          )
        # ) %>% dplyr::group_by( #Cannot do here, some are split in two files.
        #   # Basal, Consumer, or All, and Attributes
        #   Measurement, Set, Number, History, Pool, Noise, Neutral, Space
        #   ) %>% dplyr::summarise(
        #   # Summarising over Environment and Time.
        #   Avg = mean(Value, na.rm = TRUE),
        #   Mdn = median(Value, na.rm = TRUE),
        #   .groups = "drop"
        )


        return(list(
          Attr = Attributes,
          Jacs = frame
        ))
      }
    }
    parallel::stopCluster(clust)

    ## Combine Results: ########################################################
    Attributes <- dplyr::bind_rows(lapply(
      TimeJaccards, function(d) {d$Attr}
    ))
    TimeJaccards <- dplyr::bind_rows(lapply(
      TimeJaccards, function(d) {d$Jacs}
    )) %>% dplyr::group_by(
      # Do it here because files spread out, not a problem for invadability.
      # Basal, Consumer, or All, and Attributes
      Measurement, Set, Number, History, Pool, Noise, Neutral, Space
    ) %>% dplyr::summarise(
      # Summarising over Environment and Time.
      Avg = mean(Value, na.rm = TRUE),
      Mdn = median(Value, na.rm = TRUE),
      .groups = "drop"
    )

    ## Report Results: #########################################################
    print("Pool")
    print(table(Attributes$Pool))
    print("Noise")
    print(table(Attributes$Noise))
    print("Neutral")
    print(table(Attributes$Neutral))
    print("Space")
    print(table(Attributes$Space))

    # Test to make sure nothing has gone wrong.
    setNumber <- unique(Attributes$Set)
    if (length(setNumber) != 1) warning(
      paste(
        "setNumber contains", length(setNumber),
        if (length(setNumber > 1)) {
          paste(" entries:", paste0(setNumber, collapse = ", "))
        } else {"entries."}
      )
    )

    # Perhaps the riskiest bit.
    save(Attributes,
         TimeJaccards,
         file = paste0("Set", setNumber+3, "-AlTimeJacData.RData"),
         compress = "bzip2"
    )

    return(TRUE)
  }
)

# print("Successes")
# print(sets, setsLetters, success)
