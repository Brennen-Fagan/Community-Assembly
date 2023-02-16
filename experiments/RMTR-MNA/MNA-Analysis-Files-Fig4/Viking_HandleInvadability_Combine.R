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

sets <- c(1:3) # Choices: 1:3.
setsLetters <- LETTERS[sets]

# Problems with X11
options(bitmapType = "cairo")

ribbonQuantiles <- c(0.05, 0.95) # For plotting vs time ribbons.
centralEst <- function(x) mean(x, na.rm = TRUE)
# For plotting vs time central point estimate.

if (length(sets) == 0) stop("No sets chosen.")

# Manually verified the following packages are installed.
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)

# Plotting parameters.
plot_width = 8; plot_height = 6; plot_units = "in"; plot_dpi = 320 # ggsave.

# Locations.
directory <- '.'

inputLocation <- c(
  #file.path(directory, "SaveInvadability_2022-05-27"),
  #file.path(directory, "SaveInvadability_2022-06-20")
  file.path(directory, "SaveInvadability_All")
)
if (!dir.exists(inputLocation)) stop(paste(inputLocation, "does not exist."))

# Cluster.
# cargs <- as.numeric(commandArgs(TRUE))
cargs <- 4

setsFiles <- lapply(
  setsLetters,
  function(sL) {
    c(dir(path = inputLocation,
          pattern = paste0(".+Master", sL, ".+[.]RData$"),
          full.names = TRUE),
      dir(path = inputLocation,
          pattern = paste0(".+HiDisp", sL, ".+[.]RData$"),
          full.names = TRUE))
  }
)

source("Viking_HandleDiversity_HelperFunctions.R")

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
    Invadabilities <- foreach::foreach(
      file = iterators::iter(files),
      .packages = c("dplyr", "ggplot2", "Matrix"),
      .inorder = FALSE,
      .export = toExport # Technically different environments, so be explicit.
    ) %dopar%  {
      gc()

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
      } else if(!exists("Invadability")) {
        return(paste(paste(idNums, collapse = " "), "Invadability Not Found"))
      } else {
        Attributes <- extractAttributes(Invadability, idNums)

        # We neglected to include species presence absence with species invad.
        # This is a problem if we are wanting to do regional analysis, since a
        # region is invadable w.r.t. a species if the species is not present
        # anywhere but could be present somewhere.
        # This is okay for local or occupancy analysis though.
        # At a local scale, if a species is reported to be able to invade,
        # then it must not be present, must have positive per capita growth and
        # must survive and grow at least once over the trial event.
        # With only one place to check, that's fine.
        # Occupancy is similar, as it is just the number of patch-species pairs
        # that are unfulfilled but could be fulfilled.
        # It would be better if we had community size though.
        # If we keep it to raw numbers, it will at least share scales with
        # the local and regional richnesses though.

        frame <- dplyr::bind_rows(
          lapply(
            seq_along(Invadability$Invadabilities),
            function(i, d, a) {
              assignAttributes(
                i = i, a = a,
                data.frame(
                  Environment = 1:10,
                  Local = unlist(lapply(1:10, function(j) {
                    sum(d[[i]]$invadability[1, 1:100 + 100*(j - 1)])
                  }))
                ) %>% dplyr::summarise(
                  Occupancy = sum(Local),
                  LocalAvg = mean(Local),
                  LocalMdn = median(Local),
                  .groups = "drop"
                )
                )
            },
            d = Invadability$Invadabilities,
            a = Attributes
          )
        )

        return(list(
          Attr = Attributes,
          Invs = frame
        ))
      }
    }
    parallel::stopCluster(clust)

    ## Combine Results: ########################################################
    Attributes <- dplyr::bind_rows(lapply(
      Invadabilities, function(d) {d$Attr}
    ))
    Invadabilities <- dplyr::bind_rows(lapply(
      Invadabilities, function(d) {d$Invs}
    ))

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
         Invadabilities,
         file = paste0("Set", setNumber, "-AllInvadabilityData.RData"),
         compress = "bzip2"
    )

    return(TRUE)
  }
)

print("Successes")
print(sets, setsLetters, success)
