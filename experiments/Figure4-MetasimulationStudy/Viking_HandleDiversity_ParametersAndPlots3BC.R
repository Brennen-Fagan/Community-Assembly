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
sets <- c(1, 3) # Choices: 1:3.
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
time_grouping_size <- 20 # Used for initial filtering down of data, median.
time_averaging_size <- 4 # Used for averaging after grouping has been done.
time_burnin <- 1 # Discard data before this time for summarisation.
time_burnout <- 6.5 # Discard data after this time for summarisation.
time_units <- 1E4 # Inherited from diversity creation. Controls time axis label.
plot_width = 8; plot_height = 6; plot_units = "in"; plot_dpi = 320 # ggsave.

# Locations.
directory <- '.'
#outputLocation <- file.path(directory, paste0("SaveImages_", Sys.Date()))
#if (!dir.exists(outputLocation)) {
#  dir.create(outputLocation, showWarnings = FALSE)
#}

inputLocation <- c(
  # file.path(directory, "SaveDiversityBC_2022-06-17"),
  # file.path(directory, "SaveDiversityBC_2022-06-20"),
  # file.path(directory, "SaveDiversityBC_2022-06-28"),
  # file.path(directory, "SaveDiversityBC_2022-07-15"),
  # file.path(directory, "SaveDiversityBC_2022-07-25"),
  file.path(directory, "SaveDiversityBC_All")
)
if (!all(dir.exists(inputLocation))) stop(paste(inputLocation, "does not exist."))

plottingDataLocation <- file.path(
  directory,
  paste0("SavePlotData_All_", time_grouping_size, "_", time_averaging_size))
if (!dir.exists(plottingDataLocation)) {
  dir.create(plottingDataLocation, showWarnings = FALSE)
}

# Cluster.
cargs <- as.numeric(commandArgs(TRUE))
# cargs <- 14

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

source("Viking_HandleDiversity_HelperFunctionsBC.R")

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
    Diversities <- foreach::foreach(
      file = iterators::iter(rev(files)), 
      .packages = c("dplyr", "ggplot2"),
      .inorder = FALSE,
      .export = toExport # Technically different environments, so be explicit.
    ) %dopar%  {
      gc()

      prospectiveFile <- file.path(
        plottingDataLocation,
        paste0(tools::file_path_sans_ext(basename(file)),
               "_", time_grouping_size, "_", time_averaging_size, ".RData")
      )

      if (file.exists(prospectiveFile)) {
        print(paste(prospectiveFile, "exists."))
        tryCatch({
          load(prospectiveFile)

          if (!exists("DiversityMetrics")) {stop("DiversityMetrics does not exist.")}

          return(DiversityMetrics)
        }, error = function(e) {
          print(paste("Failed to load", prospectiveFile, ",", e))
          return(NULL)
        })
      } else {
        print(paste("Generating:", prospectiveFile))
        DiversityMetrics <- loadDiversity(file)
        save(DiversityMetrics, file = prospectiveFile)
        return(DiversityMetrics)
      }
    }
    parallel::stopCluster(clust)

    ## Combine Results: ########################################################
    Attributes <- dplyr::bind_rows(lapply(
      Diversities, function(d) {d$Attr}
    ))
    DiversitiesAlpha <- dplyr::bind_rows(lapply(
      Diversities, function(d) {d$Alpha}
    ))
    DiversitiesBeta <- dplyr::bind_rows(lapply(
      Diversities, function(d) {d$Beta}
    ))
    DiversitiesGamma <- dplyr::bind_rows(lapply(
      Diversities, function(d) {d$Gamma}
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

    # ## Create Plotting Data Frames: ############################################
    # DiversitiesAlphaGamma <- dplyr::full_join(
    #   DiversitiesAlpha,
    #   DiversitiesGamma,
    #   by = c(
    #     "Time", "Set", "Number", "History",
    #     "Pool", "Noise", "Neutral", "Space"
    #   ),
    #   suffix = c(", Alpha", ", Gamma")
    # )
    # 
    # ## Create Helper Labels: ###################################################
    # levelsPoolNoise <- paste(
    #   DiversitiesGamma$Pool,
    #   DiversitiesGamma$Noise,
    #   sep = ";;"
    # ) %>% unique() %>% strsplit(
    #   split = ";;", fixed = TRUE
    # )
    # 
    # spaceToDispersal <- data.frame(
    #   Space = unique(Attributes$Space)
    # ) %>% dplyr::mutate(
    #   Dispersal = 1 - exp( -9 / as.numeric(Space) ),
    #   Dispersal = paste0(
    #     formatC(Dispersal))
    # ) %>% dplyr::arrange(Space)
    # 
    # neutralToLambdas <- data.frame(
    #   Neutral = unique(Attributes$Neutral)
    # ) %>% tidyr::separate(
    #   Neutral, into = c("Immigration", "Extinction"), remove = FALSE,
    #   sep = "[,][ ]"
    # ) %>% dplyr::mutate(
    #   Immigration = paste0(Immigration, " %*% imm."),
    #   Extinction = paste0(Extinction, " %*% ext.")
    # )
    # 
    # ## Plotting: ###############################################################


    # Perhaps the riskiest bit.
    save(Attributes,
         DiversitiesAlpha,
         DiversitiesBeta,
         DiversitiesGamma,
         file = paste0("Set", setNumber, "-AllLoadedDataBC.RData"),
         compress = "bzip2"
    )

    return(TRUE)
  }
)

print("Successes")
print(sets, setsLetters, success)
