library(dplyr)

# Directory Functions and Objects: ############################################
directory <- "." # Should be "VariantExperiments"
source(file.path(directory, "TimeSpaceAndTimeSeries-0-Functions.R"))
source(file.path(directory, "TimeSpaceAndTimeSeries-9-Dictionaries.R"))

# Initial results from the 9a round suggest increasing the consumer body range
# upper bound and decrease the consumer consumption range k3.
modifiedcase <- basecase
# modifiedcase$PoolConsumerLogBodySize <- "c(-1, 1)"
modifiedcase$PoolK3ConsumerPredationRange <- 10^-1
modifiedcase$InteractionK3ConsumerPredationRange <- 10^-1
modifiedcase$InteractionK5BasalBiomass <- 30
# modifiedcase$InteractionEliminationThreshold <- 10^-5
# modifiedcase$ColonizationPropaguleSize <-
#   modifiedcase$InteractionEliminationThreshold * 4 * 10^3

# Latin Hypercube Sampling Experiments over the parameters
# There are some queries about some of the parameters I've usually taken for
# granted at this stage, such as the number of species and ratio of basals to
# consumers. There's also enduring questions about if we can play with the
# elminiation threshold.
# > runif(1)*1e8
# [1] 5001003

# Columns of interest, across three dictionaries:
paramsOfInterest <- c(
  # poolpatchDictionaryOrigin
  "BasalConsumerRatio",
  "NSpecies",
  "PoolConsumerLogBodySize",
  # dynamicsDictionaryOrigin
  "InteractionEliminationThreshold",
  # eventsDictionaryOrigin
  "ColonizationPropaguleSize"
)

modifiedcase <- modifiedcase[!names(modifiedcase) %in% paramsOfInterest]

LatinHypercube <- withRandom({
  lhs::optimumLHS(n = 50, # Number of experiments
                  k = length(paramsOfInterest),  # Number of Parameters
                  maxSweeps = 3, eps = 0.1, verbose = FALSE)
}, seed = 5001003)
LatinHypercube <- as.data.frame(LatinHypercube)
colnames(LatinHypercube) <- paramsOfInterest

LatinHypercube <- LatinHypercube %>% dplyr::mutate(
  dplyr::across(.cols = c(dplyr::starts_with("Pool"),
                          "NSpecies", "BasalConsumerRatio"),
                .fns = function(colval) {
                  colname <- dplyr::cur_column()
                  vals <- sort(unique(
                    (poolpatchDictionaryOrigin %>%
                       dplyr::filter(PoolK3ConsumerPredationRange != 10^0)
                    )[[colname]]))
                  colval <- ceiling(colval*length(vals)) # unif -> integer.
                  colval <- vals[colval]
                  return(colval)
                })
)
LatinHypercube <- LatinHypercube %>% dplyr::mutate(
  dplyr::across(.cols = dplyr::starts_with("Interaction"),
                .fns = function(colval) {
                  colname <- dplyr::cur_column()
                  vals <- sort(unique(
                    (dynamicsDictionaryOrigin %>%
                       dplyr::filter(InteractionK3ConsumerPredationRange != 10^0)
                    )[[colname]]))
                  colval <- ceiling(colval*length(vals)) # unif -> integer.
                  colval <- vals[colval]
                  return(colval)
                })
)
LatinHypercube <- LatinHypercube %>% dplyr::mutate(
  dplyr::across(.cols = dplyr::starts_with("Colonization"),
                .fns = function(colval) {
                  colname <- dplyr::cur_column()
                  vals <- sort(unique(eventsDictionaryOrigin[[colname]]))
                  colval <- ceiling(colval*length(vals)) # unif -> integer.
                  colval <- vals[colval]
                  return(colval)
                })
)

experiments <- lapply(1:nrow(LatinHypercube), function(i) {
  params <- LatinHypercube[i, ]

  list(
    ppDO =
      poolpatchDictionaryOrigin %>% dplyr::filter(
        # Standard, and only implemented, LM1996 Function and Parameters,
        PoolDispersalSpeed == 1,
        NumberEnvironments == 1,

        apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
          #             (despite documentation saying otherwise)
          .cols = dplyr::any_of(names(params)),
          .fns = function(colval) {
            colname <- dplyr::cur_column()
            params[[colname]] == colval
          }
        ), 1, all),

        apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
          #             (despite documentation saying otherwise)
          .cols = dplyr::any_of(names(modifiedcase)),
          .fns = function(colval) {
            colname <- dplyr::cur_column()
            modifiedcase[[colname]] == colval
          }
        ), 1, all)
      ),

    dynDO =
      dynamicsDictionaryOrigin %>% dplyr::filter(
        apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
          #             (despite documentation saying otherwise)
          .cols = dplyr::any_of(names(params)),
          .fns = function(colval) {
            colname <- dplyr::cur_column()
            params[[colname]] == colval
          }
        ), 1, all),

        apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
          #             (despite documentation saying otherwise)
          .cols = dplyr::any_of(names(modifiedcase)),
          .fns = function(colval) {
            colname <- dplyr::cur_column()
            modifiedcase[[colname]] == colval
          }
        ), 1, all)
      ),

    eDO =
      eventsDictionaryOrigin %>% dplyr::filter(
        EventsNumberMultiplier == 20, # Longer simulation. Keep number same while
        ImmigrationMultiplier == 0.1, # Decreasing rate of occurrence.
        ExtirpationMultiplier == 0.1,
        ExtirpationProportion == 1, # Extirpation == Extinction in a 1 patch system.
        apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
          #             (despite documentation saying otherwise)
          .cols = dplyr::any_of(names(params)),
          .fns = function(colval) {
            colname <- dplyr::cur_column()
            params[[colname]] == colval
          }
        ), 1, all),

        apply(dplyr::across(# deprecated, but if_all doesn't permit cur_column
          #             (despite documentation saying otherwise)
          .cols = dplyr::any_of(names(modifiedcase)),
          .fns = function(colval) {
            colname <- dplyr::cur_column()
            modifiedcase[[colname]] == colval
          }
        ), 1, all)
      ),

    icDO = initialConditionsDictionaryOrigin %>% dplyr::filter(
      Species == "None"
    ),

    dispDO =
      dispersalDictionaryOrigin %>% dplyr::filter(
        Configuration == "None" # Doesn't make sense for 1 patch systems.
      ),

    distDO =
      distanceDictionaryOrigin %>% dplyr::filter(
        rhofunction %in% c("rho.noop")
      ),

    aDO =
      affinityDictionaryOrigin %>% dplyr::filter(
        SpeciesAffinities %in% c("evensplit_01"),
        PatchAffinities %in% c("rep_0")
      )#,

    # iPDO =
    #   # Note: Every combination of iPDO with ppDO.
    #   interventionPatchDictionaryOrigin %>% dplyr::filter(
    #     PatchAffinities %in% c("rep_0", "rep_0.25", "rep_0.5", "rep_0.75", "rep_1"),
    #     InterventionLocation == 1,
    #     InterventionPercentage == 1
    #   ),
    #
    # iTDO =
    #   interventionTimeDictionaryOrigin %>% dplyr::filter(
    #     Method == "mean"
    #   ),
    #
    # # interventionDispersalDictionaryChoice
    # iDispChoice = "p",
    # # interventionDistanceDictionaryChoice
    # iDistChoice = "p"
  )
})
