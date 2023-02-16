# Follow up to MNA-Image-ExampleOutcome-Presence-Fig3.R
# The question is the distribution of inter-arrival times of events,
# with a particular eye for how the various types and scales interact.

# The premise email:
# What I've done (to make sure we're on the same page):
#   For each of the time series in Figure 3, for local, regional, and pool
#       (null/all neutral events) scales
#   I have extracted the appearances and disappearances of species within
#       patches as well as across the system and double-checked them against
#       the list of neutral events to categorize them.
#   I then arranged them by time, and then taken the time differences as the
#       sample(s) of interest.
#   I then fit a few distributions (exponential, gamma, Levy) via MLE
#       (directly via estimators or indirectly via optimisation) before
#       comparing fits. See attached.
# The obvious conclusions are...


library(dplyr)        # Data Manipulation
library(tidyr)        # Data Pivotting
library(ggplot2)      # 2-D Plot
library(gridExtra)    # 2-D Plotting utilities
library(fitdistrplus) # Easy fitting distributions
library(rmutil)       # Levy Distributions

# Problems with X11
options(bitmapType = "cairo")

by_for_thinning <- 1 # 1 row (time step) per ... rows.
divide_time_by <- 1 # time units
burn_in <- 1E4 # time units

# Functions: ###################################################################

load_safe_thin <- function(fname, bythin, divtime, burn) {
  loaded <- tryCatch({load(fname)},
                     error = function(e) {
                       print(fname)
                       print(e)
                       return(NA)
                     })
  if (is.na(loaded)) {
    return(NA)
  } else {
    loaded <- get(loaded)
    loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                             to = nrow(loaded$Abundance),
                                             by = bythin), ]

    toEliminate <-
      loaded$Abundance[, -1] <
      loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
    loaded$Abundance[, -1][toEliminate] <- 0

    loaded$Abundance <- loaded$Abundance[
      loaded$Abundance[, 1] > burn,
    ]

    loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divtime
    return(loaded)
  }
}

load_safe <- function(fname) {
  loaded <- tryCatch({load(fname)},
                     error = function(e) {
                       print(fname)
                       print(e)
                       return(NA)
                     })
  if (all(is.na(loaded))) {
    return(NA)
  } else {
    return(sapply(loaded, get,
                  envir = sys.frame(sys.parent(0)),
                  simplify = FALSE, USE.NAMES = TRUE))
  }
}

FindPresenceChanges <- function(loaded) {
  abund <- loaded$Abundance
  abund[, -1] <- abund[, -1] > loaded$Parameters$EliminationThreshold
  numSpecies <- (ncol(abund) - 1) / loaded$NumEnvironments

  reports <- lapply(1:numSpecies, function(i, a) {
    # Looking for a data.frame(
    #  Time, Species, Environment, Type = c(Immigration, Extirpation),
    #  Regional = c(True, False), Neutral = c(True, False)
    # )

    time <- a[, 1]
    x <- a[, i + numSpecies * (1:loaded$NumEnvironments - 1) + 1]
    xdiff <- apply(x, 2, function(y) y - dplyr::lag(y))

    # First time after a change occurs. Note this double lists.
    changelocs <- apply(xdiff, 2, function(y) list(which(y != 0)))
    if(length(unlist(changelocs)) == 0) {
      return(NULL)
    }

    reportspatch <- lapply(1:loaded$NumEnvironments, function(j, y, locs) {
      if (length(locs[[j]][[1]]) == 0) return(NULL)
      stopifnot(((as.numeric(colnames(y)[j]) - 1) %/% numSpecies) + 1 == j)

      data.frame(
        Time = time[locs[[j]][[1]] - 1],
        Species = ((i - 1) %% numSpecies) + 1,
        Environment = as.character(j),
        Type = dplyr::case_when(
          y[locs[[j]][[1]], j] ==  1 ~ "Immigration",
          y[locs[[j]][[1]], j] == -1 ~ "Extirpation",
          TRUE ~ "OOPS"
        )
      )
    }, y = xdiff, locs = changelocs) %>% dplyr::bind_rows()

    changelocs <- sort(unique(unlist(changelocs)))

    # Regional Report:
    changesregional <- lapply(changelocs, function(loc, y) {
      # For each change loc, look at the timestep before and of.
      target <- y[-1:0 + loc, ]
      # If all 0 before and any 1 of -> "Immigration",
      # If any 1 before and all 0 of -> "Extirpation"
      # Otherwise, discard.
      dplyr::case_when(
        all(target[1, ] == 0) && any(target[2, ] == 1) ~ "Immigration",
        any(target[1, ] == 1) && all(target[2, ] == 0) ~ "Extirpation",
        TRUE ~ "Discard"
      )
    }, y = x)

    # print(changelocs)
    # print(time[changelocs])
    # print(changesregional)

    reportsregion <- data.frame(
      Time = time[changelocs - 1],
      Species = ((i - 1) %% numSpecies) + 1,
      Type = unlist(changesregional),
      Regional = TRUE
    ) %>% dplyr::filter(Type != "Discard")

    reports <- dplyr::full_join(
      reportspatch, reportsregion,
      by = c("Time", "Species", "Type")
    ) %>% dplyr::mutate(
      Regional = ifelse(is.na(Regional), FALSE, Regional)
    )

    return(reports)

  }, a = abund)

  dplyr::bind_rows(reports)
}

addNeutral <- function(loaded, presenceChanges) {
  # Make Consistent
  events <- loaded$Events %>% dplyr::filter(
    Times > burn_in
  ) %>% dplyr::mutate(
    Times = Times / divide_time_by,
    Type = dplyr::case_when(
      Type == "Extinct" ~ "Extirpation",
      Type == "Arrival" ~ "Immigration",
      TRUE ~ "OOPS"
    ),
    Environment = as.character(Environment)
  ) %>% dplyr::rename(
    Time = Times
  )

  events %>% dplyr::full_join(
    presenceChanges, by = c("Time", "Species", "Environment", "Type")
  ) %>% dplyr::mutate(
    Neutral = !is.na(Success),
    Success = ifelse(is.na(Success), TRUE, Success),
    Regional = ifelse(is.na(Regional), FALSE, Regional)
  )
}

CalcInterarrivals <- function(presenceChanges) {
  presenceChanges %>% dplyr::group_by(
    Simulation
  ) %>% dplyr::arrange(
    Time
  ) %>% dplyr::mutate(
    LagTime = dplyr::lag(Time),
    WaitTime = Time - LagTime,
    LagType = dplyr::lag(Type)
  )

}

not <- function(f) function(...) !f(...)

# Implementation: ##############################################################

### Setup: #####################################################################
# All .RData
files_dat <- dir(
  path = "Data_2022-09-16",
  pattern = "MNA[-]ExampleOutcome[-].+[.]RData$",
  full.names = TRUE
)

# Separate out PoolMats
files_dat_PM <- files_dat[
  grepl(x = files_dat,
        pattern = "PoolMats",
        fixed = TRUE)
]
files_dat <- files_dat[
  !grepl(x = files_dat,
         pattern = "PoolMats",
         fixed = TRUE)
]

### Load: ######################################################################
Results <- sapply(
  files_dat,
  load_safe_thin,
  bythin = by_for_thinning,
  divtime = divide_time_by,
  burn = burn_in,
  simplify = FALSE, USE.NAMES = TRUE
)

names(Results) <- basename(names(Results))

PoolsMats <- sapply(
  files_dat_PM,
  load_safe,
  simplify = FALSE, USE.NAMES = TRUE
)

numSpecies <- table(PoolsMats[[1]]$Pool$Type)

presenceChanges <- lapply(Results, function(res) {
  pc <- FindPresenceChanges(res)
  pc <- addNeutral(res, pc)
  return(pc)
})

### Properties: ################################################################
Properties <- strsplit(names(presenceChanges), '-',
                       fixed = TRUE)
# 1st Chunk: Name, Discard
# 2nd Chunk: Iteration + Distance
# 3rd Chunk: Result or Extinction Rate (or Arrival?)
# 4th Chunk: Number of Environments
# 5th Chunk: Space Type + .RData
# Note the mix of Keyword and Location Structure (oops).
# Note also that this strsplit character is a bad decision and should be changed for next time. (D'oh.)
# (E.g. Dates DD-MM-YYYY, Decimals 1.35e-05.)
Properties <- data.frame(
  do.call(rbind, Properties),
  stringsAsFactors = FALSE
)
names(Properties)[1:6] <- c(
  "MNA", "IterANDRate", "Modifier", "EnvNum", "Space", "DistAND.RData"
)

Properties$FullName <- names(presenceChanges)

# Capture the position between the text (first group)
# and the set of numbers (somehow without the +).
# The \\K resets so that we do not capture any text.
patternString <- "((?>[a-zA-Z]+)(?=[0-9eE]))\\K"

# Split strings. Some of the trick will be to introduce
# a character to make the separation around. We use "_".
Properties <- Properties %>% dplyr::mutate(
  IterANDRate = gsub(pattern = patternString,
                     replacement = "_",
                     x = IterANDRate, perl = TRUE),
  Modifier = gsub(pattern = patternString,
                  replacement = "_",
                  x = Modifier, perl = TRUE),
  EnvNum = gsub(pattern = patternString,
                replacement = "_",
                x = EnvNum, perl = TRUE)
) %>% tidyr::separate(
  IterANDRate, into = c("Iter", "Rate"),
  sep = "[_]", fill = "right"
) %>% tidyr::separate(
  Modifier, into = c("Modifier", "ModIntensity"),
  sep = "[_]", fill = "right"
) %>% tidyr::separate(
  EnvNum, into = c("Env", "Environments"),
  sep = "[_]"
) %>% tidyr::separate(
  DistAND.RData, into = c("Distance", ".RData"),
  sep = "[.]"
) %>% dplyr::select(
  -MNA, -.RData, -Env
) %>% dplyr::mutate(
  Distance = dplyr::case_when(
    is.na(Distance) ~ "1e+00",
    TRUE ~ Distance
  )
)

### Add Properties to Loaded: ##################################################
presenceChanges <- lapply(1:length(presenceChanges),
                    function(i, df, nm) {
                      df[[i]] %>% dplyr::mutate(
                        Simulation = nm[i]
                      )
                    },
                    df = presenceChanges,
                    nm = names(presenceChanges))

presenceChanges <- dplyr::left_join(
  dplyr::bind_rows(presenceChanges),
  Properties,
  by = c("Simulation" = "FullName")
) %>% dplyr::mutate(
  # Distance = ifelse(
  #   Space == "Line" | Space == "Ring", 1, Inf
  # ),
  Distance = 10^as.numeric(Distance),
  Space = Distance,
  Dispersal = 1 - exp( -2 / as.numeric(Space) ),
  Dispersal = paste0(formatC(Dispersal)),
  Dispersal = factor(Dispersal, ordered = TRUE, levels = rev(unique(Dispersal)))
)

dispersalLevels <- levels(presenceChanges$Dispersal)

### Subsets for analysis: ######################################################

presenceChanges2 <- presenceChanges %>% dplyr::group_by(
  Time, Species, Type,
  Simulation, Iter, Rate, Modifier, ModIntensity,
  Environments, Space, Distance, Dispersal
) %>% dplyr::summarise(
  Environment = paste(sort(Environment), collapse = ", "),
  Success = any(Success),
  Regional = any(Regional),
  Neutral = any(Neutral)
) %>% dplyr::ungroup(
)

Neutral <- presenceChanges %>% dplyr::filter(
  Neutral == TRUE
) %>% CalcInterarrivals()
Regional <- presenceChanges %>% dplyr::filter(
  Regional == TRUE, Success == TRUE
  ) %>% CalcInterarrivals()
Local <- presenceChanges %>% dplyr::filter(
  Success == TRUE
  ) %>% CalcInterarrivals()

together <- list(Neutral, Regional, Local)
names(together) <- c("Neutral", "Regional", "Local")

minimumWaitTimes <-
  lapply(Results, function(x) min(diff(x$Abundance[, 1])) / 1)

Neutral <- Neutral %>% dplyr::mutate(
  WaitTime = dplyr::case_when(
    WaitTime != 0 ~ WaitTime,
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-0.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-0.RData`,
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-5.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-5.RData`,
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-Inf.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-Inf.RData`,
    TRUE ~ -1
  )
)
Regional <- Regional %>% dplyr::mutate(
  WaitTime = dplyr::case_when(
    WaitTime != 0 ~ WaitTime,
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-0.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-0.RData` * runif(dplyr::n()),
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-5.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-5.RData` * runif(dplyr::n()),
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-Inf.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-Inf.RData`* runif(dplyr::n()),
    TRUE ~ -1
  )
)
Local <- Local %>% dplyr::mutate(
  WaitTime = dplyr::case_when(
    WaitTime != 0 ~ WaitTime,
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-0.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-0.RData` * runif(dplyr::n()),
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-5.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-5.RData` * runif(dplyr::n()),
    Simulation == "MNA-ExampleOutcome-Result-Env10-Ring-Inf.RData" ~
      minimumWaitTimes$`MNA-ExampleOutcome-Result-Env10-Ring-Inf.RData` * runif(dplyr::n()),
    TRUE ~ -1
  )
)

write.csv(Neutral, file = "Neutral.csv")
write.csv(Regional, file = "Regional.csv")
write.csv(Local, file = "Local.csv")

ggplot2::ggplot(Regional, ggplot2::aes(
  x = WaitTime
)) + ggplot2::geom_density(
) + ggplot2::facet_wrap(. ~ Dispersal)

with(Regional, table(Type, LagType, Dispersal))

fitsMME <- lapply(seq_along(together), function(i, presChanges) {
  temp <- lapply(dispersalLevels, function(dLevel, pc, nm) {

    tempexp <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "exp", method = "mme"
    )
    tempgam <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "gamma", method = "mme"
    )

    cdfcomp(list(tempexp, tempgam),
            legendtext = c("exp", "gamma"),
            main = paste(nm, dLevel, "CDFs"))

    return(list("exp" = tempexp, "gamma" = tempgam))

  }, pc = presChanges[[i]], nm = names(presChanges)[i])
  names(temp) <- dispersalLevels
  temp
}, presChanges = together)


fitsMLE <- lapply(seq_along(together), function(i, presChanges) {
  temp <- lapply(dispersalLevels, function(dLevel, pc, nm) {

    tempexp <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "exp", method = "mle"
    )
    tempgam <- tryCatch(pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "gamma", method = "mle",
      start = c(list("shape" = 2), as.list(tempexp$estimate)),
      lower = c(0, 0)
    ), error = function(e) {print(e);return(NULL)})

    retval <- list("exp" = tempexp, "gamma" = tempgam)

    cdfcomp(retval[unlist(lapply(retval, not(is.null)))],
            legendtext =
              c("exp", "gamma")[unlist(lapply(retval, not(is.null)))],
            main = paste(nm, dLevel, "CDFs"))

    return(retval)

  }, pc = presChanges[[i]], nm = names(presChanges)[i])
  names(temp) <- dispersalLevels
  temp
}, presChanges = together)

fitsBoth <- lapply(seq_along(together), function(i, presChanges) {
  temp <- lapply(dispersalLevels, function(dLevel, pc, nm) {

    tempexp <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "exp", method = "mme"
    )
    tempgam <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "gamma", method = "mme"
    )

    tempexp2 <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "exp", method = "mle"
    )
    tempgam2 <- tryCatch(pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "gamma", method = "mle",
      start = as.list(tempgam$estimate),
      lower = c(0, 0)
    ), error = function(e) {print(e);return(NULL)})

    tempexp3 <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "exp", method = "qme", probs = 0.5
    )
    tempgam3 <- tryCatch(pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "gamma", method = "qme", probs = c(0.25, 0.75)
    ), error = function(e) {print(e);return(NULL)})

    tempexp4 <- pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "exp", method = "mge", gof = 'KS'
    )
    tempgam4 <- tryCatch(pc %>% dplyr::filter(
      !is.na(WaitTime),
      Dispersal == dLevel
    ) %>% dplyr::pull(WaitTime) %>% fitdistrplus::fitdist(
      distr = "gamma", method = "mge", gof = 'KS'
    ), error = function(e) {print(e);return(NULL)})

    retval <- list("expMME" = tempexp, "gammaMME" = tempgam,
                   "expMLE" = tempexp2, "gammaMLE" = tempgam2,
                   "expQME" = tempexp3, "gammaQME" = tempgam3,
                   "expMGE" = tempexp4, "gammaMGE" = tempgam4)

    # cdfcomp(retval[unlist(lapply(retval, not(is.null)))],
    #         legendtext =
    #           names(retval)[unlist(lapply(retval, not(is.null)))],
    #         main = paste(nm, dLevel, "CDFs"))

    return(retval)

  }, pc = presChanges[[i]], nm = names(presChanges)[i])
  names(temp) <- dispersalLevels
  temp
}, presChanges = together)



