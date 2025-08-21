calculateColExtMetrics <- function(sim) {
  #TODO: I forgot the case of Arrival but Already Present
  #TODO: I seem to have some legitimate approximately triple events,
  #      where the species is present, is about to drop out, and drops out
  #      when an arrival event is supposed to happen.

  # Might be easier to take the binary matrix and to say
  # 0->1 and in the successful events ~ colonization
  # 0->1 and not in the successful events ~ dispersal
  # 1->0 and in the successful events ~ neutral extirpation
  # 1->0 and not in the successful events ~ dynamic extirpation
  # But note that more than one of these can co-occur.
  binaryMatrix <-
    sim$Abundance[, -1] > sim$Parameters$EliminationThreshold
  ColExtMatrix <- # Next - Current so timing of event lines up with the change
    binaryMatrix[2:nrow(sim$Abundance), ] -
    binaryMatrix[1:(nrow(sim$Abundance) - 1), ]

  nspecies <- ncol(ColExtMatrix)

  # Reattach the times
  ColExtMatrix <- cbind(sim$Abundance[1:(nrow(sim$Abundance) - 1), 1],
                        ColExtMatrix)

  allEvents <- lapply(1:nrow(ColExtMatrix), function(i, mat, binmat, events) {
    # print(i)
    # if (i == 160) {
    #   print("stop here")
    # }

    time <- mat[i, 1]
    changes <- which(mat[i, -1] != 0)
    event <- events %>% dplyr::filter(Times == time, Success) # N.B. nrow can >1

    if (nrow(event) == 0 && length(changes) == 0) {
      # Done, nothing to report.
      return(NULL)
    }

    environments <- ((changes - 1) %/% nspecies) + 1 # 1:20 w/10 -> 10 1s, 10 2s
    species <- ((changes - 1) %% nspecies) + 1 # 1:20 w/10 -> c(1:10, 1:10)
    stopifnot(environments <= sim$NumEnvironments)

    if (nrow(event) == 1 && length(changes) == 1 &&
        event$Species == species && event$Environment == environments &&
        ((mat[i, changes + 1] == 1 && event$Type == "Arrival") ||
         (mat[i, changes + 1] == -1 && event$Type == "Extinct")
        )) {
      # Trivial, everything is already recorded
      return(event)
    }

    # Bundle information together for easier processing.
    changesdf <- data.frame(
      change = changes,
      type = mat[i, changes+1], # +1 column for time
      species = species,
      environments = environments
    ) %>% dplyr::mutate(index = dplyr::row_number())

    # Else: complicated, some combination of events or non-neutral events.
    # Proceed along events, removing any changes that are already accounted for.
    # (Essentially, we're doing a complicated implementation of a filtering join
    # where we're trying to filter out things that have already been addressed.)
    # Simultaneously, if an event occurs and the result is not registered as a
    # change, then something must have countered it. If neutral events cannot
    # occur simultaneously, then it must be a non-neutral event.
    if (nrow(event) > 0) {
      for (j in 1:nrow(event)) {
        thisEvent <- event[j, ]
        thisChange <- changesdf %>% dplyr::filter(
          species == thisEvent$Species,
          environments == thisEvent$Environment
        )

        if (nrow(thisChange) == 1 &&
            ((thisChange$type == 1 && thisEvent$Type == "Arrival") ||
             (thisChange$type == -1 && thisEvent$Type == "Extinct")
            )) {
          # This change is already in event(s) as this event.
          # Remove it from the data.frame of event(s) to add.
          changesdf <- changesdf %>% dplyr::filter(index != thisChange$index)
        } else if (
          nrow(thisChange) == 0 # then an event happened, but no change recorded.
          # result depends on whether a population is there or not.
        ) {
          speciesPresent <- with(
            thisEvent,
            binmat[i, Species + (Environment - 1) * nspecies] # No Time Col.
          )
          newEvent <- data.frame(
            Times = time,
            Species = thisEvent$Species,
            Environment = thisEvent$Environment,
            Type = if(thisEvent$Type == "Arrival" && speciesPresent) {
              "Present"
            } else if (thisEvent$Type == "Arrival" && !speciesPresent) {
              "Dynamic Loss"
            } else if (thisEvent$Type == "Extinct" && speciesPresent) {
              "Dispersal"
            } else if (thisEvent$Type == "Extinct" && !speciesPresent) {
              "Dispersal2" # This shouldn't occur, Success should be FALSE.
            } else {
              "Oops"
            },
            Success = TRUE
          )
          if (newEvent$Type == "Present") {
            # Instead, just replace the current event to present.
            # (Excess work, but clearer code!)
            event[j, ]$Type <- "Present"
          } else {
            event <- rbind(event, newEvent)
          }
        }
      }
    }

    # Proceed along changes, adding them to events.
    if (nrow(changesdf) > 0) {
      for (j in 1:nrow(changesdf)) {
        thisChange <- changesdf[j, ]
        thisEvent <- event %>% dplyr::filter(
          Species == thisChange$species,
          Environment == thisChange$environments
        )
        if (nrow(thisEvent) == 0) {
          # Event definitely not recorded, create a new one and add it.
          event <- rbind(event, data.frame(
            Times = time,
            Species = thisChange$species,
            Environment = thisChange$environments,
            Type = switch(thisChange$type + 2, # -1:1 + 2 -> 1:3
                          "Dynamic Loss",
                          "Oops",
                          "Dispersal"
            ),
            Success = TRUE
          ))
        } else if (nrow(thisEvent) == 1 &&
                   thisEvent$Species == thisChange$species &&
                   thisEvent$Environment == thisChange$environments &&
                   ((thisChange$type == 1 && thisEvent$Type == "Arrival") ||
                    (thisChange$type == -1 && thisEvent$Type == "Extinct")
                   )) {
          # Trivial, everything is already recorded.
          # Note this also should already have been caught!
          next()
        } else {
          # Event recorded, but opposite of what is observed -- undone twice.
          # Assuming that events take place at distinct times, which should be
          # guaranteed from the code we are using, this shouldn't happen.
          # Discretisation means that it can happen though, and we can have
          # something considered a false positive: arrival events are tracked as
          # positive (because the species is added to the local population) even
          # if the population is declining and about to be removed.
          if (nrow(thisEvent) == 1 &&
              thisEvent$Species == thisChange$species &&
              thisEvent$Environment == thisChange$environments &&
              (thisChange$type == -1 && thisEvent$Type == "Arrival")
          ) {
            # "Poorly Timed False Positive Arrival"
            event <- rbind(event, data.frame(
              Times = time,
              Species = thisChange$species,
              Environment = thisChange$environments,
              Type = "Dynamic Loss",
              Success = TRUE
            ))
            # We also should account for that the Arrival is actually a Present.
            event[event$Times == thisEvent$Times &
                    event$Species == thisEvent$Species &
                    event$Environment == thisEvent$Environment &
                    event$Type == thisEvent$Type ,]$Type <- "Present"
          } else {
            stop("Double-check assumptions: a 'triple event' occurred.")
          }
        }
      }
    }

    return(event)
  },
  mat = ColExtMatrix,
  binmat = binaryMatrix,
  events = sim$Events)

  allEvents <- dplyr::bind_rows(allEvents)

  endState <- which(sim$Abundance[nrow(sim$Abundance), -1] >
                      sim$Parameters$EliminationThreshold)

  allEvents <- dplyr::bind_rows(
    allEvents,
    # Add in failed colonization events, as requested by CDT.
    sim$Events %>% dplyr::filter(!Success, Type == "Arrival"),
    # Add in faux out events corresponding to remaining in the simulation.
    data.frame(
      Times = sim$Abundance[nrow(sim$Abundance), 1],
      Species = ((endState - 1) %% nspecies) + 1,
      Environment = ((endState - 1) %/% nspecies) + 1,
      Type = "EndOfSimulation",
      Success = TRUE
    )
  ) %>% dplyr::arrange(Times)

  # Adding in typings is then dependent on how we are breaking up affinity
  # and what information from the pool is desired. I think it might make
  # more sense to do that outside of this function.
  return(allEvents)
}
