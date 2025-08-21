# NOTE: Breaks idiom; should be renamed to be a verb.
SpeciesStringsToBeta <- function(
  .x, .y,
  SpeciesColumn = "SamplingNonZeroSpecies",
  AbundanceColumn = "SamplingNonZeroAbundances",
  Method = "Jaccard",
  PresenceAbsence = TRUE
) {
  stopifnot(length(.x$Patch) == length(unique(.x$Patch)))

  .x$Patch <- as.character(.x$Patch)

  flagMeaningless <- .x$Patch[.x[[AbundanceColumn]] == ""]

  uniqueSpecies <- sort(unique(unlist(strsplit(
    .x[[SpeciesColumn]],
    split = ", "))))

  comdatmat <- matrix(0,
                      nrow = length(.x$Patch),
                      ncol = length(uniqueSpecies))
  colnames(comdatmat) <- uniqueSpecies
  rownames(comdatmat) <- .x$Patch

  for(i in seq_along(.x$Patch)) {
    # print(strsplit(.x$SamplingNonZeroSpecies[i], split = ", ")[[1]])
    # print(strsplit(.x$SamplingNonZeroAbundances[i], split = ", ")[[1]])
    if (SpeciesColumn == AbundanceColumn) {
      # Abundance = Species => List of IDs with repetition for abundance.
      vals <- table(strsplit(.x[[AbundanceColumn]][i], split = ", ")[[1]])
      indices <- names(vals)
    } else {
      # Else, abundance in same order as species detected.
      vals <- as.numeric(strsplit(.x[[AbundanceColumn]][i], split = ", ")[[1]])
      indices <- strsplit(.x[[SpeciesColumn]][i], split = ", ")[[1]]
    }

    comdatmat[.x$Patch[i], indices] <- vals

  }

  if (PresenceAbsence) comdatmat <- comdatmat > 0

  # Short for Jaccard, our default.
  Jacs <- vegan::vegdist(method = tolower(Method), x = comdatmat)

  data.frame(
    Beta = Jacs[1:length(Jacs)],
    Patch1 = rep(attr(Jacs, "Labels"), (length(attr(Jacs, "Labels"))-1):0),
    Patch2 = attr(Jacs, "Labels")[
      sequence(from = seq_along(attr(Jacs, "Labels"))[-1],
               nvec = (length(attr(Jacs, "Labels")) - 1):1)
      ]
  ) %>% dplyr::mutate(
    Meaningless = Patch1 %in% flagMeaningless | Patch2 %in% flagMeaningless
  )

}
