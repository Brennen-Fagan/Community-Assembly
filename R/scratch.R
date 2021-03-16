# dplyr::mutate(na$Sequence, Type = dplyr::case_when(
#   # If everyone who was here is here and you were here, then you are present.
#   all(lag(Community) %in% Community) & Addition %in% lag(Community) ~ "Present",
#   # If everyone who was here is here and you are not here, then you failed
#   all(lag(Community) %in% Community) & !(Addition %in% Community) ~ "Type 1 (Failure)",
#   # If everyone who was here is here and you are here but were not here, you succeeded
#   # Note: case_when assumes failure of previous conditions, like if-else if.
#   all(lag(Community) %in% Community) & Addition %in% Community ~ "Type 2 (Permanent)",
#   # Otherwise (e.g. invasion caused collapse)
#   TRUE ~ "Type 3 (!Permanent)"
# )
#
# # But for a numerical system, we run into the awkward issue where a species on
# # the decline just happens to disappear when a new invasion happens, rather than
# # being in a true type 3 situation.

# Need to retrieve the length statistics
# na$Sequence$Outcome
