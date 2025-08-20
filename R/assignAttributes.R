assignAttributes <- function(i, d, a) {
  d %>% dplyr::bind_rows(
  ) %>% dplyr::mutate(
    Set = a$Set[i],
    Number = a$CaseNumber[i],
    History = a$History[i],
    Pool = a$Pool[i],
    Noise = a$Noise[i],
    Neutral = a$Neutral[i],
    Space = a$Space[i]
  )
}