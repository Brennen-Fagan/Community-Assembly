extractAttributes <- function(Diversity, idNums) {
  data.frame(
    Pool = toString(Diversity$PoolMod),
    Noise = toString(Diversity$NoiseMod),
    Neutral = toString(Diversity$NeutralMod),
    Space = toString(Diversity$SpaceMod),
    Set = idNums[1],
    CaseNumber = idNums[2],
    History = if(is.na(idNums[3])) {
      # All Histories
      1:10
    } else idNums[3],
    Part = if(is.na(idNums[4])) {
      # All Parts
      1
    } else idNums[4],
    stringsAsFactors = FALSE
  )
}