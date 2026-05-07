# Base:
# Algorithmic: first index = 100%, second index = 50%
#              (0) => 1,0,0, (0.5) => 0,1,0, (1) => 0,0,1
#              (0.25) => 0.5,0.5,0, (0.75) => 0,0.5,0.5

# Modification:
# Trying to move colors out from M and towards C and Y.
colorPaletteAlg <- function(intervention) {
  intervention <- as.numeric(strsplit(
    gsub(intervention, pattern = "[(|)]", replacement = ""),
    split = "->")[[1]])
  x <- intervention[1]
  y <- if(length(intervention) == 2) {
    intervention[2]
  } else {
    intervention[1]
  }

  DescTools::CmykToRgb(
    min(1.25*(max(0, (0.5-x)/0.5) + 0.5*max(0, (0.5-y)/0.5)), 1),
    min(max(0, (0.5 - abs(x - 0.5))/0.5)
        + 0.5*max(0, (0.5 - abs(y - 0.5))/0.5), 1),
    min(1.25*(max(0, (x-0.5)/0.5)+ 0.5*max(0, (y-0.5)/0.5)), 1),
    0.25
  )
}
