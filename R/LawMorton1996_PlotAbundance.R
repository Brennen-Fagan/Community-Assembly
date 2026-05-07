LawMorton1996_PlotAbundance <- function(
  Abundance,
  Sequence = NULL,
  guides = TRUE
) {
  colnames(Abundance) <- c("Time", format(1:(ncol(Abundance) - 1)))
  
  long <- tidyr::pivot_longer(
    as.data.frame(Abundance),
    2:ncol(Abundance),
    names_to = "Species",
    values_to = "Abundance"
  )
  
  long <- long[!is.na(long$Abundance), ]
  
  thePlot <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = Time,
      y = Abundance,
      color = Species
    )
  ) + ggplot2::geom_line(
  )
  
  if (!is.null(Sequence)) {
    thePlot <- thePlot + ggplot2::geom_vline(
      data = Sequence,
      mapping = ggplot2::aes(xintercept = Events),
      linetype = "dashed",
      color = "black"
    )
    
    timeDiff <- mean(diff(Sequence$Events), na.rm = TRUE)
    
    thePlot <- thePlot + ggplot2::geom_label(
      data = Sequence,
      mapping = ggplot2::aes(
        x = Events + timeDiff/3,
        label = Addition
      ),
      y = max(long$Abundance, na.rm = TRUE) * 0.85,
      color = "black"
    )
  }
  
  if (guides == FALSE) {
    thePlot <- thePlot + ggplot2::guides(color = "none")
  }
  
  return(invisible(thePlot))
}