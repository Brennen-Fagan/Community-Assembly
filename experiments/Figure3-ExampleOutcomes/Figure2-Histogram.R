# Histogram, saved at 260 width by 320 height
# follows on from loading a PoolMats.RData object.
# Intended for usage with 
#    Data_2022-09-16/MNA-ExampleOutcome-PoolMats-Env10.RData

tosave <- ggplot2::ggplot(
  Pool, 
  ggplot2::aes(
    x = Size, 
    fill = Type
  )
) + ggplot2::geom_histogram(
  bins = 25, color = "black"
) + ggplot2::theme_bw(
  base_size = 22, base_family = "Arial"
) + ggplot2::coord_flip(
) + ggplot2::scale_fill_discrete(
  type = c("green4", "yellow")
) + ggplot2::ylab(
  "   Count"
) + ggplot2::theme(
  legend.position = -c(0.25, 0.15), 
  legend.title = ggplot2::element_blank(), 
  legend.box.background = ggplot2::element_blank(), 
  legend.background = ggplot2::element_blank(), 
  legend.spacing.x = ggplot2::unit(0.01, 'cm')
)

ggplot2::ggsave(
    tosave,
    "Pool-260x320.png",
    units = "px",
    width = 260, height = 320, dpi = 300
)