# Check all functions used, sort them by whether they are internal to this pkg.

library(NCmisc)

# Assumed this dir is root package directory.
thisdir <- '.'
targdirs <- file.path(thisdir, "experiments", c(
  "Figure3-ExampleOutcomes", "Figure4-MetasimulationStudy", "Tests"
))

files <- dir(
  path = targdirs, pattern = "[.]R$",
  full.names = TRUE, include.dirs = TRUE,
  no.. = TRUE, ignore.case = TRUE,
  recursive = TRUE # should be excessive
)

functionsused <- lapply(files, NCmisc::list.functions.in.file)
functionsused_nots <- lapply(functionsused, function(x, keep) {
  x[names(x) %in% keep]
}, keep = "character(0)")
functionsused_nots <- sort(unique(unlist(functionsused_nots)))

notmine <- c(
  '.', 'across', 'aes', 'aes_string', 'all_of', 'annotate', 'anti_join',
  'any_of', 'arrange', 'arrangeGrob', 'bandSparse', 'bdiag', 'bind_rows',
  'case_when', 'clusterEvalQ', 'clusterExport', 'cmpfun', 'coord_cartesian',
  'coord_flip', 'Diagonal', 'distinct', 'drop0', 'element_blank',
  'element_rect', 'element_text', 'expand_grid', 'expand_limits', 'facet_grid',
  'facet_wrap', 'file_path_sans_ext', 'foreach', 'full_join',"geom_abline",
  "geom_bin2d", "geom_boxplot", "geom_density_2d", "geom_histogram",
  "geom_hline", "geom_line", "geom_path", "geom_point", "geom_raster",
  "geom_ribbon", "geom_text", "geom_vline", "ggplot", "ggplot_build",
  "ggplot_gtable", "ggsave", "ggtitle", "group_by", "group_split", 'iter',
  'labs', 'left_join', 'lsoda', 'lsode', 'makeCluster', 'Matrix', 'mutate', 'n',
  'ode', 'percent_format', 'pivot_longer', 'pivot_wider', 'plot_layout',
  'position_dodge', 'pull', 'registerDoParallel', 'rename',
  'rownames_to_column', 'runsteady', 'scale_alpha', 'scale_color_continuous',
  'scale_color_discrete', "scale_color_manual", "scale_color_viridis_c",
  "scale_fill_discrete", "scale_fill_manual", "scale_fill_viridis_c",
  "scale_x_continuous", "scale_x_discrete", "scale_y_continuous", "select",
  "semi_join", "separate", "slice_max", 'sparseMatrix', 'starts_with',
  'stat_summary', 'steady', 'stopCluster', 'summarise', 'sym', 'theme',
  'theme_bw', 'theme_void', 'ungroup', 'unit', 'unite', 'vegdist', 'wrap_plots',
  'xlab', 'ylab'

)

functionsused_nots <- functionsused_nots[
  !(functionsused_nots %in% notmine)
  ]
