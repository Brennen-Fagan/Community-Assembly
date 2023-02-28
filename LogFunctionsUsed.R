# Check all functions used, sort them by whether they are internal to this pkg.

library(NCmisc)

library(RMTRCode2)
allOurFunctions <- ls("package:RMTRCode2")

# Assumed this dir is root package directory.
thisdir <- '.'
targdirs <- file.path(thisdir, "experiments", c(
  "Figure3-ExampleOutcomes", "Figure4-MetasimulationStudy", "Tests"
))
Rdir <- file.path(thisdir, "R")

# Get files, look for functions used that are mine.
targfiles <- dir(
  path = targdirs, pattern = "[.]R$",
  full.names = TRUE, include.dirs = TRUE,
  no.. = TRUE, ignore.case = TRUE,
  recursive = TRUE # should be excessive
)

functionsused <- lapply(targfiles, NCmisc::list.functions.in.file)
names(functionsused) <- targfiles
functionsused_nots <- lapply(functionsused, function(x, keep) {
  x[names(x) %in% keep]
}, keep = c("character(0)", "package:RMTRCode2"))
functionsused_nots <- sort(unique(unlist(functionsused_nots)))

notmine <- c(
  '.', 'across', 'aes', 'aes_string', 'all_of', 'annotate', 'anti_join',
  'any_of', 'arrange', 'arrangeGrob', "as_data_frame", 'bandSparse', 'bdiag',
  'bind_rows', 'case_when', 'clusterEvalQ', 'clusterExport', 'cmpfun',
  'coord_cartesian', 'coord_flip', 'delete_edges', 'Diagonal', 'distinct',
  'drop0', 'edge_attr', 'element_blank', 'element_rect', 'element_text',
  'expand_grid', 'expand_limits', 'facet_grid', 'facet_wrap',
  'file_path_sans_ext', "FlowBasedTrophicLevel", 'foreach', 'full_join',
  "geom_abline", "geom_bin2d", "geom_boxplot", "geom_density_2d",
  "geom_histogram", "geom_hline", "geom_label", "geom_line", "geom_path",
  "geom_point", "geom_raster", "geom_ribbon", "geom_text", "geom_vline",
  "ggplot", "ggplot_build", "ggplot_gtable", "ggsave", "ggtitle",
  "graph_from_adjacency_matrix", "graph_from_data_frame", "group_by",
  "group_split", "guides", 'iter', 'labs', "layout.circle", 'left_join',
  'lsoda', 'lsode', 'makeCluster', 'Matrix', 'mutate', 'n', 'ode',
  'percent_format', 'pivot_longer', 'pivot_wider', 'plot_layout',
  "PlotWebByLevel", 'position_dodge', 'pull', 'registerDoParallel', 'rename',
  'rownames_to_column', 'runsteady', 'scale_alpha', 'scale_color_continuous',
  'scale_color_discrete', "scale_color_manual", "scale_color_viridis_c",
  "scale_fill_discrete", "scale_fill_manual", "scale_fill_viridis_c",
  "scale_x_continuous", "scale_x_discrete", "scale_y_continuous", "select",
  "semi_join", "separate", "set.vertex.attribute", "slice_max", 'sparseMatrix',
  'starts_with', 'stat_summary', 'steady', "stode", 'stopCluster', 'summarise',
  'sym', 'theme', 'theme_bw', 'theme_void', "TrophicLevels", 'ungroup', 'unit',
  'unite', 'vegdist', "vertex_attr", 'wrap_plots', 'xlab', 'ylab'

)

functionsused_nots <- functionsused_nots[
  !(functionsused_nots %in% notmine)
  ]

# Repeat for R files.
Rfiles <- dir(
  path = Rdir, pattern = "[.]R$",
  full.names = TRUE, include.dirs = TRUE,
  no.. = TRUE, ignore.case = TRUE,
  recursive = TRUE # should be excessive
)

functionsusedR <- lapply(Rfiles, NCmisc::list.functions.in.file)
names(functionsusedR) <- Rfiles

functionsusedR_nots <- lapply(functionsusedR, function(x, keep) {
  x[names(x) %in% keep]
}, keep = c("character(0)", "package:RMTRCode2"))
functionsusedR_nots <- sort(unique(unlist(functionsusedR_nots)))
functionsusedR_nots <- functionsusedR_nots[
  !(functionsusedR_nots %in% notmine)
  ]

# Compare the two to identify functions to potentially deprecate and move to
# the "Extraneous" directory.
possiblyDeprecated <- allOurFunctions[
  !(allOurFunctions %in% c(functionsused_nots, functionsusedR_nots))
  ]


# The rest of the inspection needs to be by hand.
