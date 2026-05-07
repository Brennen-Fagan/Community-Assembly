
# Editted Gemini variation of existing script
## --- R Project Clean-up & Deprecation Audit (Cross-Platform Optimized) ---
# Check all functions used, sort them by whether they are internal to this pkg.
# Goal: identify code no longer used in the production of the manuscripts.
library(here)
library(NCmisc)
library(lintr)
library(RMTRCode2)

# 1. DEFINE TARGETS
# We anchor all paths to the project root to ensure this script runs
# identically across OS types.
all_our_functions <- ls("package:RMTRCode2")

r_dirs <- c(
  here("R"),
  here("experiments", "FaganEtAl2023", "R"),
  here("experiments", "FaganEtAl2026", "R")
)
exp_dirs <- c(
  here("experiments", "FaganEtAl2023", "Figure3-ExampleOutcomes"),
  here("experiments", "FaganEtAl2023", "Figure4-MetasimulationStudy"),
  here("experiments", "FaganEtAl2026"),
  here("experiments", "FaganEtAl2026", "Tests"),
  here("experiments", "Tests")
)

# 2. FUNCTION
# Standard search for filtering both in and between calls.
get_used_functions <- function(dirs) {
  files <- list.files(
    path = dirs,
    pattern = "(?i)\\.R$", # Case insensitive
    recursive = TRUE,
    full.names = TRUE
  )

  # Some local folders need to not be included
  files <- files[!grepl("/(Deprecated|Extraneous)/", files)]

  if(length(files) == 0) return(character(0))

  usage_list <- lapply(files, NCmisc::list.functions.in.file)

  # Extract names belonging to this project's namespace.
  # character(0) captures internal calls (no package prefix).
  unlist(lapply(usage_list, function(x) {
    my_namespaces <- c("character(0)", "package:RMTRCode2")
    return(unlist(x[names(x) %in% my_namespaces]))
  }))
}

# 3. IDENTIFY ACTIVE FILES
used_in_experiments <- get_used_functions(exp_dirs)
used_in_r_files <- get_used_functions(r_dirs)
all_used <- unique(c(used_in_experiments, used_in_r_files))

# 4. FILTER EXTERNAL DEPENDENCIES
not_mine <- c(
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

actually_used <- all_used[!(all_used %in% not_mine)]

# 5. RESULTS
# setdiff for functions in the package but not called in analyses.
check_functions <- setdiff(all_our_functions, actually_used)

# 6. LINTING
r_files_to_lint <- list.files(path = r_dirs, pattern = "(?i)\\.R$", full.names = TRUE)

check_lint <- do.call(rbind, lapply(r_files_to_lint, function(f) {
  as.data.frame(lintr::lint(f, linters = lintr::object_usage_linter()))
}))
