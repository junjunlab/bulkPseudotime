#' bulkPseudotime class
#'
#' @slot sample_dist_plot list.
#' @slot sample_pca_plot list.
#' @slot pseudotime_matrix data.frame.
#' @slot pseudotime_matrix_long data.frame.
#' @slot pseudotime_anno data.frame.
#' @slot ordered_matrix_list list.
#' @slot heatmap_list list.
#' @slot heatmap_list_group list.
#'
#' @export
bulkPseudotimeClass <- setClass("bulkPseudotimeClass",
                                slots = list(sample_dist_plot = "list",
                                             sample_pca_plot = "list",
                                             pseudotime_matrix = "data.frame",
                                             pseudotime_matrix_long = "data.frame",
                                             pseudotime_anno = "data.frame",
                                             ordered_matrix_list = "list",
                                             heatmap_list = "list",
                                             heatmap_list_group = "list"))
