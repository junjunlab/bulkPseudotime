# ==============================================================================
# pseudotime_heatmap function for plot
# ==============================================================================

#' pseudotime_heatmap
#'
#' @param object bulkPseudotime object.
#' @param order which order heatmap to be choosed to plot.
#' @param timePointCol color for timePointCol.
#' @param timelineCol color for timelineCol.
#' @param ht_col color for heatmp.
#' @param heatmap_params more parameters for Heatmap function.
#' @param markGenes the gene names to be marked in heatmap.
#' @param fontsize fontsize for mark genes.
#' @param ... other args will be ignored.
#'
#' @export
setGeneric("pseudotime_heatmap",function(object,order = NULL,
                                         timePointCol = NULL,
                                         timelineCol = NULL,
                                         ht_col = NULL,
                                         heatmap_params = list(),
                                         markGenes = NULL,
                                         fontsize = 12,...) standardGeneric("pseudotime_heatmap"))


#' method for pseudotime_heatmap
#'
#' @param object bulkPseudotime.
#'
#' @return Heatmap object
#' @import methods
#' @export
setMethod("pseudotime_heatmap",
          signature(object = "bulkPseudotime"),
          function(object = NULL,
                   order = NULL,
                   timePointCol = NULL,
                   timelineCol = NULL,
                   ht_col = NULL,
                   heatmap_params = list(),
                   markGenes = NULL,
                   fontsize = 12,
                   ...){
            # ============================================================================
            # timeline annotation
            # ============================================================================
            grid_anno <- object@pseudotime_anno
            if(is.null(timePointCol)){
              time_col <- circlize::rand_color(n = length(unique(grid_anno$time_anno)))
            }else{
              time_col <- timePointCol
            }

            names(time_col) <- unique(grid_anno$time_anno)

            if(is.null(timelineCol)){
              col_fun = circlize::colorRamp2(c(0,10), c("white", "#0066CC"))
            }else{
              col_fun = circlize::colorRamp2(seq(0,10,length = length(timelineCol)), timelineCol)
            }

            if(is.null(ht_col)){
              ht_col <- circlize::colorRamp2(c(-2, 0, 2), c("#006633", "white", "#990066"))
            }else{
              ht_col <- circlize::colorRamp2(c(-2, 0, 2), ht_col)
            }

            time_anno <- HeatmapAnnotation(timePoint = grid_anno$time_anno,
                                           timeline = grid_anno$time,
                                           col = list(timePoint = time_col,
                                                      timeline = col_fun),
                                           border = T)

            mat <- object$ordered_matrix_list[[order]][,1:length(object@pseudotime_anno$time)]

            # whether mark your genes on plot
            if (!is.null(markGenes)) {
              # all genes
              rowGene <- rownames(mat)

              # tartget gene
              annoGene <- markGenes

              # get target gene index
              index <- match(annoGene, rowGene)

              # some genes annotation
              geneMark <- ComplexHeatmap::rowAnnotation(gene = ComplexHeatmap::anno_mark(at = index,
                                                                                         labels = annoGene,
                                                                                         labels_gp = grid::gpar(
                                                                                           fontface = "italic",
                                                                                           fontsize = fontsize)
              ))

              right_annotation <- geneMark
            } else {
              right_annotation <- NULL
            }

            # ============================================================================
            # heatmap
            # ============================================================================
            # ht <- Heatmap(matrix = object$ordered_matrix_list[[order]][,1:length(object$pseudotime_anno$time)],
            #               name = "zscore",
            #               border = T,
            #               col = ht_col,
            #               top_annotation = time_anno,
            #               cluster_rows = F,cluster_columns = F,
            #               show_row_names = F,show_column_names = F,
            #               column_title = paste0("Order",order))

            ht <- do.call(Heatmap,modifyList(list(matrix = mat,
                                                  name = "zscore",
                                                  border = T,
                                                  col = ht_col,
                                                  top_annotation = time_anno,
                                                  right_annotation = right_annotation,
                                                  cluster_rows = F,cluster_columns = F,
                                                  show_row_names = F,show_column_names = F,
                                                  column_title = paste0("Order",order)),
                                             heatmap_params))

            return(ht)
          }
)


# ==============================================================================
# pseudotime_line function for plot
# ==============================================================================

#' pseudotime_line
#'
#' @param object bulkPseudotime object.
#' @param genes gene names for lineplot.
#' @param ncol columns for facet plot.
#' @param ... other args will be ignored.
#'
#' @return ggplot object.
#' @export
setGeneric("pseudotime_line",function(object,genes = NULL,ncol = NULL,...) standardGeneric("pseudotime_line"))


#' pseudotime_line
#'
#' @param object bulkPseudotime.
#'
#' @return ggplot object
#' @export
setMethod("pseudotime_line",
          signature(object = "bulkPseudotime"),
          function(object = NULL,
                   genes = NULL,ncol = NULL,...){
            df <- object$pseudotime_matrix_long %>%
              dplyr::filter(gene %in% genes)

            # plot
            ggplot(df) +
              geom_line(aes(x = pseudotime,y = exp,color = pseudotime),linewidth = 1) +
              facet_wrap(~gene,scales = "free",ncol = ncol) +
              theme_bw() +
              theme(axis.text = element_text(colour = "black"),
                    strip.text = element_text(face = "bold.italic",size = rel(1.1))) +
              scale_color_viridis_c(option = "turbo") +
              ylab("Gene Relative Expression")
          }
)
