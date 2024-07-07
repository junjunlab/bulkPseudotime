globalVariables(c('order', 'timePointCol', 'timelineCol', 'ht_col', 'heatmap_params',
                  'markGenes', 'fontsize', '...','genes', 'ncol',"gene","pseudotime"))

#' bulkPseudotime function
#'
#' @param expMat Expression matrix(tpm/fpkm data) for bulk RNA-seq.
#' @param timePointCol the color for timePointCol.
#' @param timelineCol the color for timelineCol.
#'
#' @return bulkPseudotime object
#'
#' @importFrom stats dist loess
#' @import ggplot2
#' @import utils
#' @importFrom scales rescale
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom dplyr left_join mutate filter
#'
#' @export
bulkPseudotime <- function(expMat = NULL,
                           timePointCol = NULL,
                           timelineCol = NULL){
  # ==============================================================================
  # 1.get the timeline
  # ==============================================================================
  pca <- FactoMineR::PCA(t(expMat),scale.unit = T,ncp = 5,graph = F)

  # get distance for samples
  res <- pca$ind$coord
  res <- res[,1:2]

  dis <- stats::dist(res)
  dis <- as.matrix(dis)

  # sample_dist_plot
  sample_dist_plot <- Heatmap(matrix = dis,cluster_rows = F,cluster_columns = F)

  res_p <- data.frame(res)
  res_p$sample <- rownames(res_p)

  # calculate distance
  raw.timeline <- c(0,lapply(1:ncol(expMat),function(x){
    point <- 0 + dis[x - 1,x]
  }) %>% unlist() %>% cumsum())

  # rescale to c(0,10)
  new.timeline <- scales::rescale(raw.timeline,to = c(0,10))

  # ==============================================================================
  # 2.produce gene expression according timelines
  # ==============================================================================
  # do zscore
  tpm.z <- as.data.frame(t(apply(expMat,1,scale)))
  names(tpm.z) <- names(expMat)

  # produce continues point
  grid <- data.frame(time = seq(0,10,length.out = 100*ncol(expMat)))

  # make time point annotations according to colnames
  col_name <- colnames(expMat)

  cut_time <- cut(x = grid$time,
                  breaks = new.timeline,
                  labels = col_name[2:length(col_name)]) %>%
    as.character()

  cut_time[1] <- col_name[1]

  grid_anno <- data.frame(time = grid$time,time_anno = cut_time)


  # prediction method
  pseudotime.model.fun <- function(value){
    time <- new.timeline
    data <- tibble::tibble(value = value,time = time)
    model <- stats::loess(value ~ time,data)
    predict <- grid %>% modelr::add_predictions(model)
    return(predict)
  }

  # extarct value from prediction function
  get_predicted_exp <- function(mat = NULL,timeline = NULL){
    res <- apply(mat,1,pseudotime.model.fun)

    # x = 1
    results <- lapply(seq_along(res),function(x){
      tmp <- res[[x]]$pred
      return(tmp)
    }) %>% do.call("rbind",.) %>% data.frame()

    rownames(results) <- row.names(mat)
    colnames(results) <- timeline$time
    return(results)
  }

  results <- get_predicted_exp(mat = tpm.z,timeline = grid)
  results_raw <- get_predicted_exp(mat = expMat,timeline = grid)

  # ============================================================================
  # timeline annotation
  if(is.null(timePointCol)){
    time_col <- circlize::rand_color(n = length(unique(grid_anno$time_anno)))
  }else{
    time_col <- timePointCol
  }

  names(time_col) <- unique(grid_anno$time_anno)

  if(is.null(timelineCol)){
    col_fun = circlize::colorRamp2(c(0,10), c("white", "#993399"))
  }else{
    col_fun = circlize::colorRamp2(seq(0,10,length = length(timelineCol)), timelineCol)
  }

  time_anno <- HeatmapAnnotation(timePoint = grid_anno$time_anno,
                                 timeline = grid_anno$time,
                                 col = list(timePoint = time_col,
                                            timeline = col_fun),
                                 border = T)

  # ==============================================================================
  # 3.order genes
  # ==============================================================================
  # get Dim1 and Dim2 for genes
  gene.pca <- FactoMineR::PCA(results,scale.unit = T,ncp = 5,graph = F)
  res2 <- gene.pca$ind$coord
  res2 <- res2[,1:2]

  results2 <- cbind(results,res2)

  gene.pca2 <- FactoMineR::PCA(t(results),scale.unit = T,ncp = 5,graph = F)
  res3 <- gene.pca2$ind$coord
  res3 <- data.frame(res3[,1:2])
  res3$time_anno <- grid_anno$time_anno

  results2_raw <- cbind(results_raw,res2) %>%
    tibble::rownames_to_column(var = "gene")

  # long format
  results2_raw_long <- reshape2::melt(results2_raw,
                                      variable.name = "pseudotime",
                                      value.name = "exp",
                                      id.vars = c("Dim.1","Dim.2","gene")) %>%
    dplyr::mutate(pseudotime = as.numeric(as.character(pseudotime))) %>%
    dplyr::left_join(y = grid_anno,by = c("pseudotime" = "time"))

  # sample_pca_plot
  sample_pca_plot <-
    ggplot() +
    geom_point(data = res3,aes(x = Dim.1,y = Dim.2,color = time_anno)) +
    geom_path(data = res_p,aes(x = Dim.1,y = Dim.2),linewidth = 1) +
    geom_point(data = res_p,aes(x = Dim.1,y = Dim.2,color = sample),
               size = 5) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    scale_color_brewer(palette = "Set3")


  # calculate atan2 value
  results2$atan2.1 <- atan2(results2$Dim.1,results2$Dim.2)
  results2$atan2.2 <- atan2(results2$Dim.1,-results2$Dim.2)
  results2$atan2.3 <- atan2(-results2$Dim.1,results2$Dim.2)
  results2$atan2.4 <- atan2(-results2$Dim.1,-results2$Dim.2)

  # order
  order1 <- dplyr::arrange(results2,results2$atan2.1)
  order2 <- dplyr::arrange(results2,results2$atan2.2)
  order3 <- dplyr::arrange(results2,results2$atan2.3)
  order4 <- dplyr::arrange(results2,results2$atan2.4)

  order_list <- list(order1,order2,order3,order4)
  names(order_list) <- paste0("Order",1:4)

  # loop for puduce heatmap
  lapply(1:4,function(x){
    ht <- Heatmap(matrix = order_list[[x]][,1:(100*ncol(expMat))],
                  name = "zscore",
                  border = T,
                  top_annotation = time_anno,
                  cluster_rows = F,cluster_columns = F,
                  show_row_names = F,show_column_names = F,
                  column_title = paste0("Order",x))
    return(ht)
  }) -> ht_list

  names(ht_list) <- paste0("Order",1:4)

  # combine
  ht_listcb <- Reduce("+",ht_list)

  # ============================================================================
  # return data
  # ============================================================================
  # output <- list(sample_dist_plot = sample_dist_plot,
  #                sample_pca_plot = sample_pca_plot,
  #                pseudotime_matrix = results2,
  #                pseudotime_matrix_long = results2_raw_long,
  #                pseudotime_anno = grid_anno,
  #                ordered_matrix_list = order_list,
  #                heatmap_list = ht_list,
  #                heatmap_list_group = ht_listcb)

  # ============================================================================
  # create homerResult object
  # ============================================================================
  res <- methods::new("bulkPseudotime",
                      sample_dist_plot = list(sample_dist_plot),
                      sample_pca_plot = list(sample_pca_plot),
                      pseudotime_matrix = results2,
                      pseudotime_matrix_long = results2_raw_long,
                      pseudotime_anno = grid_anno,
                      ordered_matrix_list = list(order_list),
                      heatmap_list = list(ht_list),
                      heatmap_list_group = list(ht_listcb))

  return(res)
}
