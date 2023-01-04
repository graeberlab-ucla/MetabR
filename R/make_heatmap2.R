#' Make heatmap
#'
#' @param matrix Data matrix
#' @param samples Sample info sheet data
#' @param heat.color Color
#' @param cluster_samples Logical; T to cluster
#' @param width Plot width
#' @param height Plot height
#'
#' @return Heatmap
#' @export
#'
make_heatmap2 <- function (matrix, samples = samples, heat.color = normal, cluster_samples = TRUE, width = NA, height = NA, ext=NA)
{
  if(is.na(ext)){
  if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "Exp") {
    ext = "Relative Amounts"
  } else if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "MID") {
    ext = "MIDs"
  } else if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "FC") {
    ext = "Fractional Contribution"
  } else if (gsub("(.)*_|[0-9]", "", colnames(matrix))[1] == "Labeled") {
    ext = "Percent Labeled"
  }
  }
  if (exists("samples") == F)
    samples <- info
  ann <- select(samples, Condition, Cell.Number) %>% as.data.frame()
  rownames(ann) <- colnames(matrix)
  #ann$Cell.Number <- 1
  if(all(is.na(Norv))| length(Norv)==0) {
    ann_colors = list(
      Condition=colors[1:length(unique(gsub('_(.)*','',colnames(matrix))))],
      Cell.Number = c("white", "green"))
  }else{
    ann_colors = list(
      Condition=colors[1:length(unique(gsub('_(.)*','',colnames(matrix))))],
      Norvaline = c("white", "purple"),
      Cell.Number = c("white", "green"))
    ann$Norvaline <- Norv
  }
  names(ann_colors[["Condition"]]) <- unique(gsub("_(.)*",
                                                  "", colnames(matrix)))
  if (all(is.na(ann$Cell.Number))) { # Remove Cell.Number annotation if not provided
   ann$Cell.Number <- NULL
   ann_colors[["Cell.Number"]] <- NULL
  }

  matrix[is.na(matrix)] <- 0
  heatmap.title = paste(Title, "-Heatmap-", ext, ".pdf", sep = "")
  if (cluster_samples == F)
    heatmap.title = paste(Title, "-Heatmap-", ext, "-Unclustered",
                          ".pdf", sep = "")
  if (is.na(width) & is.na(height)) {
    pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = cluster_samples,
                       clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                       color = colorRampPalette(heat.color)(100), border_color = "black",
                       scale = "row", cellwidth = 20, cellheight = 10,
                       annotation = ann, annotation_colors = ann_colors,
                       show_colnames = F, main = paste(Title, ext, sep = "-"),
                       filename = heatmap.title)
  }else {
    ### Added code to remove rows with zero variance
    zero_variance <- apply(matrix, 1, var)
    zero_variance_vec <- zero_variance[zero_variance == 0]
    matrix_filt <- matrix[!rownames(matrix) %in% names(zero_variance_vec),]
    matrix_remove <- matrix[rownames(matrix) %in% names(zero_variance_vec),]
    ### Changed matrix to matrix_filt
    pheatmap::pheatmap(matrix_filt, cluster_row = T, cluster_col = cluster_samples
                       ,clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                       color = colorRampPalette(heat.color)(100), border_color = "black",
                       scale = "row", cellwidth = 20, cellheight = 10,
                       annotation = ann, annotation_colors = ann_colors,
                       show_colnames = F, main = paste(Title, ext, sep = "-"),
                       filename = heatmap.title, width = width, height = height)

    #pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = cluster_samples,
    #clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
    #color = colorRampPalette(heat.color)(100), border_color = "black",
    #scale = "row", cellwidth = 20, cellheight = 10,
    #annotation = ann, annotation_colors = ann_colors,
    #show_colnames = F, main = paste(Title, ext, sep = "-"),
    #filename = heatmap.title, width = width, height = height)
  }
  #dev.off()
}
