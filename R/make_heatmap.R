#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param A dataframe containing the relative amounts data, a data frame with information on the samples, and the heatmap color scheme.
#'
#' @return An annotated heatmap in a pdf file.
#'
#' @examples make_heatmap(matrix = RelA, samples = samples, heat.color = normal)
#'
#' @export
#'



make_heatmap <- function(matrix, samples=samples, heat.color=normal){
if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Exp') {
  ext = 'Relative Amounts'
} else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='MID') {
  ext = 'MIDs'
} else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='FC') {
  ext = 'Fractional Contribution'
} else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Labeled') {
  ext = 'Percent Labeled'
}

if (exists('samples')==F) samples <- info
ann <- select(samples, Condition, Cell.Number) %>%
  as.data.frame()
rownames(ann) <- colnames(matrix)

if (exists('Norv')==T) {
  if(all(is.na(Norv)))
    Norv <- rep(1, length(Norv))
  ann$Norvaline <- Norv
} else ann$Norvaline <- 1

ann_colors = list(
  Condition=colors[1:length(unique(gsub('_(.)*','',colnames(matrix))))],
  Norvaline = c("white", "purple"),
  Cell.Number = c("white", "green")
)
names(ann_colors[['Condition']]) <- unique(gsub('_(.)*','',colnames(matrix)))

matrix[is.na(matrix)] <- 0
heatmap.title=paste(Title, '-Heatmap-',ext,'.pdf', sep='')
pheatmap::pheatmap(matrix, cluster_row=T, cluster_col=T,
                   clustering_distance_rows='correlation',
                   clustering_distance_cols='correlation',
                   color = colorRampPalette(heat.color)(100),
                   border_color="black", scale="row",
                   cellwidth = 20, cellheight = 10,
                   annotation=ann, annotation_colors = ann_colors,
                   show_colnames = F, main=paste(Title,ext,sep='-'),
                   filename=heatmap.title)
dev.off()
}
