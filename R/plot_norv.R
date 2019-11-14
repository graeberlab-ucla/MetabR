#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param
#'
#' @return
#'
#' @examples
#'
#' @export
#'

plot_norv(data, info, title)
{
  Norv <- data %>%
    filter(Used_ID=='Norvaline') %>%
    gather(Sample, Norv,-Used_ID, -Iso) %>%
    .$Norv

  pdf(file = title, width=12, height=9, pointsize=12)
  norv_plot<-data %>%
    filter(Used_ID=='Norvaline') %>%
    gather(Sample, Value,-Used_ID, -Iso) %>%
    mutate(Sample=factor(Sample, levels=info$Sample.Name)) %>%
    dplyr::select(Sample, Value) %>%
    ggplot(., aes(Sample, Value)) +
    geom_point(size=3) +
    geom_line(aes(as.integer(Sample), Value), color='blue') +
    labs(x='Sample Name',y='Response',title='Norvaline Response') +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + scale_x_discrete(limit = info$Sample.Name)
  print(norv_plot)
  dev.off()

  return(Norv)
}
