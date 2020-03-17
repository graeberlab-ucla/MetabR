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

plot_istd <- function(Std, info, title,pdf_width=14)
{
  pdf(title, width=pdf_width, height=10)
  std_plot<-ggplot(Std, aes(Sample, log(Value,10)))+
    geom_point(size=3)+
    facet_wrap(~ID, scales='free')+
    theme_bw()+
    theme(text = element_text(face='bold',size=14),
          axis.text.x = element_text(angle=90, vjust=.3, hjust=1))+
    labs(x='',y='Relative response of ISTD (A.U., log10)') + scale_x_discrete(limits = info$Sample)
  print(std_plot)  #print statement is required to plot within if statements, loops, functions
  dev.off()
}
