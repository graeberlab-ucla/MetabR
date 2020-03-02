#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param A ggplot object, the metabolite pathway of interest, and various other graphical parameters.
#'
#' @return Relative amounts plots for the metabolites in the specified pathway (met).
#'
#' @examples bar_plot_update_manual(a, met, x, y, axis.text.x scales, type, num_cond)
#'
#' @export
#'

bar_plot_update_manual <- function(a, met, Title, x, y, axis.text.x, scales, type = NULL, num_cond=NULL,index = NULL)
{
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  col<-c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","sienna1","darkgreen","lightpink1","navy","olivedrab1",
         "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
         "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")

  if(!is.null(index))
  {
    j  <- 1
    k <- 1
    extra_qc <- c("peachpuff1", "seashell1", "wheat2", "snow1")
    res <- vector()
    for( i in 1:num_cond)
    {
      if(i %in% index[[1]])
        res <- c(res, "yellow1")
      else if(i %in% index[[2]])
        res <- c(res, "grey45")
      else if(i %in% index[[3]])
      {
        res <- c(res, extra_qc[k])
        k <- k + 1
      }
      else
      {
        res <- c(res, col[j])
        j <- j + 1
      }
    }
    col <- res
  }


  if (num_cond > 11)
  {
    a + geom_bar(position="dodge", stat="identity", width=0.9) +
      geom_bar(position="dodge", stat="identity", colour="black", width=0.9) +
      facet_wrap( ~ Name, scales=scales) +
      theme_bw() +
      labs(x=x, y=y, title=Title, fill=element_blank()) +
      theme(
        plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
        axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
        axis.text=element_text(size=11, face="bold"),
        #axis.text.x=element_blank(),
        axis.text.x=axis.text.x,
        legend.title=element_text(face="bold", size=12),
        legend.text=element_text(face="bold",size=12),                  #sets legend text
        strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
        panel.grid.major=element_blank()) +
      geom_errorbar(aes(ymin=RelAmounts_Ave, ymax=RelAmounts_Ave+RelAmounts_Std), position=position_dodge(0.9), width=.2)+
      scale_fill_manual(values = col) + geom_text(aes(label=under_50_percent), position=position_dodge(width=0.9), vjust=-0.25)+
      geom_bar(position="dodge", stat="identity", colour="black", width=0.9)
  }
  else if (type == '250k')
  {
    a + geom_bar(position="dodge", stat="identity", width=0.9) +
      geom_bar(position="dodge", stat="identity", colour="black", width=0.9) +
      facet_wrap( ~ Name, scales=scales) +
      theme_bw() +
      labs(x=x, y=y, title=Title, fill=element_blank()) +
      theme(
        plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
        axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
        axis.text=element_text(size=11, face="bold"),
        #axis.text.x=element_blank(),
        axis.text.x=axis.text.x,
        legend.title=element_text(face="bold", size=12),
        legend.text=element_text(face="bold",size=12),                  #sets legend text
        strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
        panel.grid.major=element_blank()) +
      geom_errorbar(aes(ymin=RelAmounts_Ave, ymax=RelAmounts_Ave+RelAmounts_Std), position=position_dodge(0.9), width=.2)+
      #scale_fill_brewer(palette = "Spectral")
      scale_fill_manual(values = col)
  }
  else
  {
    a + geom_bar(position="dodge", stat="identity", width=0.9) +
      geom_bar(aes(linetype=under_50_percent,colour=ifelse(under_50_percent=="X","snow3","black"), size=under_50_percent),position="dodge", stat="identity", width=0.9) +
      guides(linetype=FALSE)+scale_size_manual(values=c(0.3,1.3),guide=FALSE)+facet_wrap( ~ Name, scales=scales) +
      theme_bw() +scale_linetype_manual(values=c("solid","58"))+
      labs(x=x, y=y, title=Title, fill=element_blank()) +
      theme(
        plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
        axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
        axis.text=element_text(size=11, face="bold"),
        #axis.text.x=element_blank(),
        axis.text.x=axis.text.x,
        legend.title=element_text(face="bold", size=12),
        legend.text=element_text(face="bold",size=12),                  #sets legend text
        strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
        panel.grid.major=element_blank()) +
      geom_errorbar(aes(ymin=RelAmounts_Ave, ymax=RelAmounts_Ave+RelAmounts_Std), position=position_dodge(0.9), width=.2)+
      #scale_fill_brewer(palette = "Set1")
      scale_fill_manual(values = col)
  }

}