
#set type to '250k' for blanks+250K analysis
bar_plot_update_manual <- function(a, met, Title, x, y, axis.text.x, scales, type = NULL, num_cond=NULL,index)
{
 if (num_cond > 11)
  {
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    col<-c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","sienna1","darkgreen","khaki1","lightpink1","navy","olivedrab1",
           "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
           "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")
    col[index[[1]]]<-"yellow1"
    col[index[[2]]]<-"grey45"
    #col=sample(color, 30)  ##here 12 will pull twelve random colors and store them in 'col'
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
      geom_errorbar(aes(ymin=Norm_Av, ymax=Norm_Av+Norm_Std), position=position_dodge(0.9), width=.2)+
      #geom_text(vjust=-0.1065, color='darkblue', fontface='bold')+
      scale_fill_manual(values = col)
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
        geom_errorbar(aes(ymin=Norm_Av, ymax=Norm_Av+Norm_Std), position=position_dodge(0.9), width=.2)+
        #geom_text(vjust=-0.1065, color='darkblue', fontface='bold')+
       # scale_fill_brewer(palette = "Set1")
        scale_fill_brewer(palette = "Spectral")
    }
  else
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
      geom_errorbar(aes(ymin=Norm_Av, ymax=Norm_Av+Norm_Std), position=position_dodge(0.9), width=.2)+
     # geom_text(vjust=-0.1065, color='darkblue', fontface='bold')+
      #scale_fill_brewer(palette = "Spectral")
      scale_fill_brewer(palette = "Set1")
  }
      
} 

