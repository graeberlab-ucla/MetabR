#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#' @param metabolites list of metabolites used to run experiment
#' @param df amounts3 or amounts
#' @param repeats Not currently used
#' @param n n = num_conditions
#' @param type type = "tf" or other
#' @param index index from maven_analysis pre-determined groupings
#' @param title_type "nonpathway1", "nonpathway2", or "nonpathway3", or something else
#' @param only_M0 Not currently being used
#' @param mediumtitle mediumtitle="Plots-Relative-To-Fresh"
#' @param color_lst m_ind<-which(!grepl("unspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia",unique(amounts$Condition)))   amounts3_colors<- colors1[m_ind]
#'
#' @return a plot of the relative amounts or samples relative to unspent of specific metabolites
#' @export
#'
# @examples full_medium_plot_func(unname(unlist(lst)),a, n = num_conditions, type = "tf", title_type = "medium", mediumtitle="Plots-Relative-To-Fresh", color_lst=amounts3_colors)
#'
full_medium_plot_func<-function (metabolites, df, repeats, n, type=NULL, index = NULL, title_type,
                             only_M0 = NA,  mediumtitle=NULL,color_lst=NULL)
{
  if (length(metabolites) > 1) {
    if (title_type == "nonpathway1") {
      ending <- "other metabolites - I"
    }
    else if (title_type == "nonpathway2") {
      ending <- "other metabolites - II"
    }
    else if (title_type == "nonpathway3") {
      ending <- "other metabolites - III"
    }
    else if (title_type == "medium") {
      ending <- "metabolites found in medium"
    }
    else {
      ending <- "select metabolites"
    }
  }
  else ending = ""

  if (sum(grepl("Exp", names(df))) >= 1) {
    if (type == "tf") {
      met = subset(df, Name %in% metabolites)
      stopifnot(length(unique(met$Name)) >= 1)
      met <- mutate(met, Name = paste(Name, Sig, sep = " "),
                    Sig = "")
      Title = paste0("Amounts of ", ending, " (not corrected for blank values)")
      if(!is.null(mediumtitle)) Title<-mediumtitle
      x = ""
      y = "Amounts"
      a <- ggplot(met, aes(Condition, RelAmounts_Ave,
                           group = Condition, fill = Condition, label = Sig))
      axis.text.x = element_blank()

    }
    else {
      met = subset(df, Name %in% metabolites)
      stopifnot(length(unique(met$Name)) >= 1)
      Title = paste0("Amounts of ", ending)
      if(!is.null(mediumtitle)) Title<-mediumtitle
      x = ""
      y = "Amounts"
      a <- ggplot(met, aes(Condition, Av,
                           group = Condition, fill = Condition, label = Name))
      axis.text.x = element_blank()
    }
  }

  scales = "free"
   num_cond=n

   #  if(is.null(color_lst)){
   #    color_lst<-c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","darkgreen","lightpink1","navy","olivedrab1",
   #           "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
   #           "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")
   #  }else{
   #    color_lst<-color_lst
   #  }
   #
   # if(!is.null(index))
   # {
   #   j  <- 1
   #   k <- 1
   #   extra_qc <- c("peachpuff1", "seashell1", "wheat2", "snow1")
   #   res <- vector()
   #   for( i in 1:num_cond)
   #   {
   #     if(i %in% index[[1]])
   #       res <- c(res, "yellow1")
   #     else if(i %in% index[[2]])
   #       res <- c(res, "grey45")
   #     else if (i %in% index[[4]]){
   #       res <- c(res, "darkorange1")
   #     }
   #     else if(i %in% index[[3]])
   #     {
   #       res <- c(res, extra_qc[k])
   #       k <- k + 1
   #     }
   #     else
   #     {
   #       res <- c(res, color_lst[j])
   #       j <- j + 1
   #     }
   #   }
   #   color_lst <- res
   # }
   color_lst<-colors

    a + geom_bar(position="dodge", stat="identity", width=0.9) +
      geom_bar(aes(linetype=under_50_percent,color = under_50_percent, size=under_50_percent),position="dodge", stat="identity", width=0.9) +
      guides(linetype='none',fill=guide_legend(ncol=1))+
       scale_size_manual(values=c(0.3,0.8), guide = 'none') +
       scale_colour_manual(values = c("black", "gray29"), guide = 'none') +
      {if (type != "tf") scale_y_continuous(labels = scales::scientific)} +
      facet_wrap( ~ Name, scales=scales) +
      theme_bw() +scale_linetype_manual(values=c("solid","58"))+
      labs(x=x, y=y, title=Title, fill=element_blank()) +
      theme(
        plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
        axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
        axis.text=element_text(size=11, face="bold"),
        axis.text.x=axis.text.x,
        legend.title=element_text(face="bold", size=12),
        legend.text=element_text(face="bold",size=12),                  #sets legend text
        strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
        panel.grid.major=element_blank()) +
      {if (type != "tf") geom_errorbar(
        aes(ymin=Av-Std, ymax=Av+Std),position=position_dodge(0.9), width=.2)} +
      {if (type == "tf") geom_errorbar(
        aes(ymin=RelAmounts_Ave-RelAmounts_Std, ymax=RelAmounts_Ave+RelAmounts_Std),
        position=position_dodge(0.9), width=.2)} +
      scale_fill_manual(values = color_lst)
    #
  }

