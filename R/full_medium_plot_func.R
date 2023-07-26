#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#' @param metabolites list of metabolites used to run experiment
#' @param df amounts3 or amounts
#' @param n n = num_conditions
#' @param type type = "tf" or other
#' @param index index from maven_analysis pre-determined groupings
#' @param title_type "nonpathway1", "nonpathway2", or "nonpathway3", or something else
#' @param only_M0 Not currently being used
#' @param mediumtitle mediumtitle="Plots-Relative-To-Fresh"
#'
#' @return a plot of the relative amounts or samples relative to unspent of specific metabolites
#' @export
#'
# @examples full_medium_plot_func(unname(unlist(lst)),a, n = num_conditions, type = "tf", title_type = "medium", mediumtitle="Plots-Relative-To-Fresh", color_lst=amounts3_colors)

full_medium_plot_func<-function (metabolites, df, n, type=NULL, index = NULL, title_type,
                             only_M0 = NA,  mediumtitle = NULL)
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
  } else {
    ending = ""
  }
  fresh <- grepl("Relative-To-Fresh", mediumtitle)
  if (sum(grepl("Exp", names(df))) >= 1) {
    if (type == "tf") {
      met <- df %>%
        filter(Name %in% metabolites) %>%
        mutate(Name = paste(Name, Sig, sep = " "),
               Sig = "",
               err_ymin = ifelse(
                 (RelAmounts_Ave-RelAmounts_Std < 0) & (!fresh), 0, RelAmounts_Ave-RelAmounts_Std),
               err_ymax = RelAmounts_Ave + RelAmounts_Std)
      stopifnot(length(unique(met$Name)) >= 1)

      Title = paste0("Amounts of ", ending, " (not corrected for blank values)")
      if(!is.null(mediumtitle)) Title<-mediumtitle
      x = ""
      y = "Amounts"
      a <- ggplot(met, aes(Condition, RelAmounts_Ave,
                           group = Condition, fill = Condition, label = Sig)) +
        {if(fresh) geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.25)}

    }
    else {
      met = subset(df, Name %in% metabolites) %>%
        mutate(Name = paste(Name, Sig, sep = " "),
               Sig = "",
               err_ymin = ifelse(Av-Std < 0, 0, Av-Std),
               err_ymax = Av + Std)
      stopifnot(length(unique(met$Name)) >= 1)
      Title = paste0("Amounts of ", ending)
      if(!is.null(mediumtitle)) Title<-mediumtitle
      x = ""
      y = "Amounts"
      a <- ggplot(met, aes(Condition, Av,
                           group = Condition, fill = Condition, label = Name))
    }
  }
    a + geom_bar(position="dodge", stat="identity", width=0.9) +
      geom_bar(aes(linetype=under_50_percent,color = under_50_percent, size = under_50_percent),position="dodge", stat="identity", width=0.9) +
      guides(linetype='none',fill=guide_legend(ncol=1))+
       scale_size_manual(values=c(0.3,0.8), guide = 'none') +
       scale_colour_manual(values = c("black", "gray29"), guide = 'none') +
       scale_linetype_manual(values=c("solid","58")) +
       scale_fill_manual(values = colors) +
      {if (type != "tf") scale_y_continuous(labels = scales::scientific)} +
      facet_wrap( ~ Name, scales = "free") +
      theme_bw() +
      labs(x=x, y=y, title=Title, fill=element_blank()) +
      theme(
        plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
        axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
        axis.text=element_text(size=11, face="bold"),
        axis.text.x=element_blank(),
        legend.title=element_text(face="bold", size=12),
        legend.text=element_text(face="bold",size=12),                  #sets legend text
        strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
        panel.grid.major=element_blank()) +
      geom_errorbar(
        aes(ymin = err_ymin, ymax = err_ymax),
        position=position_dodge(0.9), width=.2)
  }

