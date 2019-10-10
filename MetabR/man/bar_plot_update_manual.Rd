\name{bar_plot_update_manual}
\alias{bar_plot_update_manual}
\title{Creates individual barplots for each metabolite defined in the pathway.}
\usage{
bar_plot_update_manual(a, met, Title, axis.text.x, scales, type = NULL, num_cond = NULL)
}
\arguments{
\item{a}{a ggplot barplot that has the conditions on the x axis and the normalized averages along the y axis.}
\item{met}{The subsetted data frame containing only the metabolites in the pathway of interest.}
\item{Title}{The title that will go on the top of each page of the output pdf.}
\item{axis.text.x}{Customization of x-axis. Default of bar_plot_update is element_blank().}
\item{scales}{Customization of scales. Default of bar_plot_update is "free".}
\item{type}{The type of quality control or analysis to be run. eg. "250k". This argument is inherited from bar_plot_update when it is called.}
\item{num_cond}{The number of conditions in the experiment}
}
\description{
This function is internally called by bar_update_manual and creates the relative amount barplots for each of the metabolites defined in the pathway of interest.
}

