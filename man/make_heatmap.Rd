\name{make_heatmap}
\alias{make_heatmap}
\title{This function takes in a dataframe with the relative amounts, a data frame with information on the samples, and the heatmap color scheme to be used. It outputs a pdf file containing an annotated heatmap.}
\usage{
make_heatmap(matrix = RelA, samples = samples, heat.color = normal, width = 8, height = 5)
}
\arguments{
\item{matrix}{A dataframe containing the relative amounts.}
\item{samples}{A dataframe with information on the samples.}
\item{heat.color}{The color scheme that is to be used in the heatmap.}
\item{width}{Optional parameter to set width dimensions of pdf output plot.}
\item{height}{Optional parameter to set height dimensions of pdf output plot.}
}
\description{
This function creates an annotated heatmap, allowing visalization of sample Relative amount metabolite values, Norvaline values, Cell Number readings, and the experimental condition that each sample was run under. The 'width' and 'height' parameters allow user to manually set pdf output dimensions. Otherwise, pheatmap::pheatmap calculates the dimensions, but sometimes annotation legends/titles are cut off.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}


