\name{make_heatmap}
\alias{make_heatmap}
\title{This function takes in a dataframe with the relative amounts, a data frame with information on the samples, and the heatmap color scheme to be used. It outputs a pdf file containing an annotated heatmap.}
\usage{
make_heatmap(matrix = RelA, samples = samples, heat.color = normal)
}
\arguments{
\item{matrix}{A dataframe containing the relative amounts.}
\item{samples}{A dataframe with information on the samples.}
\item{heat.color}{The color scheme that is to be used in the heatmap.}
}
\description{
This function creates an annotated heatmap, allowing visalization of sample Relative amount metabolite values, Norvaline values, Cell Number readings, and the experimental condition that each sample was run under.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}


