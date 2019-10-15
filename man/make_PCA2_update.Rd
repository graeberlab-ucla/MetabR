\name{make_PCA2_update}
\alias{make_PCA2_update}
\title{Creates a pdf of PCA plots and correlation Circle plot.
}
\usage{
make_PCA2_update(matrix = RelA, a = 2, b = 3, cutoff = 0.4)
}
\arguments{
\item{matrix}{The matrix containing the data to be used.}
\item{a}{The PC to be used for the x-axis (the default is 1).}
\item{b}{The PC to be used for the y-axis (the default is 2).}
\item{cutoff}{The cutoff value for correlation (metabolites with correlations greater than the 'cutoff' will be excluded from the output plots.)}
}
\description{
This function takes in a data matrix and exports a pdf file containing a scree plot, a pair plot with the first five PCs, and a PCA plot with PCs a and b.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}
