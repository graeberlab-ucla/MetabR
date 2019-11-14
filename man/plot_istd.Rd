\name{plot_istd}
\alias{plot_istd}
\title{This function creates a plot which can be used to visualize the analysis of the internal standards}
\usage{
plot_istd(Std = std, info = info_ordered, title = title)
}
\arguments{
\item{Std}{A dataframe created by make_istd().}
\item{info}{The information sheet that is ordered in the order that the samples were ran.}
\item{title}{The title of the resulting pdf.}

}
\description{
This function will create the internal standards (ISTD) plot and output it as a pdf file into the current working directory.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}
