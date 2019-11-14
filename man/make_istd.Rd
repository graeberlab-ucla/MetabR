\name{make_istd}
\alias{make_istd}
\title{This function takes in a data frame containing tracefinder data and outputs a list of two dataframes, data (a modified dataframe of the parameter data) and std (the dataframe to be used to make the ISTD plot).}
\usage{
make_istd(data = data)
}
\arguments{
\item{data}{This parameter is a dataframe containing the tracefinder data.}
}
\description{
This function performs background subtraction on the tracefinder data and creates a dataframe in the correct format create the ISTD plot.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}

