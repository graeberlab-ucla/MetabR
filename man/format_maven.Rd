\name{format_maven}
\alias{format_maven}
\title{Formats the Maven-processed data into the correct format for further analysis.}
\usage{
format_maven(data_sheet = "data.xlsx", sheets = c(1:2), output = "output.csv")
}
\arguments{
\item{data_sheet}{An excel sheet that contains the Maven processed data.}
\item{sheets}{The sheets in the excel sheet that are to be included in the output file.}
\item{output}{The name of the output file.}
}
\description{
This function is formats Maven data into the correct format for Maven processing.
}
\author{
Daniel Tan, Jordan Uyeki
}
