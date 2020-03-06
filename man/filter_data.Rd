\name{filter_data}
\alias{filter_data}
\title{Filters raw metabolomic data.}
\usage{
filter_data(data, info, type = "QC")
}
\arguments{
\item{data}{a dataframe containing raw metabolomic data. Each row is one isotopologue with the metabolite and isotopologue in the first column and the rest of the columns being samples.}
\item{info}{The sample info sheet with the conditions in a column titled 'Condition'.}
\item{type}{An optional parameter. The default outputs a dataframe with the corrected values. If set to 'QC', a vector of isotopologues to be filtered is outputted instead.}
}
\description{
This function filters the raw metabolomic data. For each isotopologue, values are changed to zero if less than half of a specific condition's replicates have values and no other condition for that isotopologue has more than half of its samples have values.
}
\author{
Jordan Uyeki
}
