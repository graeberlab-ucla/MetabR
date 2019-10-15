\name{my_make_RelAmounts_QC}
\alias{my_make_RelAmounts_QC}
\title{This function takes in a data frame and transforms it into a data frame with the relative metabolite amounts.}
\usage{
my_make_RelAmounts_QC(DF = data2)
}
\arguments{
\item{DF}{This parameter is a dataframe (data2) that contains information about the sample.  The dataframe should contain the following columns: Name, Iso, Condition, Exp, Value, Used_ID, KEGG.ID, Nr.C, Rt, and Formula columns.}
}
\description{
This function transforms the data into a form that can be used for plotting and visualization of the relative amounts of the different metabolites in the analysis.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}
