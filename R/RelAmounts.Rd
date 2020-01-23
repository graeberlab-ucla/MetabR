\name{RelAmounts}
\alias{RelAmounts}
\title{This function takes in a data frame and optional parameters and transforms it into a data frame with the relative metabolite amounts.}
\usage{
RelAmounts(DF = data2, anova = TRUE, type = '250K')
}
\arguments{
\item{DF}{This parameter is a dataframe (data2) that contains information about the sample. The dataframe should contain the following columns: Name, Iso, Condition, Exp, Value, Used_ID, KEGG.ID, Nr.C, Rt, and Formula columns.}
\item{anova}{Optional parameter (default is FAlSE). Set to true if ANOVA is to be calculated (for actual analysis).}
\item{type}{Optional parameter. Set to '250K' if running 250K QC script.}
}
\description{
This function transforms the data into a form that can be used for plotting and visualization of the relative amounts of the different metabolites in the analysis.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}
