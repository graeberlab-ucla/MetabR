\name{correct_iso}
\alias{correct_iso}
\title{Natural Abundance correction using IsoCorrectoR method for C,N,CN
}
\usage{
correct_iso(df,label,correction)
}
\arguments{
\item{df}{The isotopomer dataframe with values for correction}
\item{label}{The labelling type of the data e.g. "C","N","CN",etc}
\item{correction}{The desired molecule(s) to correct for}
}
\description{
This function takes in the isotopomer data table with uncorrected values (from Maven or Tracefinder), and conducts natural abundance correction for the desired molecule(s) using the IsocorrectoR method. This function can correct for Carbon, Nitrogen, and dual Carbon-Nitrogen labelling. 
}
\author{
Daniel J Tan
}