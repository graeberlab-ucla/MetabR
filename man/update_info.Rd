\name{update_info}
\alias{update_info}
\title{Updates sample information sheet.
}
\usage{
update_info(info_dir = "/Users/juyeki/lab_folder")
}
\arguments{
\item{info_dir}{The path to the folder that contains the information sheet that is to be updated.}
}
\description{
This function searches for the information sheet in info_dir (only excel sheet in info_dir should be information sheet), opens it, adds the necessary information, and outputs the new sheet, overriding the original.
}
\author{
Gerardo Javier, Daniel Tan, Jordan Uyeki
}
