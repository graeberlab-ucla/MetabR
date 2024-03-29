#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param A dataframe containing raw metabolomic data. A dataframe containing sample information, specifically condition. An optional parameter 'type' set to 'QC' for regular QC.
#'
#' @return A dataframe of the filtered data. Or if type = 'QC', a vector of the isotopologues to be filtered.
#'
# @examples filter_data(data = data, info = info, type = 'QC')
#'
#' @export
#'


filter_data <- function(data, info, type = "reg")
{
  cond <- as.vector(info$Condition)
  cond <- cond[-which(grepl("blank|QC.|QC-", cond))]
  replicates <- length(cond)/length(unique(cond))

  first_sample <- 2
  last_sample <- which(grepl("blank|QC.|QC-", colnames(data)))[1] - 1

  data$to_remove <- NA
  groupings <- seq(first_sample, last_sample, replicates)
  thresh <- floor(0.5*replicates)

  for(i in 1:nrow(data))
  {
    result <- TRUE
    if(all(is.na(data[i,first_sample:last_sample])))  #### deals with all NA rows  (just remove???)
      result <- FALSE
    for(j in groupings)
    {
      count <- 0
      for(k in j:(j+replicates-1))
      {
        if(!is.na(data[i,k]))
          count <- count + 1
        if(count > thresh)
          result <- FALSE
      }
    }
    data$to_remove[i] <- result
  }

  print(as.vector(data$Compound[which(data$to_remove == TRUE)]))   #show ones where values are removed
  if(type == "QC")
    return(as.vector(data$Compound[which(data$to_remove == TRUE)]))

  for(i in 1:nrow(data))
  {
    if(data$to_remove[i] == TRUE)
      data[i,first_sample:last_sample] <- NA
  }
  data$to_remove <- NULL

  return(data)
}

