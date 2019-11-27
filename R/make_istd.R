#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param A dataframe containing the tracefinder data.
#'
#' @return A list of two dataframes, data (the updated data) and Std (the dataframe to be used to make the ISTD plot)
#'
#' @examples make_istd(data)
#'
#' @export
#'

make_istd <- function(data, remove)
{
  if(remove)
  {
    #background subtraction
    if (sum(grepl('blank|water|buffer', names(data))) > 0) {
      blank <- data[,grep('blank|water|buffer', names(data))]
      data <- data[,-grep('blank|water|buffer', names(data))]
      blank[is.na(blank)] <- 0
      print("background subtraction")
    }

    if (sum(grepl('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]|200[Kk]', names(data))) > 0) {                          ##collect all standard samples and remove from regular samples
      Standard <- data[,c(1:3, grep('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]|200[Kk]', names(data)))]
      data <- data[,-grep('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]|200[Kk]', names(data))]
      print("removing standards from sample data")
    }
  }

  if(remove)
  {
    Std <- filter(data, grepl('Std', data$Iso)) %>%
      dplyr::select(-Iso) %>%
      gather(Sample, Value, -Used_ID) %>%
      filter(!(grepl('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]', Sample))) %>%
      group_by(Used_ID) %>%
      mutate(Av = mean(Value, na.rm=T),
             Std = sd(Value, na.rm=T),
             CV = Std / Av,
             ID = paste0(gsub(' Std','',Used_ID), ' ',round(CV*100,0),'%')) %>%
      ungroup()
  } else{
    Std <- filter(data, grepl('Std', data$Iso)) %>%
      dplyr::select(-Iso) %>%
      gather(Sample, Value, -Used_ID) %>%
      group_by(Used_ID) %>%
      mutate(Av = mean(Value, na.rm=T),
             Std = sd(Value, na.rm=T),
             CV = Std / Av,
             ID = paste0(gsub(' Std','',Used_ID), ' ',round(CV*100,0),'%')) %>%
      ungroup()
  }

  result <- list()
  result$data <- data
  result$Std <- Std

  return(result)
}
