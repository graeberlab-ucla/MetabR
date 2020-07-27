#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param The name of the excel sheet that contains the Maven processed data, the sheets of the excel sheet to include, and the name of the output file
#'
#' @return Outputs an excel sheet, formatted in the correct format for processing.
#'
#' @examples format_maven(data_sheet = "data.xlsx", sheets = c(1:2), output = "output.csv)
#'
#' @export
#'


format_maven <- function(data_sheet, sheets, output)
{
  require(readxl)
  combined <- data.frame()

  for(sheet in sheets)
  {
    data <- data.frame(read_excel(path = data_sheet, sheet = sheet, col_names = F))
    new_colnames <- c("Name", "Iso", as.character(data[1, 2:ncol(data)]))

    data <- rbind(c(rep(NA, ncol(data))), data, c(rep(NA, ncol(data))))

    result <- data.frame(matrix(ncol = ncol(data)+1))
    colnames(result) <- new_colnames

    for(i in 1:nrow(data))
    {
      if(grepl(pattern = "_Std", x = data[i,1]))
      {
        data[i,1] <- gsub(pattern = "_Std", replacement = "", x = data[i,1])
        data[i+1, 1] <- "Std"
      }
    }

    index <- 1
    for(i in 1:sum(is.na(data[,1]) ))
    {
      if(i != sum(is.na(data[,1])))
      {
        current_NA <- which(is.na(data[,1]))[i]
        next_NA <- which(is.na(data[,1]))[i + 1]
        num_isotopes <- next_NA - current_NA - 2

        metabolite <- data[current_NA+1, 1]
        result[index:(index + num_isotopes - 1), 1] <- metabolite

        result[index:(index + num_isotopes - 1), 2] <- data[(current_NA + 2):(current_NA + 1 + num_isotopes), 1]
        result[index:(index + num_isotopes - 1), 3:(ncol(data)+1)] <- data[(current_NA + 2):(current_NA + 1 + num_isotopes), 2:ncol(data)]
        index = index + num_isotopes
      }
    }
    combined <- rbind(combined, result)
  }

  write.csv(combined, output, row.names = F)
}
