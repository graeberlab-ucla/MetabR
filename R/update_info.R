#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param The path to the original information excel sheet of the analysis.
#'
#' @return An updated excel sheet with sample run orders and the quality control blanks, 250ks, and samples.
#'
#' @examples update_info(info_dir = "/Users/juyeki/lab_folder")
#'
#' @export
#'

update_info <- function(info_dir, folder = NULL )
{
  library(readxl)
  library(xlsx)
  library(dplyr)

  ICS <- ifelse (grepl("ICS", info_dir), ICS <- TRUE, ICS <- FALSE)
  Medium <- ifelse (grepl("Medium", info_dir), Medium <- TRUE, Medium <- FALSE)
  Cond_Medium <- ifelse (grepl("Cond Medium", info_dir), Cond_Medium <- TRUE, Cond_Medium <- FALSE)

  if (is.null(folder))
  {
    if (ICS)
    {
      order_dir <- paste0(info_dir, "/ICS/data")
    } else if (Cond_Medium){
      order_dir <- paste0(info_dir, "/AA/data")
    } else if (Medium) {
      order_dir<-paste0(info_dir,"/Medium/data")
    }else {
      order_dir <- paste0(info_dir, "/AA/data")
    }
  } else{
    order_dir <- paste0(info_dir, folder)
  }

  #info sheet
  setwd(info_dir)
  file_name <- list.files()[grep('.xls[x]?',list.files())]
  info <- read_excel(file_name)
  info <- info[!is.na(info$Sample),] #removes non-sample rows (sometimes there is a skipped line)
  info <- info[!grepl("QC", info$Sample),] #removes extra QCs (for now)

  #Run order
  setwd(order_dir)
  file_data <- file.info(list.files(pattern = "*.raw"))
  file_data <- file_data[order(as.POSIXct(file_data$mtime)), ]

  run_order <- row.names(file_data)
  run_order <- gsub(x = run_order, pattern = ".raw", replacement = "")
  run_order <- as.data.frame(run_order)
  colnames(run_order)[1] <- "samples"
  run_order$Run.Order <- rownames(run_order)

  #make info$sample namings match run_order$Sample
  samples <- as.vector(run_order$samples)
  qc <- samples[grepl("QC", samples)]
  qc <- qc[order(qc)]
  samples <- samples[!grepl("QC", samples)]
  samples <- samples[order(samples)]
  samples <- c(samples, qc)

  #adding empty rows to info for QC's
  to_add <- length(samples) - nrow(info)
  for (i in 1:to_add)
    info[nrow(info)+1,] <- NA
  info <- cbind(samples, info)
  info$Sample <- NULL

  #adding run order to info
  info <- inner_join(x = info, y = run_order, by = "samples")
  info <- info[,c(1,  ncol(info), 2:(ncol(info)-1))] #reordering
  colnames(info)[1] <- "Sample"

  #filling condition column if missing (the QCs)
  for (i in 1:nrow(info))
    if (is.na(info$Condition[i]))
      info$Condition[i] <- as.character(info$Sample[i])

  info$Run.Order<-(as.numeric(info$Run.Order))  #added this line for proper run order sorting
  #writing out new info sheet
  setwd(info_dir)
  write.xlsx(x = info, file = file_name, row.names = F, showNA = F)
}
