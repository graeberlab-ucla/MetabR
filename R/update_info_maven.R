
#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param The path to the original information excel sheet of the analysis.
#'
#' @return An updated excel sheet with sample run orders and the quality control blanks, 250ks, pools, and samples. Used for Maven only.
#'
# @examples update_info(info_dir = "/Volumes/MavenData ALMOHI/XXX lab/file_folder")
#'
#' @export
#'
update_info_maven <- function(info_dir, folder = NULL ){

library(readxl)
library(xlsx)
library(dplyr)
#info sheet

setwd(info_dir)
file_name <- list.files(pattern='.xls[x]?')[!grepl("\\$|~|raw data|orig|Vanq|Accucor|Metabo",list.files(pattern='.xls[x]?'))][1]
info <- read_excel(file_name)
info <- info[!is.na(info$Sample),] #removes non-sample rows (sometimes there is a skipped line)
info <- info[!grepl("QC", info$Sample),] #removes extra QCs (for now)
#setwd("../")
#Run order
order_dir<-paste0(info_dir,"/RAW files")
setwd(order_dir)
file_data <- file.info(list.files(pattern = "*.raw"))
file_data <- file_data[order(as.POSIXct(file_data$mtime)), ]

run_order <- row.names(file_data)
run_order <- gsub(x = run_order, pattern = ".raw", replacement = "")
run_order <- as.data.frame(run_order)
colnames(run_order)[1] <- "samples"
run_order$Run.Order <- rownames(run_order)

#make info$sample namings match run_order$Sample
#Alexzandra added 1/25/22 - lists blanks before 250ks
samples <- as.vector(run_order$samples)
qc_blank <- samples[grepl("blank", samples)]
qc_blank <- qc_blank[order(qc_blank)]

qc_250k <- samples[grepl("250K", samples) | grepl("250k", samples) | grepl("50k", samples) | grepl("50K", samples)]
qc_250k <- qc_250k[order(qc_250k)]

qc_pool <- samples[grepl("pool|Pool", samples)]
qc_pool <- qc_pool[order(qc_pool)]

samples <- samples[!(samples %in% qc_blank) & !(samples %in% qc_250k) &!(samples %in% qc_pool)]
samples <- samples[order(samples)]
samples <- c(samples, qc_blank)
samples <- c(samples, qc_pool)
samples <- c(samples, qc_250k)

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
setwd("../")
setwd(info_dir)
write.xlsx(x = info, file = file_name, row.names = F, showNA = F)
}


