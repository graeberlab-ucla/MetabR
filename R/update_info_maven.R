#' Update Maven info sheet
#'
#' Updates excel sheet with sample run orders and the quality control blanks, 250ks, pools, and samples. Used for Maven only.
#'
#' @import dplyr
#' @param info_dir string; path to directory containing original info excel sheet for dataset
#'
#' @return None; writes updated Excel sheet info sheet in specified directory
#' @export
#'
update_info_maven <- function(info_dir){
#info sheet
setwd(info_dir)
file_name <- list.files(pattern='.xls[x]?')[!grepl("\\$|~|raw data|orig|Vanq|Accucor|Metabo",list.files(pattern='.xls[x]?'))][1]
info <- openxlsx::read.xlsx(file_name)
info <- info[!is.na(info$Sample),] #removes non-sample rows (sometimes there is a skipped line)
info <- info[!grepl("QC", info$Sample, ignore.case = T),] #removes extra QCs (for now)
#setwd("../")
#Run order
order_dir<-paste0(info_dir,"/RAW files")
if(dir.exists(order_dir)){
  setwd(order_dir)
}else{
  stop(paste0("RAW files folder does not exists"))
}
raw_files <- list.files(pattern = "*.raw")

if (length(raw_files) == 0) stop(paste0("RAW files folder doesn't contain any files"))

file_data <- file.info(raw_files)
file_data <- file_data[order(as.POSIXct(file_data$mtime)), ]

run_order <- row.names(file_data)
run_order <- gsub(x = run_order, pattern = ".raw", replacement = "")
run_order <- as.data.frame(run_order)
colnames(run_order)[1] <- "samples"
run_order$Run.Order <- rownames(run_order)

#make info$sample namings match run_order$Sample
#Alexzandra added 1/25/22 - lists blanks before 250ks
samples <- as.vector(run_order$samples)

qc_blank <- samples[grepl("blank", samples, ignore.case = T) & !grepl("PB", samples, ignore.case = T)]
qc_blank <- qc_blank[order(qc_blank)]

qc_250k <- samples[grepl("250K", samples) | grepl("250k", samples) | grepl("50k", samples) | grepl("50K", samples) | grepl("25k", samples) | grepl("25K", samples)]
qc_250k <- qc_250k[order(qc_250k)]

qc_pool <- samples[grepl("pool|Pool", samples, ignore.case = T)]
qc_pool <- qc_pool[order(qc_pool)]

#processing blanks
qc_other<-samples[grepl("QC&|PB", samples, ignore.case = T)]
qc_other <- qc_other[order(qc_other)]

remaining_samples <- samples[!(samples %in% qc_blank) & !(samples %in% qc_250k) &!(samples %in% qc_pool) &!(samples %in% qc_other)]

#other samples such as "ISTD" or STD or ImP that have different concentrations - usually standards
other<-remaining_samples[which(!grepl("^[A-Z]{2}[-.][0123456789]{3}", remaining_samples))]
#other <- other[other %in% info$Sample] #may be necessary if the ordering is incorrect of the "other" samples

samples<-remaining_samples[which(grepl("^[A-Z]{2}[-.][0123456789]{3}", remaining_samples))]

samples <- samples[order(samples)]
sample.length <- length(samples)
samples <- c(samples, other)
samples <- c(samples, qc_other) #processing blank
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
for (i in (sample.length+1):nrow(info))
  #if (is.na(info$Condition[i]))
  info$Condition[i] <- as.character(info$Sample[i])

# Remove extra whitespace before and after Condition and Sample names
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
info$Condition <- trim(info$Condition)
info$Sample <- trim(info$Sample)

info$Run.Order<-(as.numeric(info$Run.Order))  #added this line for proper run order sorting
#writing out new info sheet
setwd("../")
setwd(info_dir)
openxlsx::write.xlsx(x = info, file = file_name, rowNames = F, showNA = F)
}


