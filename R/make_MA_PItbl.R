#' Make MetaboAnalyst compatible peak intensity table CSV
#'
#' Uses "MavenData" file and reformats to be used with MetaboAnalyst (https://www.metaboanalyst.ca/)
#'
#' @param data_dir Path to network data folder; from pipeline_meta
#' @param info Sample info sheet dataframe; from pipeline_meta
#' @param Abbrev Abbreviation dataframe; from pipeline_meta
#' @param abbr Logical; 'T' to replace names with abbreviations
#'
#' @return Peak intensity CSV in data_dir / network folder
#' @export
make_MA_PItbl <- function(data_dir, info, Abbrev, abbr = F) {
  md <- readxl::read_excel(
    list.files(data_dir, pattern = "^MavenData", full.names = T))
  md <- md[complete.cases(md),] # Remove "A" rows
  md[,2:ncol(md)] <- sapply(md[,2:ncol(md)], as.numeric) # Numeric type
  names(md) <- sub("_neg.mzXML|_pos.mzXML", "", names(md)) # Fix colnames
  names(md)[1] <- "Sample"
  md <- md[!grepl(" Std", md$Sample),] # Remove stds
  # Sum metabolite parent with isotopologue if labeled
  if(labeled == "full" | labeled == "partial" | grepl("PARENT",md[1,1])){
    md$Sample <- sub(" C\\d\\d[ |-][P|l].*", "", md$Sample) # Carbon labeled
    md$Sample <- sub(" N\\d\\d[ |-][P|l].*", "", md$Sample) # Nitrogen labeled
    md <- md %>%
      group_by(Sample) %>%
      summarise(across(where(is.numeric), sum))
  }
  if (abbr) { # Replace name with abbrev
    # First pass
    md_abbr <- right_join(Abbrev[,c("Used_ID", "Abb")], md, by = c("Used_ID" = "Sample"))
    md_abbr[is.na(md_abbr$Abb),]$Abb <- Abbrev[Abbrev$Synonyms %in% md_abbr[is.na(md_abbr$Abb),]$Used_ID,]$Abb # Pass through synonyms
    md_abbr$Used_ID <- NULL
    names(md_abbr)[1] <- "Sample"
    md <- md_abbr
  }
  if (!all(is.na(info$Cell.Number))){ # Add Cell.Number row if provided
    info$Cell.Number <- ifelse(is.na(info$Cell.Number), "1", info$Cell.Number)
    md <- rbind(c("Cell.Number", info$Cell.Number), md)
  }
  md <- rbind(c("Condition", info$Condition), md) # Add labels as first row
  md <- md %>% select(-tidyr::starts_with('QC-250'))   # Filter out QC-250ks
  md[1,] <- ifelse(grepl("QC", md[1,]), "QC", md[1,]) # Change blank condition to "QC"
  # Replace 0s with NA disregarding QC-blanks
  md_NA <- md[2:nrow(md),2:(which(grepl("QC", md))[1] - 1)]
  md_NA <- as.data.frame(sapply(md_NA, as.numeric))
  md_NA[md_NA == 0] <- NA
  md_NA <- as.data.frame(sapply(md_NA, as.character))
  md[2:nrow(md),2:(which(grepl("QC", md))[1] - 1)] <- md_NA

  write.csv(md, file = paste0(data_dir, "/", Title, "_MetaboAnalystPeakIntensityTbl.csv"),
            row.names = F)
}
