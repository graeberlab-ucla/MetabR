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
    md$Sample <- sub(" C\\d\\d[ |-][P|l].*", "", md$Sample)
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
  md <- rbind(c("Condition", info$Condition), md) # Add labels as first row
  md <- md %>% select(-starts_with('QC'))   # Filter out QC
  write.csv(md, file = paste0(data_dir, "/", Title, "_MetaboAnalystPeakIntensityTbl.csv"),
            row.names = F)
}
