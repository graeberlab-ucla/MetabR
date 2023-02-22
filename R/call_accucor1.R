#' call_accucor1
#' Calls Accucor Natural Isotope Correction for single labeled sets (e.g C13, N15, D).
#' Converts mid_output data to format for Accucor1.
#' install.packages("accucor")
#' @param mid_output "isotopologue data uncorrected.csv".
#' @param abbrev The downloaded google sheets Abbrev_New2 sheet. Columns Abb, Neutral.Formula, and KEGG.ID are used in function.
#' @param Title Title to use for naming the corrected file
#' @param info Sample info sheet, commonly used when running Metabolomics data. Only the info$Condition column is used (make sure this contains factors).
#'
#' @return corr2 NIC corrected version. Not yet in the format for "MIDS corrected output"
#'
#' @importFrom dplyr select arrange mutate rename
#' @importFrom tidyr unite gather spread separate
#' @export
#'
#' @examples temp<-call_accucor1(mid_output, Abbrev, Title, info)
#'
call_accucor1<-function(mid_output, abbrev, Title, info){
  # require(accucor)
  # require(dplyr)
  # require(stringr)

  resolution= 70000
  if(dim(subset(mid_output, startsWith(as.character(Iso), "C13-")))[1]!=0) label="C13"
  if(dim(subset(mid_output, startsWith(as.character(Iso), "N15-")))[1]!=0){
    label="N15"
    resolution=140000
  }
  if(dim(subset(mid_output, startsWith(as.character(Iso), "D-")))[1]!=0) label = "D"

  corr_success<-TRUE

  info$Condition <- factor(info$Condition, levels=as.character(unique(info$Condition)))
  mid_output$Condition<- factor(mid_output$Condition, levels = levels(info$Condition))

  uncorr<-mid_output%>%
    select(Name, Condition, Iso, which(grepl("Exp", names(mid_output)))) %>%
    gather(Exp, Value, -Name, -Condition, -Iso) %>%
    arrange(Name, Condition) %>%
    unite(Condition_Exp, c(Condition, Exp), sep='_') %>%
    mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
    spread(Condition_Exp, Value)

  uncorr$Iso<-gsub("-", "-label-", uncorr$Iso)

  m_ind<-match(uncorr$Name, abbrev$Abb)
  uncorr$Formula<-abbrev$Neutral.Formula[m_ind]
  uncorr<- uncorr %>% select(Name, Formula, Iso, everything())
  names(uncorr)[names(uncorr) == "Name"] <- "Compound"
  names(uncorr)[names(uncorr) == "Iso"] <- "IsotopeLabel"

  formula_na<-uncorr[which(is.na(uncorr$Formula)),] #save rows with no Neutral Formula to incorporate back in later

  if(dim(formula_na)[1]!=0) {
    print(paste(list(unique(formula_na$Compound)), "do not have Neutral Formulas on the Abbrev sheet and will NOT be corrected"))
    uncorr2<-uncorr[-which(is.na(uncorr$Formula)),] #remove rows that have no Neutral Formula as it causes problems with natural_abundance_correction()
    }else{
    uncorr2<-uncorr
  }


  print("Performing Natural Isotope Correction...")
  corrected <- accucor::natural_abundance_correction(
    path = uncorr2,
    resolution = resolution, output_base = TRUE)

  file.rename("TRUE_corrected.xlsx", paste0(Title,"_Accucor1_function_output.xlsx"))

  corr<-corrected$Corrected

  if(dim(formula_na)[1]!=0){
    #may need to state this as an unsuccessful correction
  isochart<-as.data.frame(stringr::str_split_fixed(formula_na[,3], "\\-", 3))
  isochart <- stringr::replace(isochart, isochart=='', '0')
  formula_na$IsotopeLabel<-isochart[,3]
  formula_na$Formula<-NULL
  names(formula_na)[2]<-names(corr)[2]
  if(all(names(formula_na)==names(corr))){
    corr<-rbind(corr, formula_na)
  }else{
    print("Compounds that didn't contain a Neutral Formula in the Abbrev sheet have been removed from data")
    corr_success==FALSE
  }
  }
  #adjust uncorr to be in the same format as corr
  uncorr$IsotopeLabel<-gsub("C12 PARENT", "0", uncorr$IsotopeLabel)
  uncorr$IsotopeLabel<-gsub(paste0(label,"-label-"), "", uncorr$IsotopeLabel)
  uncorr$IsotopeLabel<-paste0(stringr::str_pad(uncorr$IsotopeLabel, 2, pad = "0"))
  uncorr$Comp_Iso<-paste(uncorr$Compound, uncorr$IsotopeLabel, sep="_")

  names(corr)[grepl("Label",names(corr))]<-"Iso"
  corr$Iso<-paste0(stringr::str_pad(corr$Iso, 2, pad = "0"))#makes it so all numbers have two digits. 00, 01, 02, ... keeps order correct.
  corr<-corr%>% mutate(Comp_Iso=paste(Compound, Iso, sep="_"))

  #Get rid of extra rows that were added during correction (Accucor generates all possible isotopologues).
  #Only keeping isotopes from original uncorrected, filtered data.
  corr<-corr[corr$Comp_Iso %in% uncorr$Comp_Iso,]

  corr2<-corr%>%select(-Comp_Iso)%>%
    gather(Condition_Exp, Value, -Compound, -Iso) %>%
    separate(Condition_Exp, c('Condition', 'Exp'), sep='_')%>%
    rename(Name=Compound)

  corr2$Iso<-ifelse(corr2$Iso=='00', "C12 PARENT", paste0(label,"-", corr2$Iso))
  m_ind<-match(corr2$Name, abbrev$Abb)
  corr2$KEGG.ID<-abbrev$KEGG.ID[m_ind]

  corr2["Norm_Av"] <- corr2["Norm_Std"] <- corr2["CV"] <- corr2["Av"] <- corr2["Nr.C.x"] <-
    corr2["MID1"] <- corr2["MID2"] <- corr2["MID3"] <- corr2["Nr.C.y"] <- corr2["Old Uncorrected Value"] <-
    corr2["ANOVA"] <- corr2["Sig"] <- NA
  if(corr_success==TRUE){
    print("Correction completed.")
  }else{
    print("Correction was UNsuccessful. Review Compounds with NA in the Neutral Formula Column.")
  }
  return(corr2)
}

