#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param the isotopomer dataframe, labelling of dataset, and desired molecule(s) for correction
#'
#' @return isotopomer dataframe with a column for corrected values.
#'
#' @examples correct_iso(isotopomer_df, "C","C")
#'
#' @export
#'


library(qdap)
library(tidyverse)
library(dplyr)
library(data.table)
library(tidyr)
library(IsoCorrectoR)
library(writexl)
library(stringr)
library(MetabR)
library(readr)
library(plyr)

correct_iso<-function(df,label,correction){
  df_name<-paste(deparse(substitute(df)))
  df<-data.frame(df)
  
  if ("Exp" %in% names(df) == FALSE){
    df<-gather(df,key="Exp",value="Value",colnames(df)[grep("Exp",colnames(df))])
  }
  
  df2<-df
  colnames(df)<-gsub("\\.","-",colnames(df))
  df$Condition_Exp <- paste(df$Condition,df$Exp,sep = "_")
  df<-df[, c("Name", "Iso","Condition_Exp","Value")]
  df<-df %>% spread(Condition_Exp,Value)
  colnames(df)[1]<-"Measurements/Samples"
  df$`Measurements/Samples`<-gsub("\\+","plus",df$`Measurements/Samples`)
  df$Iso<-gsub("C12 parent","C12 PARENT",df$Iso)
  
  non_carbon_metabolites_present=FALSE
  non_carbon<-grep(paste(c("PO4","HS2O3","H3PO4","H4P2O7","PPi"),collapse="|"),df$`Measurements/Samples`)
  if (length(non_carbon)!=0){non_carbon_metabolites_present=TRUE}
  
  
  
  
  if (label=="C" & correction=="C"){ #from tracefinder 
    numeric_iso<-parse_number(df$Iso)
    df$Iso<-numeric_iso
    df$`Measurements/Samples`<-paste(df$`Measurements/Samples`,df$Iso,sep="_")
    df<-df[,-c(2)]
    rownames(df)<-1:nrow(df)
    
  } else if (label=="N" & correction=="N"){ #from maven
    df$Iso<-gsub("C12 PARENT", 0,df$Iso)
    df$Iso<-str_remove(df$Iso,"N15-")
    
    df$`Measurements/Samples`<-paste(df$`Measurements/Samples`,df$Iso,sep="_")
    df<-df[,-c(2)]
    rownames(df)<-1:nrow(df)
    
    
  } else if (label=="CN" | (label=="N" & correction=="CN")){ #from maven
    CN_carbons<-  str_split_fixed(df$Iso, "-", 3)[,2]
    CN_nitrogens<-str_split_fixed(df$Iso, "-", 3)[,3]
    
    df<-cbind(df,CN_carbons,CN_nitrogens)
    
    df$CN_carbons<-as.numeric(as.character(df$CN_carbons))
    df$CN_nitrogens<-as.numeric(as.character(df$CN_nitrogens))
    df$CN_nitrogens[is.na(df$CN_nitrogens)]<-0
    df$CN_carbons[is.na(df$CN_carbons)]<-0
    
    nitrogen_only_iso<-which(!grepl(paste(c("C13","C12"),collapse="|"),df[,2] ))
    #fixing nitrogen column
    for (i in 1:length(nitrogen_only_iso)){
      df[nitrogen_only_iso[i],grep("CN_nitrogens",colnames(df))]<-str_remove(df[nitrogen_only_iso[i],2],"N15-")
      df[nitrogen_only_iso[i],grep("CN_carbons",colnames(df))]<-0
    }
    
    df$CN_temp<-NA
    data("molecule_CN") 
    for (i in 1:nrow(df)){
      x<-which(molecule_CN$Molecule==df[i,1])
      c_only<-!grepl("N",molecule_CN[x,2])
      
      if (c_only){
        df[i,ncol(df)]<-paste("C",df[i,grep("CN_carbons",colnames(df))],sep="")
      }
      if(!c_only){
        c_iso<- paste("C",df[i,grep("CN_carbons",colnames(df))],sep="")
        n_iso<-paste("N",df[i,grep("CN_nitrogens",colnames(df))],sep="")
        cn_iso<-paste(c_iso,n_iso,sep=".")
        df[i,ncol(df)]<-cn_iso
        
      }}
    
    df$`Measurements/Samples`<-paste(df$`Measurements/Samples`,df$CN_temp,sep="_")
    df_save_for_later<-df[,1:2]
    
    df<-df[,-c(2,ncol(df),ncol(df)-1,ncol(df)-2)]}
  #USE THIS TO JOIN BACK ALL WITH MAVEN NAMES LATER
  
  if (non_carbon_metabolites_present==TRUE){
    dont_correct<-df[non_carbon,]
    df<-df[-non_carbon,]}
  
  rownames(df)<-1:nrow(df)
  
  df<-data.frame(df)
  colnames(df)[1]<-"Measurements/Samples"
  df_file_name<-paste(df_name,"converted.xlsx",sep = "_")
  write_xlsx(df,df_file_name)
  
  data("element_file")
  write_xlsx(element_file,"element_file.xlsx")
  dir.create("isocorrector_dir")
  
  if (label=="N" & correction =="N"){
    data("molecule_N")
    write_xlsx(molecule_N,"molecule_file_N.xlsx")
    corrected<-IsoCorrection(MeasurementFile=df_file_name, ElementFile="element_file.xlsx", MoleculeFile="molecule_file_N.xlsx",
                             CorrectTracerImpurity=FALSE, CorrectTracerElementCore=TRUE,
                             CalculateMeanEnrichment=TRUE, UltraHighRes=FALSE,
                             DirOut='./isocorrector_dir', FileOut='result', FileOutFormat='xls',
                             ReturnResultsObject=TRUE, CorrectAlsoMonoisotopic=FALSE,
                             CalculationThreshold=10^-8, CalculationThreshold_UHR=8,verbose=FALSE, Testmode=FALSE)
    unlink("molecule_file_N.xlsx")
    
  } else if (label=="CN" | (label =="N" & correction =="CN")){
    data("molecule_CN")
    write_xlsx(molecule_CN,"molecule_file_CN.xlsx")
    
    corrected<-IsoCorrection(MeasurementFile=df_file_name, ElementFile="element_file.xlsx", MoleculeFile="molecule_file_CN.xlsx",
                             CorrectTracerImpurity=FALSE, CorrectTracerElementCore=TRUE,
                             CalculateMeanEnrichment=TRUE, UltraHighRes=TRUE,
                             DirOut='./isocorrector_dir', FileOut='result', FileOutFormat='xls',
                             ReturnResultsObject=TRUE, CorrectAlsoMonoisotopic=FALSE,
                             CalculationThreshold=10^-8, CalculationThreshold_UHR=8,verbose=FALSE, Testmode=FALSE)
    
    unlink("molecule_file_CN.xlsx")
    
  } else if (correction=="C"){
    data("molecule_C")
    write_xlsx(molecule_C,"molecule_file_C.xlsx")
    
    corrected<-IsoCorrection(MeasurementFile=df_file_name, ElementFile="element_file.xlsx", MoleculeFile="molecule_file_C.xlsx",
                             CorrectTracerImpurity=FALSE, CorrectTracerElementCore=TRUE,
                             CalculateMeanEnrichment=TRUE, UltraHighRes=FALSE,
                             DirOut='./isocorrector_dir', FileOut='result', FileOutFormat='xls',
                             ReturnResultsObject=TRUE, CorrectAlsoMonoisotopic=FALSE,
                             CalculationThreshold=10^-8, CalculationThreshold_UHR=8,verbose=FALSE, Testmode=FALSE)
    
    unlink("molecule_file_C.xlsx")
    
  }
  
  unlink("element_file.xlsx")
  unlink(df_file_name)
  unlink("isocorrector_dir",recursive=TRUE)
  
  corrected<-corrected[["results"]][["Corrected"]]
  corrected<-setDT(as.data.frame((corrected)), keep.rownames = TRUE)[]
  colnames(corrected)[1]<-"Measurements/Samples"
  colnames(corrected)<-gsub("\\.","-",colnames(corrected))
  
  if (non_carbon_metabolites_present==TRUE){
    colnames(corrected)<-colnames(dont_correct)
    corrected<-rbind(corrected,dont_correct)
  }
  corrected<-separate(corrected,"Measurements/Samples",into=c("Name","Iso"),sep="_")
  if (label=="N" & correction=="N"){
    
    corrected$Iso<-gsub(0,"C12 PARENT",corrected$Iso)
    for (i in 1:nrow(corrected)){
      if (corrected[i,2]!="C12 PARENT"){
        corrected[i,2]<-paste("N15-",corrected[i,2],sep="")
      }
    }
    
  }else if (correction=="C"){
    
    corrected$Iso<-paste("M",corrected$Iso,sep="")
    
    
    
    
  } else if (label=="CN" | (label=="N" & correction=="CN")){
    corrected<-join_all(list(corrected,df_save_for_later), by = c('Measurements/Samples'), type = "left", match = "all")
    corrected<-separate(corrected,"Measurements/Samples",into=c("Name","Iso2"),sep="_")
    corrected <- subset(corrected, select = -c(Iso2))
    
    corrected <- corrected %>%
      select("Name", "Iso", everything())
  }
  
  colnames(corrected)<-gsub("\\.","-",colnames(corrected))
  gather_names<-colnames(corrected)[3:ncol(corrected)]
  corrected<-gather(corrected, key="Condition",value="Val",gather_names)
  
  corrected$Name<-gsub("plus","\\+",corrected$Name)
  corrected$Exp<-str_split_fixed(corrected$Condition, "_",n=2)[,2]
  corrected$Condition<-str_split_fixed(corrected$Condition, "_",n=2)[,1]
  colnames(corrected)[4]<-"Corrected_Value"
  
  
  test<-corrected %>%
    mutate(dummy=TRUE) %>%
    left_join(df2 %>% mutate(dummy=TRUE)) %>%
    filter(Name==Name, Iso==Iso, Condition==Condition,Exp==Exp) %>%
    select(-dummy)
  
  colnames(test)[grep("^Value$",colnames(test))]<-"Old Uncorrected Value"
  
  return(test)}

#test_version2<-test %>% spread(Exp,Corrected_Value)
#test <- test[c("Name", "Condition", "Iso","Norm_Av","Norm_Std","CV","Av","Nr.C.x","MID1","MID2","MID3","KEGG.ID","Nr.C.y","Exp","Value","Corrected_Value","ANOVA","Sig")] #reorder



