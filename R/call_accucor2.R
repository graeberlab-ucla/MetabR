#' call_accucor2
#' #devtools::install_github("wangyujue23/AccuCor2",INSTALL_opts = "--no-multiarch")
#' Calls Accucor Natural Isotope Correction for DOUBLE labeled sets.
#' Iso column options are: C13N15-x-x, C13D-x-x, C12 PARENT, C13-, N15-, D-
#' Converts mid_output data to format for Accucor2.
#' @param mid_output Uncorrected Mid output for DOUBLE labeled set.
#' @param abbrev The downloaded Google sheets: Abbrev_New2. Columns Abb, Neutral.Formula,Compound, and KEGG.ID are used in function.
#' @param output_dir path
#' @param Title Title to use for naming the corrected file
#' @return corr3 NIC corrected version. Not yet in the format for "MIDS corrected output"
#' @importFrom dplyr select arrange mutate rename relocate group_by ungroup
#' @importFrom tidyr unite gather spread separate
#' @export
#'
#' @examples temp<-call_accucor2(mid_output, abbrev, output_dir, Title)
#'
call_accucor2<-function(mid_output, abbrev, output_dir, Title){

  # library("xlsx")
  # library(accucor2)
  # library(dplyr)
  # library(tidyr)
  # library(stringr)

  sheet_name<-"Sheet1"
  resolution=140000

  if(dim(subset(mid_output, startsWith(as.character(Iso), "C13N15")))[1]!=0){
    label= "CN"
    IsoID1="C13"
    IsoID2="N15"
    IsoIDs<-c("13C#", "15N#")

  } else if(dim(subset(mid_output, startsWith(as.character(Iso), "C13D")))[1]!=0 | dim(subset(mid_output, startsWith(as.character(Iso), "D")))[1]!=0){
    label= "CH"
    IsoID1="C13"
    IsoID2="D"
    IsoIDs<-c("13C#", "2H#")
  }

  uncorr<-mid_output%>%
    select(Name, Condition, Iso, which(grepl("Exp", names(mid_output)))) %>%
    gather(Exp, Value, -Name, -Condition, -Iso) %>%
    unite(Condition_Exp, c(Condition, Exp), sep='_') %>%
    mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
    spread(Condition_Exp, Value)

  y=as.data.frame(uncorr$Iso)
  y[,1]<-gsub("C12 PARENT|M0", "0", y[,1])
  y[,1]<-gsub(paste0("^", IsoID2, "-"), paste0(IsoID2, "-0-"), y[,1])

  isochart<-as.data.frame(stringr::str_split_fixed(y[,1], "\\-", 3))
  isochart <- replace(isochart, isochart=='', '0')
  isochart[,1]<-NULL

  names(isochart)<-IsoIDs
  isochart[,1]<-as.numeric(isochart[,1])
  isochart[,2]<-as.numeric(isochart[,2])

  uncorr<-cbind(uncorr, isochart)
  uncorr$Iso<-NULL

  abbrev_og<-abbrev
  abbrev<-abbrev[,c('Abb', 'Neutral.Formula')]
  names(abbrev)[1]<-"compound"
  names(abbrev)[2]<-"formula"
  abbrev$charge<- -1
  abbrev$charge<-as.integer(abbrev$charge)
  abbrev2<-na.omit(abbrev)

  setwd(output_dir)

  write.csv(abbrev2, paste0("Abbrev_Accucor2_input.csv"), row.names = F)

  sub<-uncorr[uncorr$Name %in% abbrev2$compound,]

  others1<-uncorr[!(uncorr$Name %in% abbrev2$compound),]

  print(paste(unique(others1$Name), "are not included in isotope correction as they don't match any metabs in Abbrev."))#check!!!

  sub$parent<-1
  sub$Expected<-1

  sub<-sub%>%select(Name, parent, which(grepl("#", names(sub))), Expected, everything() ) #grabs C,N, or H columns

  xlsx::write.xlsx(sub, paste0("Uncorr_Accucor2_input_",label, ".xlsx"), row.names = F)


  corrected <- accucor2::dual_correction(paste0("Uncorr_Accucor2_input_",label, ".xlsx"),
                               sheet_name,paste0("Abbrev_Accucor2_input",  ".csv"),label, Resolution = resolution) #Isotopes can be "CN" or "CH"

  filename<-list.files(pattern = paste0("Corrected_", gsub( '-', '', Sys.Date())))
  file.rename(filename, paste0(Title,"_Accucor2_function_output.xlsx"))

  corr<-corrected$Corrected

  corr[,2]<-paste0(stringr::str_pad(corr[,2], 2, pad = "0")) #C13
  corr[,3]<-paste0(stringr::str_pad(corr[,3], 2, pad = "0")) #N15 or H2
  names(corr)[2]<-"IsoID1"
  names(corr)[3]<-"IsoID2"

  corr2<-corr%>%gather(Condition_Exp, Value, -Compound, -IsoID1, -IsoID2) #add 2H


  corr2<-corr2%>%separate(Condition_Exp, c('Condition', 'Exp'), sep='_')

  names(corr2)[which(grepl("Compound",names(corr2)))]<-"Name"

  corr2 <- corr2 %>%
    group_by(Name, Condition, Exp) %>%
    arrange(Name, Condition, Exp) %>%
    mutate(Sum = sum(Value, na.rm = T), Fraction = Value * 100 / Sum) %>%
    ungroup()

  corr2$IsoID1[which(grepl("00", corr2$IsoID1) & grepl("00", corr2$IsoID2))] <-"C12 PARENT"
  corr2$IsoID2[which(grepl("C12 PARENT", corr2$IsoID1))]<-NA

  #C13-
  corr2$IsoID1[which(grepl("00", corr2$IsoID1))] <-paste0(paste0(IsoID2,"-"), corr2$IsoID2[which(grepl("00", corr2$IsoID1))])
  corr2$IsoID2[which(grepl(paste0(IsoID2,"-"), corr2$IsoID1))]<-NA

  #N15- or H2
  corr2$IsoID1[which(grepl("00", corr2$IsoID2))] <-paste0(paste0(IsoID1,"-"), corr2$IsoID1[which(grepl("00", corr2$IsoID2))])
  corr2$IsoID2[which(grepl(paste0(IsoID1,"-"), corr2$IsoID1))]<-NA

  #C13N15- or C13H2
  corr2$IsoID1[which(!is.na(corr2$IsoID2))]<-paste0(paste0(IsoID1,IsoID2,"-"), corr2$IsoID1[which(!is.na(corr2$IsoID2))], "-", corr2$IsoID2[which(!is.na(corr2$IsoID2))])
  corr2$IsoID2[which(!is.na(corr2$IsoID2))]<-NA
  corr2$IsoID2<-NULL
  names(corr2)[2]<-"Iso"


  corr2<-corr2%>%select(Name, Iso, Condition, Exp, Value)
  corr3<-corr2

  indx <- match(corr3$Name,abbrev_og$Abb)
  corr3$KEGG.ID <- abbrev_og$KEGG.ID[indx]

  corr3["Norm_Av"] <- corr3["Norm_Std"] <- corr3["CV"] <- corr3["Av"] <- corr3["Nr.C.x"] <-
    corr3["MID1"] <- corr3["MID2"] <- corr3["MID3"] <- corr3["Nr.C.y"] <- corr3["Old Uncorrected Value"] <-
    corr3["ANOVA"] <- corr3["Sig"] <- NA

  corr3 <- corr3 %>% relocate(KEGG.ID, .after = MID3)


  return(corr3)
}

