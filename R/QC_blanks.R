library(devtools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("KEGGREST")
#install.packages("tidyverse")
#install_github("danielbraas/MetabFUN")
#install_github("juyeki/MetabR")
library(MetabFUN)
library(MetabR)
library(tidyr)
library(readxl)
library(dplyr)
#if using new PCA Plot V2
library(ggrepel)

data_dir <- "Y:/4.0/Projects/Tissue/Scafoglio Lab/PS-10312019_10-23_Vanq"
output_dir <- "Y:/4.0/Projects/Tissue/Scafoglio Lab/PS-10312019_10-23_Vanq"
anno_dir <- "C:/Users/djtan/Documents/metabolobics scripts"
unlabelled <- TRUE
ICS <- ifelse (grepl("ICS",data_dir), ICS <- TRUE, ICS <- FALSE)
ifelse (grepl("Medium", data_dir), Medium <- TRUE, Medium <- FALSE)

#Annotation data
setwd(anno_dir)
if (ICS)
{
  abbrev <- read.csv('abbreviations_ICS.csv', header = T)
} else
{
  abbrev <- read.csv('abbreviations2.csv', header = T)
}

#Title
setwd(data_dir)
Title <- gsub('.xls[x]?','', list.files(pattern='.xls[x]?'))
Title <- gsub('_sample info sheet|_sample info', '', Title)

#Info
info <- read_excel(list.files()[grep('.xls[x]?',list.files())])
info$Condition <- gsub("^\\s+|\\s+$", "", info$Condition)      
info$Condition <- gsub('/', '-', info$Condition)
info$Condition <- gsub('_', '-', info$Condition)

#Modifying info 
if (ICS)
{
  info <- info[-grep(pattern = "QC-50K", x = info$Sample),]
} else
{
  info <- info[-grep(pattern = "QC-250K", x = info$Sample),]
}
info$Condition[grep(pattern = "QC-blank", x = info$Sample)] <- "blank"
#info <- info[-22,]
info <- info[!grepl(pattern = "QC-0", x = info$Sample),] #should remove rest of non-blank QC
num_conditions <- length(unique(info$Condition))

#Adding sample names to info
info$Condition <- factor(info$Condition, levels = unique(info$Condition))
freq <- data.frame(table(info$Condition))
sample.name <- vector()
for (i in 1:length(unique(info$Condition))){
  for (j in 1:freq$Freq[i])
  {
    sample.name <- append(sample.name, paste(freq$Var1[i], j, sep='_'))
  }
}
info$sample.name <- sample.name
samples <- info

#Tracefinder Metabolite Data
setwd(data_dir)
folder.name <- gsub('(.)*Projects/','',getwd()) %>% gsub('/','-',.)
if (ICS)
{
  tracefinder_files <- "ICS/ICS.csv"
} else if (Medium) {
  tracefinder_files<-"Medium/Medium.csv"
}else
{
  tracefinder_files <- c("AA/AA.csv", "CoAs/CoAs.csv", "Currency/Currency.csv", "dNTPs/dNTPs.csv", "FA/FA.csv", "Glycolysis/Glycolysis.csv", "Hexosamine/Hexosamine.csv", "NTPs/NTPs.csv", "PPP/PPP.csv", "TCA/TCA.csv","Urea/Urea.csv")
  #tracefinder_files <- c("AA/AA-2.csv", "CoAs/CoAs-2.csv", "Currency/Currency-2.csv", "dNTPs/dNTPs-2.csv", "FA/FA-2.csv", "Glycolysis/Glycolysis-2.csv", "Hexosamine/Hexosamine-2.csv", "NTPs/NTPs-2.csv", "PPP/PPP-2.csv", "TCA/TCA-2.csv","Urea/Urea-2.csv")
}

data <- lapply(tracefinder_files, read.csv, header=T, na.strings = c("N/F",'N/A'))
data <- do.call(rbind, data)
data$Experiment <- folder.name
data <- data[!grepl("QC-250K", x = data$Filename), ] #Remove 250ks from Metabolite Data
data <- data[!grepl("QC-50K", x = data$Filename), ] #Remove 250ks from Metabolite Data (if ICS) 
data <- data %>% dplyr::select(Compound, Filename, Area) %>% spread(Filename, Area)
#data$`QC-001` <- NULL
data <- data[, (!grepl("QC-0", x = colnames(data)))] #ADDED THIS 10/9: should remove rest of non-blank QCs


#make sure the blanks are the last conditions, move the columns if necesary

#Adding M0's to compound if it is missing an 'M'
data$Compound <- as.character(data$Compound)
for (i in 1:length(data$Compound))
{                  
  if (grepl(" M[0-9]+$", data$Compound[i])==FALSE) 
    data$Compound[i] <- paste(data$Compound[i],"M0")
}
#splitting Compound column
Used_ID <- gsub(" M[0-9]+$| Std", "", data$Compound) #removes M's
Iso <-regmatches(data$Compound, regexpr("M[0-9]+$|Std", data$Compound))
data <- data.frame(Used_ID, Iso, data[,2:length(data)])

# #matching sample ordering of info and data ##not necessary for this script
# info$Sample <- gsub('-', '.', info$Sample)
#info <- info[match(colnames(data)[3:ncol(data)], info$Sample), ] #info changed to match data 
#info <- info[c(25:27, 1:24),]

#Remove Standards from the data
data <- filter(data, !(grepl('Std', data$Iso)))
#Reordering data by M's
data$Used_ID <- as.character(data$Used_ID)
data$Iso <- factor(data$Iso, levels=paste0('M',0:50))
data <- arrange(data, Used_ID, Iso)
#adding sample names to data
####
new_colnames<-info$Sample
non_blanks_names<-new_colnames[!grepl(paste0("blank",collapse="|"),new_colnames)]
blanks_names<-new_colnames[grepl(paste0("blank",collapse="|"),new_colnames)]
col.order<-c("Used_ID","Iso",non_blanks_names,blanks_names)
col.order<-gsub("-",".",col.order)
data<-data[,col.order]
# check ordering!!
####
colnames(data) <- append(colnames(data[1:2]), info$sample.name)

#Norvaline
Norv <- data %>% filter(Used_ID=='Norvaline') %>% gather(Sample, Norv,-Used_ID, -Iso) %>% .$Norv
if (length(Norv) == 0)
{
  Norv <- rep(1,nrow(samples))
  print("No norvaline")
}
norm <- data[!data$Used_ID=="Norvaline",] 

#Data 2
data2 <- norm %>%
  gather(Exp, Value, -Used_ID, -Iso) %>%
  separate(Exp, c('Condition', 'Exp'), sep='_') %>%
  left_join(., abbrev, by = 'Used_ID') %>% 
  dplyr::rename(Name = Abb) %>% 
  dplyr::select(Name, Iso, Condition, Exp, Value, KEGG.ID, Nr.C, everything())

data2$Condition <- factor(data2$Condition, levels=as.character(unique(info$Condition)))
data2$Exp <- paste0('Exp', data2$Exp)
NA_Names <- data2 %>%
  dplyr::group_by(Name) %>%
  dplyr::mutate(NA_list=sum(is.na(Value)),
                NA_potential=dplyr::n()) %>%
  filter(NA_list==NA_potential) %>%
  distinct()

data2 <- data2 %>% filter(!(Name %in% NA_Names$Name))
data2$Name <- as.character(data2$Name)
data2 <- data2 %>% drop_na(Name)

#Filtering if experiment is unlabelled
if(unlabelled)
  data2 <- subset(data2, Iso == 'M0' | Iso == 'M1')
setwd(output_dir)
write.csv(data2, file = paste0(Title, "-all data unnormalized plus blanks.csv"), row.names=FALSE)

#load my_make_RelAmounts_QC to workspace
select <- dplyr::select
data4 <- my_make_RelAmounts_QC(data2)
index_blank<-grep("blank",unique(info$Condition))
index_250k<-grep("250K",unique(info$Condition))
indices<-list(index_250k,index_blank)
names(indices)<-c("250k","blank")

setwd(output_dir)
plotname=paste(Title,"-RelAmounts plots unnormalized plus blanks.pdf", sep='')
pdf(file = plotname, width=14, height=10, pointsize=12)
bar_update_manual('glycolysis',data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("TCA",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("PPP",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Curr",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Cys",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("AA",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("FA",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Hex",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Adenine",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Cytosine",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Guanine",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Thymine",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual("Uracil",data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual('Fru',data4, n = num_conditions, type = "blank",index=indices)
bar_update_manual('CoAs',data4, n = num_conditions, type = "blank",index=indices)

dev.off()


