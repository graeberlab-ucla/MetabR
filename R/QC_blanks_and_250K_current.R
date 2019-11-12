library(MetabR)
library(MetabFUN)
library(tidyr)
library(readxl)
library(dplyr)
#if using new PCA Plot V2
library(ggrepel)
library(RColorBrewer)

data_dir <- "Y:/4.0/Projects/Tissue/Scafoglio Lab/PS-10312019_10-23_Vanq"
output_dir <- "Y:/4.0/Projects/Tissue/Scafoglio Lab/PS-10312019_10-23_Vanq"
anno_dir <- "C:/Users/djtan/Documents/metabolobics scripts"
unlabelled <- TRUE  #TRUE if experiment is unlabelled and FALSE if it is labelled
ifelse (grepl("ICS",data_dir),ICS<-TRUE,ICS<-FALSE)
ifelse (grepl("Medium", data_dir), Medium <- TRUE, Medium <- FALSE)

#Annotation Data
setwd(anno_dir)
if (ICS)
{
  abbrev <- read.csv('abbreviations_ICS.csv', header = T)
} else
{
  abbrev <- read.csv('abbreviations2.csv', header = T)
}

setwd(data_dir)
Title <- gsub('.xls[x]?','', list.files(pattern = '.xls[x]?'))
Title <- paste("QC-Blank Analysis", Title)
Title <- gsub('_sample info|_sample info sheet', '', Title) 
#adjust title if necessary

#make sure sample info sheet is updated to include blanks and 250Ks

# enter cell names and cell numbers
#info <- read_excel(list.files()[grep('.xls[x]?',list.files())])
info <- read_excel(list.files()[grep('.xls[x]?',list.files())])
info$Condition <- gsub("^\\s+|\\s+$", "", info$Condition)        #remove leading or trailing white space
info$Condition <- gsub('/', '-', info$Condition)
info$Condition <- gsub('_', '-', info$Condition)
num_conditions <- length(unique(info$Condition))
# reorder sample info sheet if necessary
info$Condition <- factor(info$Condition, levels = unique(info$Condition))    #important when rearranging the order of samples

Freq <- data.frame(table(info$Condition))

Sample.Name <- vector()
for (i in 1:length(levels(info$Condition))){
  for (j in 1:Freq$Freq[i]){
    Sample.Name <- append(Sample.Name, paste(Freq$Var1[i], j, sep='_'))
  }
}
if (is.numeric(info$Cell.Number)==F) print('Problem: You need to convert the Cell.Number column')
#warning is okay if not normalizing


#make sure sample.name ordering matches info ordering (fix if it doesn't)
info$Sample.Name <- Sample.Name
samples <- info
numbers <- samples$Cell.Number

#Read in tracefinder metabolite peak area data
folder.name <- gsub('(.)*Projects/','', data_dir) %>% gsub('/','-',.)
if (ICS)
{
  files <- "ICS/ICS.csv"
} else if (Medium) {
  files<-"Medium/Medium.csv"
}else
{
  files <- c("AA/AA.csv", "CoAs/CoAs.csv", "Currency/Currency.csv", "dNTPs/dNTPs.csv", "FA/FA.csv", "Glycolysis/Glycolysis.csv", "Hexosamine/Hexosamine.csv", "NTPs/NTPs.csv", "PPP/PPP.csv", "TCA/TCA.csv","Urea/Urea.csv")
  #files <- c("AA/AA-2.csv", "CoAs/CoAs-2.csv", "Currency/Currency-2.csv", "dNTPs/dNTPs-2.csv", "FA/FA-2.csv", "Glycolysis/Glycolysis-2.csv", "Hexosamine/Hexosamine-2.csv", "NTPs/NTPs-2.csv", "PPP/PPP-2.csv", "TCA/TCA-2.csv","Urea/Urea-2.csv")
}

# files <- list.files(pattern='.csv', recursive=T)
# files <- files[!grepl('Suggested|Retention time', files)]
# files <- files[grep('(.)*/.(.)*.csv', files)]
#check files for only the 11 tracefinder output files
#or the one ICS/ICS.csv for ICS experiments


if (sum(grepl('ICS',files)) > 0) Title <- paste0('ICS-', Title)
if (sum(grepl('[Mm]edium|[Mm]edia',files)) > 0) Title <- paste0('Footprint-',Title)

data <- lapply(files, read.csv, header=T, na.strings = c("N/F",'N/A'))
data <- do.call(rbind, data)
###use rbind.fill from the plyr package (load in plyr first and then dplyr if necessary) if the method/batch files have differing columns, make sure the dropped columns are not important or relvant to downstream data
#data=do.call(rbind.fill, data)
data$Experiment <- folder.name

setwd(output_dir)
write.csv(data, 'all data raw.csv', row.names=F)

data <- data %>%
  dplyr::select(Compound, Filename, Area) %>%
  spread(Filename, Area)


data$Compound = as.character(data$Compound)
for (i in 1:length(data$Compound)){                   #loop through IDs and fill up with Ms
  if (grepl(" M[0-9]+$", data$Compound[i])==FALSE) data$Compound[i]=paste(data$Compound[i],"M0")}
Used_ID <- gsub(" M[0-9]+$| Std", "", data$Compound)
name_split <- regexpr("M[0-9]+$|Std",data$Compound)
Iso <-regmatches(data$Compound, name_split)
data <- data.frame(Used_ID, Iso, data[,2:length(data)])

#data <- filter(data, !(grepl('Std', data$Iso)))

data$Used_ID <- as.character(data$Used_ID)
data$Iso <- factor(data$Iso, levels=paste0('M',0:50))
data <- arrange(data, Used_ID, Iso)

#make sure 250FKs and blanks are added onto the info sheet and they are order the same as in the raw data
###
new_colnames<-info$Sample
non_blanks_names<-new_colnames[!grepl(paste0("QC",collapse="|"),new_colnames)]
blanks_names<-new_colnames[grepl(paste0("QC",collapse="|"),new_colnames)]
col.order<-c("Used_ID","Iso",non_blanks_names,blanks_names)
col.order<-gsub("-",".",col.order)
data<-data[,col.order]

# check ordering
###
for (i in 1:length(Sample.Name)) colnames(data)[i+2]=info$Sample.Name[i]

Norv <- data %>%
  filter(Used_ID=='Norvaline') %>%
  gather(Sample, Norv,-Used_ID, -Iso) %>%
  .$Norv

if (length(Norv)==0){
  Norv <- rep(1,nrow(samples))
  print("No norvaline")
}
norm <- data

norm <- norm[!norm$Used_ID=="Norvaline",] 

data2 <- norm %>%
  gather(Exp, Value, -Used_ID, -Iso) %>%
  separate(Exp, c('Condition', 'Exp'), sep='_') %>%
  left_join(., Abbrev, by='Used_ID') %>% 
  dplyr::rename(Name = Abb) %>% 
  dplyr::select(Name, Iso, Condition, Exp, Value, KEGG.ID, Nr.C, everything())

data2$Condition <- factor(data2$Condition, levels=as.character(unique(samples$Condition)))
data2$Exp <- paste0('Exp', data2$Exp)
NA_Names <- data2 %>%
  dplyr::group_by(Name) %>%
  dplyr::mutate(NA_list=sum(is.na(Value)),
         NA_potential=dplyr::n()) %>%
  filter(NA_list==NA_potential) %>%
  distinct()
data2 <- data2 %>%
  filter(!(Name %in% NA_Names$Name))
data2 <- data2[!is.na(data2$Name),]
data2$Name <- as.character(data2$Name)

#if analysis is unlabeled make sure to filter for only M0 and M1 isotopomers
if(unlabelled)
  data2 <- subset(data2, Iso == 'M0' | Iso == 'M1')

##########
write.csv(data2, file=paste0(Title,"-all data unnormalized plus blanks and 250K.csv"), row.names=FALSE)

select <- dplyr::select
data4 <- my_make_RelAmounts_QC(data2) 

#condition order
info <- info[order(info$Run.Order),]
condition_order <- info$Condition
condition_order <- as.vector(unique(condition_order))

#change the factor level of the data4$Conditions to match the new condition_order 
data4$Condition <- factor(data4$Condition, levels = condition_order)

index_blank<-grep("blank",unique(info$Condition))
index_250k<-grep("250K",unique(info$Condition))
indices<-list(index_250k,index_blank)
names(indices)<-c("250k","blank")

#load in bar_update_manual and bar_plot_update_manual from manual x axis barplots
plotname=paste(Title,"-RelAmounts plots unnormalized plus blanks and 250K.pdf", sep='')
pdf(file = plotname, width=14, height=10, pointsize=12)
bar_update_manual('glycolysis', data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("TCA",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("PPP",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Curr",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Cys",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("AA",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("FA",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Hex",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Adenine",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Cytosine",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Guanine",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Thymine",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual("Uracil",data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual('Fru',data4, n = num_conditions, type = '250k',index=indices)
bar_update_manual('CoAs',data4, n = num_conditions, type = '250k',index=indices)
dev.off()
