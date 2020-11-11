#' This function produces a data frame with mass isotopologue distribution (MID) data that is corrected for naturally occurring 13C.
#' @author Daniel Braas
#' @param DF The input data. This should be a data frame with Name, Iso, Condition, Exp, Value, Used_ID, KEGG.ID, Nr.C, Rt and Formula columns.
#' @param isotope_correction_remove_small_values_flag logical flag for wheter to run the remove small values code
#' @param isotope_correction_remove_small_values_threshold threshold used in the remove small values code
#' @return A data frame with MID data.
#' @export

make_MID2 <- function(DF,isotope_correction_remove_small_values_flag = FALSE, isotope_correction_remove_small_values_threshold = 0.1){                  # DF <- data2
  if (exists('Title')==F) stop('Title not specified')
  percent_label <- DF %>%
    select(Name, KEGG.ID, Condition, Iso, Nr.C, Exp, Value) %>%
    spread(Exp, Value)
  # INPUT: Max carbons to calculate and #iterations (2 is min number)
  na =.01109
  MaxCarbonCalculate = 50
  MaxIter = 5000
  
  # Make metabolite name column factors (was not before in test set)
  #percent_label[,1] = as.factor(percent_label[,1])
  percent_label$Name <- factor(percent_label$Name)
  
  #Duplicate data frame to input corrected values into
  percent_label_corrected=percent_label
  
  #Make empty data frame for abundance correction parameters (calculated abundance, and %error)
  correctiondata=data.matrix(matrix(0,nrow(percent_label),16))
  colnames(correctiondata)=c("Metabolite",
                             "rep1 corrected","rep2 corrected","rep3 corrected",
                             "rep1 na added","rep2 na added","rep3 na added",
                             "rep 1 %error","rep 2 %error","rep 3 %error",
                             "rep 1 numiter","rep 2 numiter","rep 3 numiter",
                             "rep 1 errorminfromiter","rep 2 errorminfromiter","rep 3 errorminfromiter")
  
  correctiondata[,1]=paste(as.character(percent_label[,1]),as.character(percent_label[,3]),as.character(percent_label[,4]))
  
  #iterate through conditions, assumes conditions is in column 3
  for(c in 1:length(levels(percent_label[,3]))){
    
    #iterate through metabolites, assumes metabolites are in column 1
    for(n in 1:length(levels(percent_label[,1]))){
      
      #see what metabolite you are currently calculating
      print(levels(percent_label[,1])[n])
      
      #find the rows that satisfy condition c and metabolite n
      index=which((percent_label[,3]==levels(percent_label[,3])[c])&((percent_label[,1]==levels(percent_label[,1])[n])))
      
      #Only continue if there are datapoints that satisfy condition c and metabolite n
      if(length(index)!=0){
        
        #read in the condition + metabolite into matrix temp
        temp=percent_label[index,]
        #find Cmax for the metabolite
        Cmax=temp[1,5]
        
        #limit the max number of carbons to calculate.
        if(Cmax>MaxCarbonCalculate){
          Cmax=MaxCarbonCalculate
        }
        
        #Iterate correction across all replicates starting from column 6 where replicate 1 is supposed to be.
        for(z in 6:ncol(temp)){
          
          #Read natural abundance isotopomers of replicate z from temp and set NA to 0
          Iob=rep(0,Cmax+1)
          Iob[1:nrow(temp)]=temp[1:nrow(temp),z]
          Iob[is.na(Iob)]=0
          
          #Set observed = Ina(Ina will be Iob updated with Icalc values where Iob=0 later)
          Ina=Iob
          
          #variables to initalize
          iter=0
          error=1e11
          errornext=1e10
          
          I=rep(0,Cmax+1)
          Iprev=rep(0,Cmax+1)
          Icalcprev=rep(0,Cmax+1)
          Icalc=rep(0,Cmax+1)
          erroritermin=0
          
          #Re-update Ina with Iob+Icalc(where Iob=0) as long as error is decreasing (Iob-Icalc)
          while(errornext<error&&iter<=MaxIter){
            
            if(iter>=2){
              
              erroritermin=error-errornext+erroritermin
              
            }
            
            #save previous rounds values
            error=errornext
            Iprev=I
            Icalcprev=Icalc
            
            #calculate error minimized from iterations
            
            
            #Subtract natural abundance Ina-->I
            for(i in 0:Cmax){
              
              #initialize toright(carbons that are naturally labeled from the true isotopomer and decrease the measured amount) and fromleft (added
              #from natural abundance labeling of lower isotope, increases the measured amount)
              
              toright = rep(0,Cmax+1)
              fromleft= rep(0,Cmax+1)
              
              if((i+1)<=Cmax){
                for(k in (i+1):Cmax){
                  toright[k+1]=(choose(Cmax-i,k-i))*(1-na)^(Cmax-k)*na^(k-i)
                }
              } else {
                toright[k+1]=0
              }
              x=0
              while(x<i){
                fromleft[x+1]=I[x+1]*(choose(Cmax-x,i-x))*(1-na)^(Cmax-i)*na^(i-x)
                x=x+1
              }
              I[i+1]=0
              I[i+1]=(Ina[i+1]-sum(fromleft))/(1-sum(toright))
              
            }
            
            #Flatten  + renormalize I to Iob
            I[I<0]=0;
            I=I/(sum(I))*sum(Iob[1:(Cmax+1)]);
            
            #Add back natural abundance I-->Icalc
            Icalc=rep(0,Cmax+1)
            for(i in 0:Cmax){
              
              toright = rep(0,Cmax+1)
              fromleft= rep(0,Cmax+1)
              
              if((i+1)<=Cmax){
                for(k in (i+1):Cmax){
                  toright[k+1]=(choose(Cmax-i,k-i))*(1-na)^(Cmax-k)*na^(k-i)
                }
              } else {
                toright[k+1]=0
              }
              x=0
              while(x<i){
                fromleft[x+1]=I[x+1]*(choose(Cmax-x,i-x))*(1-na)^(Cmax-i)*na^(i-x)
                x=x+1
              }
              
              Icalc[i+1]=I[i+1]*(1-sum(toright))+sum(fromleft)
            }
            
            #Calculate error
            errornext=sum(abs(Iob[1:(Cmax+1)]-Icalc))
            
            #update Ina with Iob and Icalc for 0 values in Iob
            Itemp=Iob
            replacement=which(Iob[1:(Cmax+1)]==0)
            Itemp[replacement]=Icalc[replacement]
            
            if(is.na(sum(Itemp))){
              Itemp=0
            }
            
            if(sum(Itemp)!=0){
              Itemp=Itemp/(sum(Itemp))*sum(Iob[1:(Cmax+1)])
            }
            
            Ina=Itemp
            if (iter !=0 & iter/1000 == round(iter/1000)) {print(iter)}
            iter=iter+1;
            
            #if error =NA b/c zero values, then stop iterating
            if(is.na(errornext)){
              errornext=2*error
            }
          }
          
          #Update values depending on #measurements > or < Cmax
          if(length(index)<Cmax){
            start=index[1]
            end=index[length(index)]
          } else {
            start=index[1]
            end=index[1]+Cmax
          }
          
          
          #input error and re-calculated natural abundance MID into correction data
          #correctiondata[start:end,z+3]=as.numeric(error/sum(Iob[1:(Cmax+1)]))*100
          #correctiondata[start:end,z]=Icalcprev[1:(end-start+1)]
          #correctiondata[start:end,z+6]=iter
          #correctiondata[start:end,z+9]=as.numeric(erroritermin/sum(Iob[1:(Cmax+1)]))*100
          
          #update new values into isotopomer matrix
          percent_label_corrected[start:end,z]=Iprev[1:(end-start+1)]
        }
      }
    }
  }
  
  
  #change 0 values to NA
  percent_label_corrected[percent_label_corrected==0]=NA
  label=percent_label_corrected
  
  #removing small values that meet three criteria
  #                     1) not MO values, and 
  #                     2) the value is below the inputted threshold pqrameter (i.e., only small values removed), and 
  #                     3) the original uncorrected data had no measurement (refered to as 'the no uncorr measuremnet logical test' below)
  #since the critria involve both corrected and uncorrected data, need to create a merged datarame
  label_uncorr=percent_label
  colnames(label_uncorr) <- gsub("^Exp(\\d+)$", "uncorr.Exp\\1", colnames(label_uncorr))
  label_uncorr$KEGG.ID <- NULL
  label_uncorr$Nr.C <- NULL
  label.temp <- merge(label,label_uncorr,by=c("Name","Condition","Iso"))
  #debugging   label.temp <- label.temp[label.temp$Name=="Malonyl-CoA",]
  
  #step 1 of 3 - reorganizing the uncorrected data, defining logical test result based on if the uncorrected data had a value
  MID.temp <- label.temp %>%
    gather(paired.Exp, Value, starts_with("uncorr.Exp")) %>%
    group_by(Name, Condition, paired.Exp) %>%
    arrange(Name, Condition, paired.Exp) %>%
    mutate(test = !(is.na(Value) | Value==0), uncorr=Value) %>%   # test = "original uncorr not empty/zero"
    ungroup() %>%
    select(Name, Condition, Iso, paired.Exp, test, Nr.C, starts_with("Exp"))
    #select(Name, Condition, Iso, paired.Exp, test, uncorr, Nr.C, starts_with("Exp"))
   
  mids <- colnames(MID.temp)[grepl("^Exp\\d+$", colnames(MID.temp))]
  n.replicates = length(mids)
  
  #step 2 - reorganizing the corrected data, calculating statistical metrics, and utilizing the no uncorr measuremnet logical test from step 1
  #removing small values that meet the three criteria summarized above 
  #debugging   MID.temp <- MID.temp[MID.temp$Name=="Malonyl-CoA",]    isotope_correction_remove_small_values_threshold=0.1
  #MID.temp <- MID.temp[MID.temp$Name=="HMG-CoA",]
  MID <- MID.temp %>%
    gather(Exp, Value, starts_with("Exp")) %>%
    group_by(Name, Condition, Exp) %>%
    arrange(Name, Condition, Exp) %>%
    mutate(Sum = sum(Value, na.rm=T)/n.replicates,
           #Fraction = Value*100/Sum) %>%
           #Fraction = ifelse((Iso=="M0" | test | format(Value*100/Sum, scientific=F)>=as.numeric(isotope_correction_remove_small_values_threshold)), # test = "original uncorr not empty/zero"
           #                  Value*100/Sum,NA)) %>%
           Fraction = ifelse((Iso=="M0" | test | Value*100/Sum>=as.numeric(isotope_correction_remove_small_values_threshold)), # test = "original uncorr not empty/zero"
                             Value*100/Sum,NA)) %>%
    ungroup() %>%
    group_by(Name, Condition, Iso) %>%
    mutate(Norm_Av=mean(Fraction, na.rm=T),
           Norm_Std=sd(Fraction, na.rm=T),
           CV=sd(Fraction, na.rm=T)/mean(Fraction, na.rm=T),
           Av = Norm_Av) %>%
    ungroup() %>%
    select(Name, Condition, Iso, Exp, paired.Exp, Fraction, Norm_Av, Norm_Std, CV, Av, Nr.C, starts_with("uncorr.Exp"))
    #select(Name, Condition, Iso, Exp, paired.Exp, test, Fraction, Norm_Av, Norm_Std, CV, Av, Nr.C, starts_with("uncorr.Exp"))
    #debug MID[(MID$Name=="Malonyl-CoA" & MID$Iso=="M4" & (MID$Exp=="MID1" | MID$Exp=="Exp1") & MID$paired.Exp=="uncorr.Exp1"),]
    #debug MID_ver1[(MID_ver1$Name=="Malonyl-CoA" & MID_ver1$Iso=="M4" & (MID_ver1$Exp=="MID1" | MID_ver1$Exp=="Exp1")),]
  
  #step 3 - cleaning up, discarding unneeded rows
  # MIDpre <- MID
  MID$paired.Exp = gsub("uncorr.", "", MID$paired.Exp)
  MID <- MID[MID$Exp == MID$paired.Exp,]
  
  MID$paired.Exp <- NULL
  MID$Exp <- gsub('Exp','MID', MID$Exp)
  
  
  #This part needs to be inserted at some point to account for NAs
  #test1 <- split(MID, MID[c('Iso','Condition', 'Name')])
  #NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  #new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Fraction)))
  
  MID$Fraction[is.na(MID$Fraction)] <- 0
  # Anova analysis instead of Ttest
  #data8=split(MID, MID[,c(3,1)], drop=TRUE)
  data8=split(MID, MID[,c("Iso","Name")], drop=TRUE)
  #ANOVA=sapply(data8, function(x) {if(sum(x$Fraction, na.rm=T)==0) return(NA) else {anova(aov(x$Fraction~x$Condition))$Pr[1]}})
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Fraction~x$Condition))$Pr[1]))
  
  data3 <- MID %>%
    spread(Exp, Fraction) %>%
    inner_join(label, ., by = c("Name", "Condition", "Iso","Nr.C")) %>%
    arrange(Name, Iso)
  
  #add indicator of significance
  data3$ANOVA=rep(ANOVA,each=length(unique(info$Condition)))
  for (i in 1:nrow(data3)){
    if (data3$ANOVA[i] == "NaN") data3$Sig[i]=""
    else if (data3$ANOVA[i] <= 0.001) data3$Sig[i]="***"
    else if (data3$ANOVA[i] <= 0.01) data3$Sig[i]="**"
    else if (data3$ANOVA[i] <= 0.05) data3$Sig[i]="*"
    else data3$Sig[i]=""
  }
  
  write.csv(data3, file=paste0(Title,"-Isotopomer data.csv"), row.names=FALSE)
  save(data3, file='MID.rdata')
  return(data3)
}
