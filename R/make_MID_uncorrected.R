#' This function produces a data frame with mass isotopologue distribution (MID) data that is NOT corrected for naturally occurring 13C.
#' @author Daniel Braas, Thomas Graeber
#' @param DF The input data. This should be a data frame with Name, Iso, Condition, Exp, Value, Used_ID, KEGG.ID, Nr.C, Rt and Formula columns.
#' @return A data frame with MID data.
#' @export

make_MID_uncorrected <- function(DF){
  if (exists('Title')==F) stop('Title not specified')
  percent_label <- DF %>%
    select(Name, KEGG.ID, Condition, Iso, Nr.C, Exp, Value) %>%
    spread(Exp, Value)
# INPUT: Max carbons to calculate and #iterations (2 is min number)
  na =.01109
  MaxCarbonCalculate = 50
  #MaxIter = 5000
  MaxIter = 2
  
# Make metabolite name column factors (was not before in test set)
  #percent_label[,1] = as.factor(percent_label[,1])
  percent_label$Name <- factor(percent_label$Name)

  #change 0 values to NA
  percent_label[percent_label==0]=NA
  label_uncorr=percent_label
  
  MID_uncorr <- label_uncorr %>%
    gather(Exp, Value, -Name, -Condition, -Iso, -KEGG.ID, -Nr.C) %>%
    group_by(Name, Condition, Exp) %>%
    arrange(Name, Condition, Exp) %>%
    mutate(Sum = sum(Value, na.rm=T),
           Fraction = Value*100/Sum) %>%
    ungroup() %>%
    group_by(Name, Condition, Iso) %>%
    mutate(Norm_Av=mean(Fraction, na.rm=T),
           Norm_Std=sd(Fraction, na.rm=T),
           CV=sd(Fraction, na.rm=T)/mean(Fraction, na.rm=T),
           Av = Norm_Av) %>%
    ungroup() %>%
    select(Name, Condition, Iso, Exp, Fraction, Norm_Av, Norm_Std, CV, Av, Nr.C)
  MID_uncorr$Exp <- gsub('Exp','MID', MID_uncorr$Exp)
  
  #This part needs to be inserted at some point to account for NAs
  #test1 <- split(MID_uncorr, MID_uncorr[c('Iso','Condition', 'Name')])
  #NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  #new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Fraction)))
  
  MID_uncorr$Fraction[is.na(MID_uncorr$Fraction)] <- 0
  # Anova analysis instead of Ttest
  data8_uncorr=split(MID_uncorr, MID_uncorr[,c(3,1)], drop=TRUE)
  #ANOVA=sapply(data8_uncorr, function(x) {if(sum(x$Fraction, na.rm=T)==0) return(NA) else {anova(aov(x$Fraction~x$Condition))$Pr[1]}})
  ANOVA=suppressWarnings(sapply(data8_uncorr, function(x) anova(aov(x$Fraction~x$Condition))$Pr[1]))
  
  data3_uncorr <- MID_uncorr %>%
    spread(Exp, Fraction) %>%
    inner_join(label_uncorr, ., by = c("Name", "Condition", "Iso","Nr.C")) %>%
    arrange(Name, Iso)
  
  #add indicator of significance
  data3_uncorr$ANOVA=rep(ANOVA,each=length(unique(info$Condition)))
  for (i in 1:nrow(data3_uncorr)){
    if (data3_uncorr$ANOVA[i] == "NaN") data3_uncorr$Sig[i]=""
    else if (data3_uncorr$ANOVA[i] <= 0.001) data3_uncorr$Sig[i]="***"
    else if (data3_uncorr$ANOVA[i] <= 0.01) data3_uncorr$Sig[i]="**"
    else if (data3_uncorr$ANOVA[i] <= 0.05) data3_uncorr$Sig[i]="*"
    else data3_uncorr$Sig[i]=""
  }
  
  write.csv(data3_uncorr, file=paste0(Title,"-Isotopomer data-uncorrected.csv"), row.names=FALSE)
  save(data3_uncorr, file='MID_uncorr.rdata')
  
  return(data3_uncorr)
}
