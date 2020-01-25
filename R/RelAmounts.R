#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param a data frame containing sample data, an optional anova parameter (set to TRUE if ANOVA is to be calculated), and an optional type parameter (set to '250K').
#'
#' @return a data frame that is in the correct format for the bar_update_manual function to create the Relative Amounts plot.
#'
#' @examples my_make_RelAmounts_QC(data2, anova = TRUE, type = '250K')
#'
#' @export

RelAmounts <- function(DF, anova = F, type = 'NULL')
{
  if (exists('Title')==F) stop('Title not specified')
  data4 <- DF %>%
    select(Name, Condition, Iso, KEGG.ID, Exp, Value) %>%
    mutate(Exp = factor(Exp, levels = unique(Exp))) %>%
    group_by(Name, Condition, Exp) %>%
    mutate(Amount=sum(Value, na.rm=T)) %>%
    ungroup() %>%
    filter(Iso=='M0') %>%
    select(Name, Condition, KEGG.ID, Exp, Amount) %>%
    spread(Exp, Amount)

  ATP_ADP=try(cbind(Name="ADP/ATP",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                    data4[data4$Name=="ADP",4:length(data4)]/data4[data4$Name=="ATP",4:length(data4)]),
              silent=T)
  if (exists('ATP_ADP')==T & class(ATP_ADP) != 'try-error') data4 <- rbind(data4, ATP_ADP)
  ATP_AMP=try(cbind(Name="AMP/ATP",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                    data4[data4$Name=="AMP",4:length(data4)]/data4[data4$Name=="ATP",4:length(data4)]),
              silent = T)
  if (exists('ATP_AMP')==T & class(ATP_AMP) != 'try-error') data4 <- rbind(data4, ATP_AMP)
  GSH_GSSG=try(cbind(Name="GSH/GSSG",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                     data4[data4$Name=="GSH",4:length(data4)]/data4[data4$Name=="GSSG",4:length(data4)]),
               silent = T)
  if (exists('GSH_GSSG')==T & class(GSH_GSSG) != 'try-error') data4 <- rbind(data4, GSH_GSSG)
  Creatine_PCreatine=try(cbind(Name="Creatine/P-Creatine",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                               data4[data4$Name=="Creatine",4:length(data4)]/data4[data4$Name=="P-Creatine",4:length(data4)]),
                         silent = T)
  if (exists('Creatine_PCreatine')==T & class(Creatine_PCreatine) != 'try-error') data4 <- rbind(data4, Creatine_PCreatine)

  data4 <- data4 %>%
    gather(Exp, Amount, -Name, -Condition,-KEGG.ID)
  data4$Amount <- suppressMessages(mapvalues(data4$Amount, c('Inf','NaN'), c(NA,NA)))
  data4$Amount[data4$Amount==0] <- NA
  test1 <- split(data4, data4[c('Name','Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Amount)))

  data4 <- suppressWarnings(data4 %>%
                              arrange(Condition, Name) %>%
                              mutate(Amount=new.Value) %>%
                              group_by(Name, Condition) %>%
                              mutate(Av=mean(Amount, na.rm=T),
                                     Std=sd(Amount, na.rm=T),
                                     CV=sd(Amount, na.rm=T)/mean(Amount, na.rm=T)) %>%
                              ungroup())

  if(anova)
  {
    data8=split(data4, data4[,1])
    ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
    ANOVA=rep(ANOVA,1,each=length(levels(data4$Condition)))

    data4 <- data4 %>%
      select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
      spread(Exp, Amount) %>%
      arrange(Name, Condition) %>%
      cbind('ANOVA'=ANOVA, 'Sig'=NA) %>%
      group_by(Name) %>%
      mutate(Norm_Av=Av/Av[1],
             Norm_Std=Std/Av[1]) %>%
      ungroup()

    for (i in 1:nrow(data4))
    {
      if (is.na(data4$ANOVA[i])==T) data4$Sig[i]=""
      else if (data4$ANOVA[i] <= 0.001) data4$Sig[i]="***"
      else if (data4$ANOVA[i] <= 0.01) data4$Sig[i]="**"
      else if (data4$ANOVA[i] <= 0.05) data4$Sig[i]="*"
      else data4$Sig[i]=""
    }
  }
  else if(type == '250K')
  {
    data4 <- data4 %>%
      select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
      spread(Exp, Amount) %>%
      arrange(Name, Condition)
    temp <- split(data4, data4['Name'])

    Rel_function <- function(x)
    {
      if(x[1] == 0)
        x = x / 1
      else
        x = x / x[1]
    }
    Rel_Av <- as.vector(sapply(temp, function(x) Rel_function(x$Av)))
    Std_Av <- as.vector(sapply(temp, function(x) Rel_function(x$Std)))

    data4 <- data4 %>% mutate(RelAmounts_Ave = Rel_Av, RelAmounts_Std = Std_Av)
  }
  else
  {
    data4 <- data4 %>%
      select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
      spread(Exp, Amount) %>%
      arrange(Name, Condition) %>%
      group_by(Name) %>%
      mutate(RelAmounts_Ave = Av/Av[1],
             RelAmounts_Std = Std/Av[1]) %>%
      ungroup()
  }

  return(data4)
}