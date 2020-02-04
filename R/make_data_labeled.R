#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param The dataframe containing the output from the MIDs function.
#'
#' @return A data frame with percent labeled calculations.
#'
#' @examples make_data_labeled(DF = data3)
#'
#' @export
#'
#'

make_data_labeled <- function(DF)
{
  if (exists('Title')==F) stop('Title not specified')
  data_labeled <- DF %>%
    select(Name, KEGG.ID, Condition, Iso, starts_with('MID')) %>%
    filter(Iso=='M0') %>%
    gather(Exp, Value, -Name, -KEGG.ID, -Condition, -Iso) %>%
    mutate(Labeled=(1-Value/100)*100) %>%
    group_by(Name, Condition) %>%
    mutate(Norm_Av=mean(Labeled, na.rm=T),
           Norm_Std=sd(Labeled, na.rm=T),
           CV=Norm_Std/Norm_Av,
           Av = Norm_Av) %>%
    select(-Iso, -Value) %>%
    ungroup()
  data_labeled$Exp <- gsub('MID','Labeled', data_labeled$Exp)
  test1 <- split(data_labeled, data_labeled[c('Name', 'Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Labeled)))

  data_labeled <- data_labeled %>%
    arrange(Condition, Name) %>%
    mutate(Labeled=new.Value)

  data8=split(data_labeled, data_labeled[,1])
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Labeled~x$Condition))$Pr[1]))
  ANOVA=rep(ANOVA,1,each=length(unique(info$Condition)))

  data_labeled <- spread(data_labeled, Exp, Labeled) %>%
    arrange(Name)

  data_labeled <- cbind(data_labeled, ANOVA, 'Sig'=NA)
  for (i in 1:nrow(data_labeled)){
    if (data_labeled$ANOVA[i] == "NaN") data_labeled$Sig[i]=""
    else if (data_labeled$ANOVA[i] <= 0.001) data_labeled$Sig[i]="***"
    else if (data_labeled$ANOVA[i] <= 0.01) data_labeled$Sig[i]="**"
    else if (data_labeled$ANOVA[i] <= 0.05) data_labeled$Sig[i]="*"
    else data_labeled$Sig[i]=""
  }

  write.csv(data_labeled, file=paste0(Title,'-labeled data.csv'), row.names=F)
  return(data_labeled)
}
