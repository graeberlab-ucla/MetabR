#' This function produces a matrix for hierarchical clustering or PCA
#'
#' @author Daniel Braas
#' @param MS_data The data frame to be used for the matrix. Can be relative amounts, fractional contribution or percent labeled.
#' @param Type Which type of data is used. Right now this is not a used variable, but will become important if I add something for isotopologues. 
#' @param anova A cutoff value: usually a p-value of some sort
#'
#' @return
#' @export
#'
#' @examples
make_matrix <- function(MS_data, Type, anova=0.05){
  if (sum(grepl('MID', names(MS_data))) >= 1) {
    data9 <- MS_data %>%
      filter(ANOVA <= anova) %>%
      select(Name, Condition, Iso, contains('MID')) %>%
      gather(Exp, Value, -Name, -Condition, -Iso) %>%
      arrange(Name, Condition) %>%
      #mutate(Condition_Exp = Condition:Exp, Condition=NULL, Exp=NULL) %>%
      #mutate(Name_Iso = Name:Iso, Name=NULL, Iso=NULL) %>%
      unite(Condition_Exp, c(Condition, Exp), sep='_') %>%
      mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
      unite(Name_Iso, c(Name, Iso), sep='_') %>%
      spread(Condition_Exp, Value)
    #data9 <- data9[,!(colSums(data9[2:length(data9)], na.rm=T)==0)]
    data9 <- data9[!(rowSums(data9[2:length(data9)], na.rm=T)==0),]
    
  } else {
    
    data9 <- MS_data %>%
      filter(ANOVA <= anova) %>%
      select(Name, Condition, grep('Exp|FC|Labeled', names(MS_data))) %>%
      gather(Exp, Value, -Name, -Condition) %>%
      arrange(Name, Condition) %>% 
      #mutate(Exp = as.numeric(gsub('Exp|FC|Labeled','',.$Exp))) %>%
      group_by(Name) %>%
      mutate(Nr.NA = sum(is.na(Value)),
             Nr.Samples = n()) %>%
      filter(Nr.NA < Nr.Samples - 1) %>%
      ungroup() %>%
      select(-Nr.NA, -Nr.Samples) %>%
      #arrange(Name, Condition, Exp) %>%
      unite(Condition_Exp, c(Condition, Exp), sep='_') %>%
      mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
      spread(Condition_Exp, Value)
  }
  
  data9[is.na(data9)] <- 0
  data5=as.matrix(data9[2:length(data9)])
  rownames(data5) <- data9$Name
  data5 <- data5[!(rowSums(data5))==0,]
  data5 <- data5[,!(colSums(data5))==0]
  
  return(data5)
}