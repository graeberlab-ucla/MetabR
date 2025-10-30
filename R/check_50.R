#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param df a metabolomics data frame
#'
#' @return a new data frame with indication of the cases in which less than 50% of samples have the metabolite measured
#'
# @examples check_50(data4)
#'
#' @export
#'
check_50<-function(df)
{

  df <- data.frame(df)

  text <- ""
  if(sum(grepl('MID', names(df))) >= 1)
  {
    text <- "MID"
  } else if(sum(grepl('FC', names(df))) >= 1)
  {
    text <- "FC"
  }else if(sum(grepl('Exp', names(df))) >= 1)
  {
    text <- "Exp"
  }else if(sum(grepl('Labeled', names(df))) >= 1)
  {
    text  <- "Labeled"
  }


  num_exp<-length(grep(text,colnames(df)))
  thresh<-round(length(grep(text,colnames(df)))/2)
  condition_samples<-grep(text,colnames(df))

  for (i in condition_samples)
  {
      df[which(is.na(df[,i])), i] <- 0
  }

  under_50_percent<-c()
  for (n in 1:nrow(df)){
    value_exist=0
    for (x in 1:num_exp){
      if ((df[n,condition_samples[x]])!= 0){
        value_exist=value_exist+1
      }}
    if (value_exist<thresh & value_exist != 0 & !grepl("blank|250K|QC-0|50K|PB", df$Condition[n], ignore.case = T)){
      under_50_percent<-c(under_50_percent,"X")
    }
    else{under_50_percent<-c(under_50_percent,"")}

  }

  new <- data.frame(cbind(df,under_50_percent))
  return (new)
}

