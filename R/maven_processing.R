#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param maven_raw a metabolomics data frame
#'
#' @return a new data frame with data reformatted
#'
# @examples maven_processing(maven_raw)
#'
#' @export
#'
maven_processing<-function(maven_raw){
  maven_raw<-data.frame(maven_raw)

  a<-as.character(maven_raw[1,2:ncol(maven_raw)])
  new_colnames<-c("Name","Iso",a)

  empty_first_vector<-c(rep(NA,12))
  maven_raw<-rbind(empty_first_vector,maven_raw)


  w<-data.frame(matrix(ncol=ncol(maven_raw)+1,nrow=nrow(maven_raw)))
  colnames(w)<-new_colnames

  w_index=1
  for (k in 1:(length(which(is.na(maven_raw$...1)))-1)){
    current_NA<-which(is.na(maven_raw$...1))[k]


    metabolite<-maven_raw[[current_NA+1,1]]
    next_na<-which(is.na(maven_raw$...1))[k+1]
    num_isotopes<-next_na-current_NA-2
    w[w_index:(w_index+num_isotopes-1),1]<-metabolite
    w[w_index:(w_index+num_isotopes-1),2]<-maven_raw[(current_NA+2):(current_NA+1+num_isotopes),1]
    w[w_index:(w_index+num_isotopes-1),3:ncol(w)]<-maven_raw[(current_NA+2):(current_NA+1+num_isotopes),2:ncol(maven_raw)]
    w_index=w_index+num_isotopes
  }

  w<-na.omit(w)

  return(w)
}
