#' info_with_subset_blanks
#'
#' @param infoA data frame that has columns Sample, Run.Order, Condition, etc.
#' @description This function takes a subset sample info sheet and only keeps the two QC-Blanks that were run directly before/ after the samples.
#' The purpose of this function is to prevent unwanted filtering before calling filter.R,
#'  that can arise when we include QC-blank values that are not associated with this subset of samples.
#' @return An info sheet that only contains the two QC-blanks that pertain to this subset.
#' @export
#'
#' @examples
info_with_subset_blanks<-function(info){
#only use on subsets
samples<-info[which(!grepl("QC|[Pp]ool|blank", info$Sample)),] %>%
  arrange(Run.Order)

first_run_sample<-samples[1,]
last_run_sample<-samples[dim(samples)[1],]

qc_blanks<-info[which(grepl("QC[-.][Bb]lank",info$Sample)),] %>%
  arrange(Run.Order)

minqc_diff<-as.data.frame(cbind(qc_blanks$Sample, as.numeric(abs(first_run_sample$Run.Order-qc_blanks$Run.Order))))
minqc<-min(as.numeric(minqc_diff$V2))
blank1<-minqc_diff%>%filter(V2==minqc)
blank1<-blank1$V1 #QC-Blank run right before samples

lastqc_diff<-as.data.frame(cbind(qc_blanks$Sample, as.numeric(abs(last_run_sample$Run.Order-qc_blanks$Run.Order))))
lastqc<-min(as.numeric(lastqc_diff$V2))
blank2<-lastqc_diff%>%filter(V2==as.character(lastqc))
blank2<-blank2$V1 #QC-Blank run right after samples

if(dim(qc_blanks)[1]>2){
qc_to_remove<-qc_blanks[which(!grepl(blank1, qc_blanks$Sample) & !grepl(blank2, qc_blanks$Sample)),]
m_ind<-match(qc_to_remove$Sample, info$Sample)
info<-info[-m_ind,]
}

return(info)
}




