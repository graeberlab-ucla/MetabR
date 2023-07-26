#' assign_sample_colors
#'
#' @param info info sheet Condition and Samples columns required
#'
#' @return info sheet with colors column that are assigned to corresponding conditions
#' @export
#'
#' @examples info<-assign_sample_colors(info)
#'
assign_sample_colors<-function(info){
  color_lst<-c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","darkgreen","lightpink1","navy","olivedrab1",
               "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
               "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")

  pool_color_lst<-c("darkorange","darkorange2", "darkorange3", "darkorange4")
  extra_qc <- c("peachpuff1", "seashell1", "wheat2", "snow1")

  info$Condition <- factor(info$Condition, levels = unique(info$Condition))
  Freq <- data.frame(table(info$Condition))
  Freq$color<-NA

  #samples
  samples<-Freq[which(!grepl("QC|blank|250|pool", Freq$Var1)), ]
  Freq[which(Freq$Var1 %in% samples$Var1),'color']<-color_lst[1:nrow(samples)]

  #pools
  Freq[which(grepl("[Pp]ool", Freq$Var1)), 'color' ]<-"darkorange"

  #uncomment if orange color range across pools is needed
  #pools<-Freq[which(grepl("[Pp]ool", Freq$Var1)), ]
  #Freq[which(Freq$Var1 %in% pools$Var1),'color']<-pool_color_lst[1:nrow(pools)]

  #blank and 250ks
  Freq[which(grepl("QC[-.]blank", Freq$Var1)), 'color' ]<-"grey45"
  Freq[which(grepl("QC[-.]250[Kk]|50[Kk]", Freq$Var1)),'color' ]<-"yellow1"

  #extra blanks
  Freq[(which(is.na(Freq$color) & grepl('blank', Freq$Var1, ignore.case = T))),'color']<-"grey20"

  #make it work for exceptions:
  Freq[which(is.na(Freq$color)),'color']<-extra_qc[1:nrow(Freq[which(is.na(Freq$color)),])]

  m_ind<-match(info$Condition, Freq$Var1)
  info$color<-Freq$color[m_ind]
  return(info)
}
