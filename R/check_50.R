check_50<-function(df,type){
 
  df<-data.frame(df)
  

  text=""
  if (type=="FC"){
    text="FC"
  }
  else if (type=="data4"){
    text="Exp"
  }
  else if (type=="MID"){
    text="MID"
  }
  
num_exp<-length(grep(text,colnames(df)))
thresh<-round(length(grep(text,colnames(df)))/2)
condition_samples<-grep(text,colnames(df))

df[which(is.na(df$Exp1)),condition_samples[1]]<-0
df[which(is.na(df$Exp2)),condition_samples[2]]<-0
under_50_percent<-c()
for (n in 1:nrow(df)){
  value_exist=0
  for (x in 1:num_exp){
    
    if ((df[n,condition_samples[x]])!= 0){
      value_exist=value_exist+1
    }}
  if (value_exist<thresh){
    under_50_percent<-c(under_50_percent,"X")
  } 
  else{under_50_percent<-c(under_50_percent,"")}
  
}

new<-data.frame(cbind(df,under_50_percent))
return (new)}
