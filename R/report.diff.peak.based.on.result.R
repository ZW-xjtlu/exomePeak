

# comparing two experimental conditions
.report.diff.peak.based.on.result <- function(TOTAL_PEAK_RESULT,DIFF,DIFF.consis,PARAMETERS) {
  
  # get dir
  dir=paste(PARAMETERS$OUTPUT_DIR,PARAMETERS$EXPERIMENT_NAME,sep='/')
  
  # get sig digits
  lg.p=signif(DIFF$DIFF$log.p/log(10), digits = 3)
  lg.fdr=signif(DIFF$DIFF$log.fdr/log(10), digits = 3)
  log2.fc=signif(DIFF$DIFF$log.fc/log(2), digits = 3)
  
  # get the differential information
  diff_score=cbind(lg.fdr, lg.p, log2.fc)
  colnames(diff_score) = c("diff.lg.fdr","diff.lg.p","diff.log2.fc");
  diff_report_total = cbind(TOTAL_PEAK_RESULT,diff_score)
  if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {write.table(diff_report_total,file=paste(dir,"diff_peak.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)}
  all_diff_peak = diff_report_total
  temp=diff_report_total[,1:12]
  names(temp)[1]=paste("#",names(temp)[1])
  if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {write.table(temp,file=paste(dir,"diff_peak.bed",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)}
  
  # get significant peaks
  if (PARAMETERS$DIFF_PEAK_CUTOFF_TYPE=="PVALUE") {
    ID1= (DIFF$DIFF$log.p < log(PARAMETERS$DIFF_PEAK_CUTOFF_PVALUE))
  } else {
    ID1= (DIFF$DIFF$log.fdr < log(PARAMETERS$DIFF_PEAK_CUTOFF_FDR))
  }
  ID2=(abs(DIFF$DIFF$log.fc) > log(PARAMETERS$DIFF_PEAK_ABS_FOLD_CHANGE))
  ID=which(ID1 & ID2)
  if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {write.table(diff_report_total[ID,],file=paste(dir,"sig_diff_peak.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)}
  sig_diff_peak = diff_report_total[ID,]
  
  # get significant bed
  temp=diff_report_total[ID,1:12]
  temp[,5]=signif(exp(diff_report_total[ID,17]),2)
  names(temp)[1]=paste("#",names(temp)[1])
  if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {write.table(temp,file=paste(dir,"sig_diff_peak.bed",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)}

  # consistent
  no_rep=length(DIFF.consis)
  ID3=ID1
  for (irep in 1:no_rep){
    # p-value
    temp1=  (DIFF.consis[[irep]]$DIFF$log.p < log(PARAMETERS$DIFF_PEAK_CONSISTENT_CUTOFF_PVALUE))
    ID1 = ID1 & temp1
    # abs fc
    temp2  = (abs(DIFF.consis[[irep]]$DIFF$log.fc) > log(PARAMETERS$DIFF_PEAK_CONSISTENT_ABS_FOLD_CHANGE))
    ID2 = ID2 & temp2
    # fc sign the same
    temp3= ((DIFF$DIFF$log.fc * DIFF.consis[[irep]]$DIFF$log.fc) > 0)
    ID3 = ID3 & temp3
  }
  ID=which(ID1 & ID2 & ID3)
  if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {write.table(diff_report_total[ID,],file=paste(dir,"con_sig_diff_peak.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)}
  con_diff_peak = diff_report_total[ID,]
    
  # get consistent significant bed
  temp=diff_report_total[ID,1:12]
  temp[,5]=signif(exp(diff_report_total[ID,17]),2)
  names(temp)[1]=paste("#",names(temp)[1])
  if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {write.table(temp,file=paste(dir,"con_sig_diff_peak.bed",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)} 
  
  RESULT=list(all_diff_peak,sig_diff_peak,con_diff_peak)
  return(RESULT)
}