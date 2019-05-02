
.report.peak.based.on.result <- function(PEAK,ANNOTATION,READS_COUNT,SAMPLE_ID,PARAMETERS,
                                         ANNOTATION_BATCH_ID){
  # generate dir
  dir=paste(PARAMETERS$OUTPUT_DIR,PARAMETERS$EXPERIMENT_NAME,sep='/')
  
  # output result 0:merged peak
  print("Get all the peaks ...")
  peak_result=.get.table.peak.result(PEAK,ANNOTATION,READS_COUNT,SAMPLE_ID,PARAMETERS,
                                     ANNOTATION_BATCH_ID,PEAK$loci2peak_merged)
  # new change
  storage.mode(peak_result[[2]]) <- "integer"
  storage.mode(peak_result[[3]]) <- "integer"
  
  merged_peak_result = peak_result
  if (length(peak_result) > 0) {
    TOTAL_PEAK_RESULT=peak_result
    bed_result=peak_result[,1:12];
    names(bed_result)[1]=paste("#",names(bed_result)[1])
    
    if (length(SAMPLE_ID$treated_ip)*length(SAMPLE_ID$treated_input) == 0) {
      # write.result
      if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {
        write.table(peak_result,file=paste(dir,"peak.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)
        write.table(bed_result,file=paste(dir,"peak.bed",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)}
    }
    
  } else {
    stop("No peak detected, try less strigent parameter, e.g., use a larger PEAK_CUTOFF_PVALUE and smaller FOLD_ENRICHMENT")
  }
  
  # output result 1: consistent peak
  print("Get the consistent peaks ...")
  peak_result=.get.table.peak.result(PEAK,ANNOTATION,READS_COUNT,SAMPLE_ID,PARAMETERS,
                                     ANNOTATION_BATCH_ID,PEAK$loci2peak_consistent)
  con_peak_result = peak_result
  if (length(peak_result) > 0) {
    bed_result=peak_result[,1:12];
    dir=paste(PARAMETERS$OUTPUT_DIR,PARAMETERS$EXPERIMENT_NAME,sep='/')
    names(bed_result)[1]=paste("#",names(bed_result)[1])
    
    if (length(SAMPLE_ID$treated_ip)*length(SAMPLE_ID$treated_input) == 0) {
      # write.result
      if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {
      write.table(peak_result,file=paste(dir,"con_peak.xls",sep="/"),sep="\t",row.names =FALSE,quote = FALSE)
      write.table(bed_result,file=paste(dir,"con_peak.bed",sep="/"),sep="\t",row.names =FALSE,quote = FALSE)}
      }} else {
        temp="No peak detected, try less strigent parameter, e.g., use a larger CONSISTENT_PEAK_CUTOFF_PVALUE and smaller CONSISTENT_PEAK_FOLD_ENRICHMENT"
        if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {
        write.table(peak_result,file=paste(dir,"con_peak.xls",sep="/"),sep="\t",row.names =FALSE,quote = FALSE)
        write.table(bed_result,file=paste(dir,"con_peak.bed",sep="/"),sep="\t",row.names =FALSE,quote = FALSE)}
      }
  
  PEAK_RESULT =  list(merged_peak_result,con_peak_result)
  
  # return
  return(PEAK_RESULT)
}