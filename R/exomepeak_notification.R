.exomepeak_notification <- function(PARAMETERS,PEAK_RESULT,DIFF_PEAK_RESULT,dir) {
  
  # save(list =c("PARAMETERS","PEAK_RESULT","DIFF_PEAK_RESULT","dir"), file="test.Rdata")
  
  
  #####################################################
  # notify the input samples
  print("---------------------------------")
  {print("The bam files used:");
  print(paste(length(PARAMETERS$IP_BAM),"IP replicate(s)"))
  print(paste(length(PARAMETERS$INPUT_BAM),"Input replicate(s)"))
  
   respect_mode = "peak"
  if (length(PARAMETERS$TREATED_IP_BAM)*length(PARAMETERS$TREATED_INPUT_BAM)>0) {
  print(paste(length(PARAMETERS$TREATED_IP_BAM),"TREATED IP replicate(s)"))
  print(paste(length(PARAMETERS$TREATED_INPUT_BAM),"TREATED Input replicate(s)"))
  respect_mode = "diff"
  }}
  
  
  #####################################################
  # for peak calling only
  if (respect_mode == "peak" ) {
    print("---------------------------------")
    print("Peak calling result: ")
    print(paste(length(PEAK_RESULT[[1]][,1]),"peaks detected on merged data."))
    if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {print(paste("Please check 'peak.bed/xls' under",dir))}
    print(paste(length(PEAK_RESULT[[2]][,1]),"consistent peaks detected on every replicates. (Recommended list)"))
    if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {print(paste("Please check 'con_peak.bed/xls' under",dir))}
  }
  
  #####################################################
  # peak calling and differential analysis
  if (respect_mode == "diff" ) {
    print("---------------------------------")
    print("Peak calling and differential analysis result: ")
    print(paste(length(DIFF_PEAK_RESULT[[1]][,1]),"peaks detected."))
    
    if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {print(paste("Please check 'diff_peak.bed/xls' under",dir))}
    print("---------------------------------")
    print(paste(length(DIFF_PEAK_RESULT[[2]][,1]),"significantly differential methylated peaks are detected."))
    
    if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {print(paste("Please check 'sig_diff_peak.bed/xls' under",dir))}
    print("---------------------------------")
    print(paste(length(DIFF_PEAK_RESULT[[3]][,1]),"consistent significantly differential methylated peaks are detected.(Recommended list) "))
    if (PARAMETERS$SAVE_RESULT_ON_DISK==TRUE) {print(paste("Please check 'con_sig_diff_peak.bed/xls' under",dir))}
    print("---------------------------------") 
  }
  
  #####################################################
  # convert the peak result
  
  if (respect_mode == "peak" ) {
    
    all_peaks = NA
    con_peaks = NA
    
    if (length(PEAK_RESULT[[1]][,1:12][,1]) > 0) {
    all_peaks = .bed2grangeslist(PEAK_RESULT[[1]][,1:12])
    mcols(all_peaks) <- PEAK_RESULT[[1]][,13:15]}
    
    if (length(PEAK_RESULT[[2]][,1:12][,1]) > 0) {
    con_peaks = .bed2grangeslist(PEAK_RESULT[[2]][,1:12])
    mcols(con_peaks) <- PEAK_RESULT[[2]][,13:15] }  

    RESULT = list(all_peaks=all_peaks, con_peaks=con_peaks)
    
  }
  
  #####################################################
  # convert the differential analysis result 
  
  if (respect_mode == "diff" ) {
    
    diff_peaks = NA
    sig_siff_peaks = NA
    con_sig_diff_peaks = NA
    
    if (length(DIFF_PEAK_RESULT[[1]][,1]) > 0) {
    diff_peaks = .bed2grangeslist(DIFF_PEAK_RESULT[[1]][,1:12])
    mcols(diff_peaks) <- DIFF_PEAK_RESULT[[1]][,13:18]}
    
    if (length(DIFF_PEAK_RESULT[[2]][,1]) > 0) {
    sig_siff_peaks = .bed2grangeslist(DIFF_PEAK_RESULT[[2]][,1:12])
    mcols(sig_siff_peaks) <- DIFF_PEAK_RESULT[[2]][,13:18]}   
    
    if (length(DIFF_PEAK_RESULT[[3]][,1]) > 0) {
    con_sig_diff_peaks = .bed2grangeslist(DIFF_PEAK_RESULT[[3]][,1:12])
    mcols(con_sig_diff_peaks) <- DIFF_PEAK_RESULT[[3]][,13:18]}  
    
    RESULT = list(diff_peaks=diff_peaks, 
                  sig_siff_peaks=sig_siff_peaks,
                  con_sig_diff_peaks=con_sig_diff_peaks)
    
  }  
  
  # result
  return(RESULT)
}