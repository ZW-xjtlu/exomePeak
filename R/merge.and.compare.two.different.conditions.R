
# merge and compare two conditions
.merge.and.compare.two.different.conditions <- function(PEAK,READS_COUNT,SAMPLE_ID,PARAMETERS){
  
  peak_reads_count = .get.peak.reads.count(PEAK,READS_COUNT,SAMPLE_ID,PARAMETERS)
  
  # get reads count
  if (length(SAMPLE_ID$untreated_ip)>1) {untreated_ip=rowSums(peak_reads_count[,SAMPLE_ID$untreated_ip])} else {untreated_ip=(peak_reads_count[,SAMPLE_ID$untreated_ip])}
  if (length(SAMPLE_ID$untreated_input)>1) {untreated_input=rowSums(peak_reads_count[,SAMPLE_ID$untreated_input])} else {untreated_input=(peak_reads_count[,SAMPLE_ID$untreated_input])}
  if (length(SAMPLE_ID$treated_ip)>1) {treated_ip=rowSums(peak_reads_count[,SAMPLE_ID$treated_ip])} else {treated_ip=(peak_reads_count[,SAMPLE_ID$treated_ip])}
  if (length(SAMPLE_ID$treated_input)>1) {treated_input=rowSums(peak_reads_count[,SAMPLE_ID$treated_input])} else {treated_input=(peak_reads_count[,SAMPLE_ID$treated_input])} 
  
  # get total
  sample_total=colSums(READS_COUNT)
  untreated_ip_total=sum(sample_total[SAMPLE_ID$untreated_ip])
  untreated_input_total=sum(sample_total[SAMPLE_ID$untreated_input])
  treated_ip_total=sum(sample_total[SAMPLE_ID$treated_ip])
  treated_input_total=sum(sample_total[SAMPLE_ID$treated_input])
  
  # normalize based on the overlapping window
  sample_total = round (sample_total * PARAMETERS$SLIDING_STEP / PARAMETERS$WINDOW_WIDTH);
  untreated_ip_total = round (untreated_ip_total * PARAMETERS$SLIDING_STEP / PARAMETERS$WINDOW_WIDTH);
  untreated_input_total = round (untreated_input_total * PARAMETERS$SLIDING_STEP / PARAMETERS$WINDOW_WIDTH);
  treated_ip_total = round (treated_ip_total * PARAMETERS$SLIDING_STEP / PARAMETERS$WINDOW_WIDTH);
  treated_input_total = round (treated_input_total * PARAMETERS$SLIDING_STEP / PARAMETERS$WINDOW_WIDTH);
  
  # do test
  if (PARAMETERS$DIFF_PEAK_METHOD == "bltest")  {
  DIFF=bltest(untreated_ip,untreated_input,treated_ip,treated_input,
                                untreated_ip_total,untreated_input_total,treated_ip_total,treated_input_total )}
  if (PARAMETERS$DIFF_PEAK_METHOD == "rhtest")  {
    DIFF=rhtest(untreated_ip,untreated_input,treated_ip,treated_input,
                untreated_ip_total,untreated_input_total,treated_ip_total,treated_input_total )}
  
  # result
  result =list()
  result$DIFF = DIFF
  result$peak_reads_count = list(untreated_ip=untreated_ip,
                                 untreated_input=untreated_input,
                                 treated_ip=treated_ip,
                                 treated_input=treated_input,
                                 untreated_ip_total=untreated_ip_total,
                                 untreated_input_total=untreated_input_total,
                                 treated_ip_total=treated_ip_total,
                                 treated_input_total=treated_input_total,
                                 sample_peak_reads_count=peak_reads_count,
                                 sample_total_reads_count = sample_total,
                                 sample_id=SAMPLE_ID)
  return(result)
  
}