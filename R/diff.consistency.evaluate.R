


# evaluate whether the difference is consistent on all the replicates
.diff.consistency.evaluate <- function(PEAK,READS_COUNT,SAMPLE_ID,PARAMETERS){
  
  # total 7 tests
  no_untreated_ip=length(SAMPLE_ID$untreated_ip);
  no_treated_ip=length(SAMPLE_ID$treated_ip);
  no_untreated_input=length(SAMPLE_ID$untreated_input);
  no_treated_input=length(SAMPLE_ID$treated_input);
  
  # consistent
  diff1=list()
  diff2=list()
  diff3=list()
  diff4=list()
    
  # for untreated ip
  for (i in 1:no_untreated_ip) {
    sample_id=SAMPLE_ID
    sample_id$untreated_ip=SAMPLE_ID$untreated_ip[i]
    DIFF=.merge.and.compare.two.different.conditions(PEAK,READS_COUNT,sample_id,PARAMETERS)
    diff1[[i]]=DIFF
  }
  
  # for treated ip
  for (i in 1:no_treated_ip) {
    sample_id=SAMPLE_ID
    sample_id$treated_ip=SAMPLE_ID$treated_ip[i]
    DIFF=.merge.and.compare.two.different.conditions(PEAK,READS_COUNT,sample_id,PARAMETERS)
    diff2[[i]]=DIFF
  } 
  
  # for untreated ip
  for (i in 1:no_untreated_input) {
    sample_id=SAMPLE_ID
    sample_id$untreated_input=SAMPLE_ID$untreated_input[i]
    DIFF=.merge.and.compare.two.different.conditions(PEAK,READS_COUNT,sample_id,PARAMETERS)
    diff3[[i]]=DIFF
  } 
  
  # for treated ip
  for (i in 1:no_treated_input) {
    sample_id=SAMPLE_ID
    sample_id$treated_input=SAMPLE_ID$treated_input[i]
    DIFF=.merge.and.compare.two.different.conditions(PEAK,READS_COUNT,sample_id,PARAMETERS)
    diff4[[i]]=DIFF
  }  
  
  DIFF=c(diff1,diff2,diff3,diff4)
  
  # return result
  return(DIFF)
}