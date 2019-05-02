
# differential analysis
rhtest <- function(untreated_ip,untreated_input,treated_ip,treated_input,
                                     untreated_ip_total,untreated_input_total,treated_ip_total,treated_input_total,
                                         minimal_count_fdr = 10) {
  # parameters
  # untreated_ip_total: an integer, sequencing depth of untreated ip sample 
  # untreated_input_total: ... untreated input sample
  # treated_ip_total: ... treated ip sample
  # treated_input_total: ... treated input sample
  
  # untreated_ip: a vector of integers, vector length equal to the number of binding sites, each element contains the number of reads within a binding site for the untreated ip sample
  # untreated_input: ... of the untreated input sample
  # treated_ip: ... of the treated ip sample
  # treated_input: ... of the treated input sample
  
  # rescaled
  if ((untreated_ip_total*treated_input_total) > (untreated_input_total*treated_ip_total)){
    # down-size untreated_ip_total*treated_input_total
    if (untreated_ip_total > treated_input_total) {
      temp=(untreated_input_total*treated_ip_total)/treated_input_total;
      untreated_ip=round(untreated_ip*temp/untreated_ip_total);
    } else {
      temp=(untreated_input_total*treated_ip_total)/untreated_ip_total
      treated_input=round(treated_input*temp/treated_input_total);
    }
  }
  
  if ((untreated_ip_total*treated_input_total) < (untreated_input_total*treated_ip_total)) {
    if (untreated_input_total>treated_ip_total) {
      temp=round(untreated_ip_total*treated_input_total/treated_ip_total)
      untreated_input=round(untreated_input*temp/untreated_input_total)   
    } else {
      temp=(untreated_ip_total*treated_input_total)/untreated_input_total
      treated_ip=round(treated_ip*temp/treated_ip_total)
    }
  }
  
  # replace 0 with 1 to avoid ridiculous p-values
  untreated_ip=pmax(1,untreated_ip)
  untreated_input=pmax(1,untreated_input)
  treated_ip=pmax(1,treated_ip)
  treated_input=pmax(1,treated_input)
  
  # test the difference
  q=treated_ip
  m=treated_ip+untreated_ip
  n=untreated_input+treated_input
  k=treated_ip+treated_input
  log.hyper.p = phyper(q-1, m, n, k, lower.tail = FALSE, log.p = TRUE)
  log.hypo.p = phyper(q, m, n, k, lower.tail = TRUE, log.p = TRUE)
  log.fc=log((treated_ip/treated_input)/(untreated_ip/untreated_input))
  log.p=pmin(log.hyper.p,log.hypo.p)
  
  # calculate p  
  pvalues=exp(log.p)
  
  # calculate fdr
  log.fdr=log(p.adjust(pvalues, method = "fdr"))
  # adjust only sig reads count testing result
  m = untreated_ip + untreated_input + treated_ip + treated_input
  ID= which(m > minimal_count_fdr)
  log.fdr_sig=log(p.adjust(pvalues[ID], method = "fdr"))
  log.fdr[ID] = log.fdr_sig
  
  # remove infinity
  log.fdr=pmax(log.fdr,-1000)
  log.p=pmax(log.p,-1000)
  
  # save result
  DIFF=list(log.fdr=log.fdr,log.p=log.p,log.fc=log.fc)
  return(DIFF)
}