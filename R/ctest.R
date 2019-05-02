ctest <-  function(IP,INPUT,TOTAL_IP,TOTAL_INPUT,FOLD=1,minimal_counts_in_fdr = 10) {
  
  # input check
  if (length(IP) != length(INPUT)) { stop("The IP and INPUT of ctest must be of the same length.", call. = TRUE, domain = NULL) }

  # replace 0 with 1
  IP=pmax(IP,1)
  INPUT=pmax(INPUT,1)
  
  # calculate p
  a=TOTAL_IP*FOLD
  b=TOTAL_INPUT
  p=a/(a+b)
  
  # get total observation
  total=IP+INPUT
  
  # cdf
  log.p=pbinom(IP-1, total, p, lower.tail = FALSE, log.p = TRUE)
  
  # calculate p  
  pvalues=exp(log.p)
  
  # calculate fdr
  log.fdr=log(p.adjust(pvalues, method = "fdr"))
  
  # with significant number of reads only
  ID=which( (IP+INPUT) > minimal_counts_in_fdr)
  log.fdr_sig=log(p.adjust(pvalues[ID], method = "fdr"))
  log.fdr[ID]=log.fdr_sig
  
  # fold enrichment
  log.fc=log((pmax(1,IP)/sum(IP))/(pmax(1,INPUT)/sum(INPUT)))
  
  # output result
  PW=list(log.p=log.p,log.fdr=log.fdr,log.fc=log.fc)
  return(PW)

}