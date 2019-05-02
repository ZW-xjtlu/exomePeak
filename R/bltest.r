
# differential analysis with LRT
bltest <- function(untreated_ip,untreated_input,treated_ip,treated_input,
                   untreated_ip_total,untreated_input_total,treated_ip_total,treated_input_total,
                   minimal_count_fdr = 10) {
  # parameters
  
  # parameters
  # untreated_ip_total: an integer, sequencing depth of untreated ip sample   --m0
  # untreated_input_total: ... untreated input sample                         --n0
  # treated_ip_total: ... treated ip sample                                   --m1
  # treated_input_total: ... treated input sample                             --n1
  
  # untreated_ip: y0 a vector of integers, vector length equal to the number of binding sites, each element contains the number of reads within a binding site for the untreated ip sample
  # untreated_input: x0 ... of the untreated input sample
  # treated_ip: y1 ... of the treated ip sample
  # treated_input: x1 ... of the treated input sample
  
    # updated 2015-8-8  
    untreated_ip_temp <- untreated_ip
    treated_input_temp <- treated_input
    untreated_input_temp <- untreated_input
    treated_ip_temp <- treated_ip
    if ((untreated_ip_total*treated_input_total) > (untreated_input_total*treated_ip_total)){
    # down-size untreated_ip_total*treated_input_total
    if (untreated_ip_total > treated_input_total) {
      temp=(untreated_input_total*treated_ip_total)/treated_input_total;
      untreated_ip_temp=round(untreated_ip_temp*temp/untreated_ip_total);
    } else {
      temp=(untreated_input_total*treated_ip_total)/untreated_ip_total
      treated_input_temp=round(treated_input_temp*temp/treated_input_total);
    }
  }
  
  if ((untreated_ip_total*treated_input_total) < (untreated_input_total*treated_ip_total)) {
    if (untreated_input_total>treated_ip_total) {
      temp=round(untreated_ip_total*treated_input_total/treated_ip_total)
      untreated_input_temp=round(untreated_input_temp*temp/untreated_input_total)   
    } else {
      temp=(untreated_ip_total*treated_input_total)/untreated_input_total
      treated_ip_temp=round(treated_ip_temp*temp/treated_ip_total)
    }
  }
  
    # replace 0 with 1 to avoid ridiculous p-values
  untreated_ip_temp=pmax(1,untreated_ip_temp)
  untreated_input_temp=pmax(1,untreated_input_temp)
  treated_ip_temp=pmax(1,treated_ip_temp)
  treated_input_temp=pmax(1,treated_input_temp)
  
  untreated_input[untreated_input==0] <- 1
  treated_input[treated_input==0] <- 1
  untreated_ip[untreated_ip==0] <- 1
  treated_ip[treated_ip==0] <- 1
  
  t0 <- untreated_input + untreated_ip #t0
  t1 <- treated_input + treated_ip
  
  # Likelihood-Ratio Test Statistic
  # Null hypothesis
  a <-  (untreated_input + treated_input) *untreated_ip_total *treated_ip_total
  #a1 <- a
  #a[a==0] <- 1;
  b1 <- (untreated_input-treated_ip)*untreated_ip_total*treated_input_total
  b2 <- (treated_input-untreated_ip)*treated_ip_total*untreated_input_total
  b <- b1 + b2
  c <- -(untreated_ip + treated_ip) * untreated_input_total * treated_input_total
  #a1[a==0] <- 1;
  p.hat.ml <- (-b + sqrt(b^2-4*a*c))/(2*a)
  
  p0.bar <- untreated_input_total/(untreated_input_total+p.hat.ml*untreated_ip_total)
  p1.bar <- treated_input_total/(treated_input_total+p.hat.ml*treated_ip_total)
  
  # Alternative hypothesis
  p0.hat <- untreated_input/t0
  p1.hat <- treated_input/t1
  
  # Original
  # log.fc=log((treated_ip/treated_input)/(untreated_ip/untreated_input))
  # updated 2014-9-10  
  # log.fc=log(((treated_ip/treated_ip_total)/(treated_input/treated_input_total))/((untreated_ip/untreated_ip_total)/(untreated_input/untreated_input_total)))

  log.fc=log((treated_ip_temp/treated_input_temp)/(untreated_ip_temp/untreated_input_temp))
  
  D <- 2 * (untreated_input*log(p0.hat/p0.bar)+untreated_ip*log((1-p0.hat)/(1-p0.bar))+treated_input*log(p1.hat/p1.bar)+treated_ip*log((1-p1.hat)/(1-p1.bar)))
#  DIFF[a==0] <- min(DIFF)
  # test the difference
  # test the difference
  
  # calculate p  
  #pvalues=exp(log.p)
  pvalues = 1-pchisq(D, 1)
  log.p = log(pvalues)
    
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
  
  # save result
  #DIFF=list(D)
  #return(DIFF)
}