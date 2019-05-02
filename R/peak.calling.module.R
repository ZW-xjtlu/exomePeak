.peak.calling.module <- function(READS_COUNT,SAMPLE_ID,PARAMETERS){

# check points comparison
if (length(SAMPLE_ID$ip)>1) {IP=rowSums(READS_COUNT[,SAMPLE_ID$ip])} else {IP=READS_COUNT[,SAMPLE_ID$ip]}
if (length(SAMPLE_ID$input)>1) {INPUT=rowSums(READS_COUNT[,SAMPLE_ID$input])} else {INPUT=READS_COUNT[,SAMPLE_ID$input]}
FOLD=PARAMETERS$POISSON_MEAN_RATIO
PW = ctest(IP,INPUT,sum(IP),sum(INPUT),FOLD)
if (PARAMETERS$PEAK_CUTOFF_TYPE =="FDR") {ID_stat= (PW$log.fdr < log(PARAMETERS$PEAK_CUTOFF_FDR))}
if (PARAMETERS$PEAK_CUTOFF_TYPE =="PVALUE") {ID_stat= (PW$log.p < log(PARAMETERS$PEAK_CUTOFF_PVALUE))}
ID_fc= (PW$log.fc > log(PARAMETERS$FOLD_ENRICHMENT))
global.sig.check.point.id=ID_stat & ID_fc
global.pw=PW

# get sample consistency
consistent.sid = global.sig.check.point.id
noSample = length(SAMPLE_ID$ip)
for (iSample in 1:length(SAMPLE_ID$ip)) {
  IP=READS_COUNT[,SAMPLE_ID$ip[iSample]]
  PW = ctest(IP,INPUT,sum(IP),sum(INPUT),FOLD)
  ID_stat= (PW$log.fdr < log(PARAMETERS$CONSISTENT_PEAK_CUTOFF_PVALUE))
  ID_fc= (PW$log.fc > log(PARAMETERS$CONSISTENT_PEAK_FOLD_ENRICHMENT))
  consistent.sid=consistent.sid & ID_stat & ID_fc
}

# output peaks
PEAK=list(PW=global.pw, Merged=global.sig.check.point.id, Consistent=consistent.sid)
ID=PEAK$Merged
PEAK$loci2peak_merged=.get.peak.from.loci(READS_COUNT,ID,PARAMETERS)
ID=PEAK$Consistent
PEAK$loci2peak_consistent=.get.peak.from.loci(READS_COUNT,ID,PARAMETERS)

return(PEAK)}