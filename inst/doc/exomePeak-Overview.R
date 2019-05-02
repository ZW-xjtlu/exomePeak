### R code from vignette source 'exomePeak-Overview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: input bams files
###################################################
library("exomePeak")
gtf <- system.file("extdata", "example.gtf", package="exomePeak")
f1 <- system.file("extdata", "IP1.bam", package="exomePeak")
f2 <- system.file("extdata", "IP2.bam", package="exomePeak")
f3 <- system.file("extdata", "IP3.bam", package="exomePeak")
f4 <- system.file("extdata", "IP4.bam", package="exomePeak")
f5 <- system.file("extdata", "Input1.bam", package="exomePeak")
f6 <- system.file("extdata", "Input2.bam", package="exomePeak")
f7 <- system.file("extdata", "Input3.bam", package="exomePeak")


###################################################
### code chunk number 3: Peak Calling
###################################################
result <- exomepeak(GENE_ANNO_GTF=gtf, 
                   IP_BAM=c(f1,f2,f3,f4), 
                   INPUT_BAM=c(f5,f6,f7))
names(result)


###################################################
### code chunk number 4: extract peaks
###################################################
recommended_peaks <- result$con_peaks # consistent peaks (Recommended!)
peaks_info <- mcols(recommended_peaks) # information of the consistent peaks
head(peaks_info)


###################################################
### code chunk number 5: extract peaks 2
###################################################
all_peaks <- result$all_peaks # get all peaks
peaks_info <- mcols(all_peaks) # information of all peaks
head(peaks_info)


###################################################
### code chunk number 6: input bams files 2
###################################################
library("exomePeak")
gtf <- system.file("extdata", "example.gtf", package="exomePeak")
f1 <- system.file("extdata", "IP1.bam", package="exomePeak")
f2 <- system.file("extdata", "IP2.bam", package="exomePeak")
f3 <- system.file("extdata", "IP3.bam", package="exomePeak")
f4 <- system.file("extdata", "IP4.bam", package="exomePeak")
f5 <- system.file("extdata", "Input1.bam", package="exomePeak")
f6 <- system.file("extdata", "Input2.bam", package="exomePeak")
f7 <- system.file("extdata", "Input3.bam", package="exomePeak")
f8 <- system.file("extdata", "treated_IP1.bam", package="exomePeak")
f9 <- system.file("extdata", "treated_Input1.bam", package="exomePeak")


###################################################
### code chunk number 7: Peak Calling and Differential Analysis
###################################################
result <- exomepeak(GENE_ANNO_GTF=gtf, 
                   IP_BAM=c(f1,f2,f3,f4), 
                   INPUT_BAM=c(f5,f6,f7),
                   TREATED_IP_BAM=c(f8), 
                   TREATED_INPUT_BAM=c(f9))


###################################################
### code chunk number 8: ScanBamParam
###################################################
names(result)
is.na(result$con_sig_diff_peaks) # no reported consistent differnetial peaks


###################################################
### code chunk number 9: ScanBamParam 2
###################################################
diff_peaks <- result$diff_peaks # consistent differential peaks (Recommended!)
peaks_info <- mcols(diff_peaks) # information of the consistent peaks
head(peaks_info[,1:3]) # peak calling information
head(peaks_info[,4:6]) # differential analysis information


###################################################
### code chunk number 10: Peak Calling and Differential Analysis (eval = FALSE)
###################################################
## result <- exomepeak(GENOME="hg19", 
##                    IP_BAM=c(f1,f2,f3,f4), 
##                    INPUT_BAM=c(f5,f6,f7),
##                    TREATED_IP_BAM=c(f8), 
##                    TREATED_INPUT_BAM=c(f9))


###################################################
### code chunk number 11: session
###################################################
sessionInfo()


