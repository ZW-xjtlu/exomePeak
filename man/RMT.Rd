\name{RMT}
\alias{RMT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract a combinatorial RNA methylome from multiple biological conditions or replicates.
}
\description{
This function will extract a combinatorial RNA methylome (tens of thousands of RNA methylation sites) from multiple biological conditions or replicates and count with MeRIP-Seq datasets.
}
\usage{
RMT(INPUT_BAM, 
    IP_BAM, 
    INPUT2IP = NA, 
    GENE_ANNO_GTF = NA, 
    GENOME = NA, 
    UCSC_TABLE_NAME = "knownGene",
    TXDB = NA,
    EXOME_OUTPUT_DIR = NA, 
    RMT_OUTPUT_DIR = NA, 
    RMT_EXPERIMENT_NAME = "RMT_Result")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{INPUT_BAM}{ a vector of strings, which specifies a number of Input control samples from the untreated condition in bam format
%%     ~~Describe \code{INPUT_BAM} here~~
}
  \item{IP_BAM}{ a vector of strings, which specifies a number of IP samples from the untreated condition in bam format
%%     ~~Describe \code{IP_BAM} here~~
}
  \item{INPUT2IP}{ a list of number, the default is NA, if it is NA, the length between input and ip is same, and each input bam corresponds one ip bam; otherwise, the length of INPUT2IP equals the number of IP_BAM, and each elment specifies the indexes of the input bams that are corresponding to a specific ip bam.
%%     ~~Describe \code{INPUT2IP} here~~ 
}
  \item{GENOME}{ a string,such as "hg19" or "mm9", which specifies the genome assembly used. If a gene annotation file is provided, the exomepeak will call peaks with it; otherwise, exomepeak will 
     	download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
%%     ~~Describe \code{GENE_ANNO_GTF} here~~
}
  \item{UCSC_TABLE_NAME}{ a string, which specifies the gene annotation used from UCSC, default: "knownGene". Please use function: supportedUCSCtables() to check available tables.
	Some tables may not be available for all genomes, and the "refGene" table doesn't work correctly due to multiple occuences of the same transcript on the same chromosome.
%%     ~~Describe \code{GENE_ANNO_GTF} here~~
}
  \item{GENE_ANNO_GTF}{ a string, which specifies a gene annotation GTF file if available, default: NA
%%     ~~Describe \code{GENE_ANNO_GTF} here~~
}
  \item{TXDB}{ an optional TxDb object for gene annotation information used in the analysis, default: NA. The exomepeak function will 
	first look at TXDB, then GENE_ANNO_GTF, and then GENOME for gene annnotation information. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
%%     ~~Describe \code{GENE_ANNO_GTF} here~~
}
  
  \item{EXOME_OUTPUT_DIR}{ a string, which specifies folder name generated in the output directory that contains all the results, default: "exomePeak_output"
%%     ~~Describe \code{EXOME_OUTPUT_DIR} here~~
}
  \item{RMT_OUTPUT_DIR}{ a string, which specifies the GO output directory, default: 
getwd()
%%     ~~Describe \code{GO_OUTPUT_DIR} here~~
}
  \item{RMT_EXPERIMENT_NAME}{ a string, which specifies folder name generated in the GO output directory that contains all the results, default: "GO_result"
%%     ~~Describe \code{GO_EXPERIMENT_NAME} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Results will be saved two folders (default: "respect_output" and "GO_result") in the specified directory (default: current working directory).
1.In the exomPeak folder, store the information of all peak in bed and xls format
2.In the GO_result folder, store the information of reads counts, rpkm in ip bams and input bams, fold change with ip and input in xls format
}
\references{
Liu, Lian, Shao-Wu Zhang, Yu-Chen Zhang, Hui Liu, Lin Zhang, Runsheng Chen, Yufei Huang, and Jia Meng. "Decomposition of RNA methylome reveals co-methylation patterns induced by latent enzymatic regulators of the epitranscriptome." Molecular BioSystems 11, no. 1 (2015): 262-274.
}
\author{
Lian Liu <liulian19860905@163.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# the RESPECT R-package has two main functions:
# 1. exome-based peak
# 2. compute reads counts, rpkm and fold change
# please feel free to contact liulian19860905@163.com for any questions

# load library and specify the parameters
gene_anno_gtf <- system.file("extdata", "example.gtf", package="exomePeak")
f1 <- system.file("extdata", "IP1.bam", package="exomePeak")
f2 <- system.file("extdata", "IP2.bam", package="exomePeak")
f3 <- system.file("extdata", "IP3.bam", package="exomePeak")
f4 <- system.file("extdata", "IP4.bam", package="exomePeak")
ip_bam <- c(f1,f2,f3,f4)
f1 <- system.file("extdata", "Input1.bam", package="exomePeak")
f2 <- system.file("extdata", "Input2.bam", package="exomePeak")
f3 <- system.file("extdata", "Input3.bam", package="exomePeak")
input_bam <- c(f1,f2,f3)
input2bam <- list(4)
input2bam[[1]] <- c(1,2) # the 1st ip sample uses the 1st & 2nd input samples as control
input2bam[[2]] <- c(1,2) # the 2nd ip sample uses the 1st & 2nd input samples as control
input2bam[[3]] <- c(1,2) # the 3rd ip sample uses the 1st & 2nd input samples as control
input2bam[[4]] <- c(3)   # the 4th ip sample uses the 3rd input sample as control

# Extract the combinatorial RNA methylome
# unfortunately, this function has not been optimized for parallel processing.
# This part will take a really long time on real data
RMT(INPUT_BAM=input_bam, IP_BAM=ip_bam, INPUT2IP=input2bam, GENE_ANNO_GTF=gene_anno_gtf)

# RMT result will be generated in RMT_result folder under current working directoy, including:
# the merged peaks or the combinatorial RNA methylome from all replicates
# the bed file, which can be viewed under IGV
# reads count of all peaks for both the ip and input samples
# RPKM of all the peaks for both the ip and input samples (for your reference only).

# Alternatively, exomePeak package can automatically download the complete transcriptome.
# And then scan the entire transcriptome for RNA methylation sites
# It will take a long time
# RMT(INPUT_BAM=input_bam, IP_BAM=ip_bam, INPUT2IP=input2bam, GENOME="hg19")


}

% R documentation directory.
\keyword{MeRIP-Seq}
