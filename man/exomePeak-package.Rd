\name{exomePeak-package}
\alias{exomePeak-package}
\alias{exomePeak}
\docType{package}
\title{
exomePeak
}
\description{
The package is developed for the analysis of affinity-based epitranscriptome shortgun sequencing data from MeRIP-seq (maA-seq). It was built on the basis of the exomePeak MATLAB package (Meng, Jia, et al. "Exome-based analysis for RNA epigenome sequencing data." Bioinformatics 29.12 (2013): 1565-1567.) with new functions for differential analysis of two experimental conditions to unveil the dynamics in post-transcriptional regulation of the RNA methylome. The exomePeak R-package accepts and statistically supports multiple biological replicates, internally removes PCR artifacts and multi-mapping reads, outputs exome-based binding sites (RNA methylation sites) and detects differential post-transcriptional RNA modification sites between two experimental conditions in term of percentage rather the absolute amount. The package is still under active development, and we welcome all biology and computation scientist for all kinds of collaborations and communications. Please feel free to contact Dr. Jia Meng <jia.meng@hotmail.com> if you have any questions.
}
\details{
\tabular{ll}{
Package: \tab exomePeak\cr
Type: \tab Package\cr
Version: \tab 1.9.1\cr
Date: \tab 2015-06-25\cr
License: \tab GPL-2 \cr
}

}

\references{

For peak calling or general purpose:
Meng, Jia, Xiaodong Cui, Manjeet K. Rao, Yidong Chen, and Yufei Huang. "Exome-based analysis for RNA epigenome sequencing data." Bioinformatics 29, no. 12 (2013): 1565-1567.

For differential analysis with rhtest:
Meng, Jia, Zhiliang Lu, Hui Liu, Lin Zhang, Shaowu Zhang, Yidong Chen, Manjeet K. Rao, and Yufei Huang. "A protocol for RNA methylation differential analysis with MeRIP-Seq data and exomePeak R/Bioconductor package." Methods 69, no. 3 (2014): 274-281.

For combinatorial RNA methylome with RMT:
Liu, Lian, Shao-Wu Zhang, Yu-Chen Zhang, Hui Liu, Lin Zhang, Runsheng Chen, Yufei Huang, and Jia Meng. "Decomposition of RNA methylome reveals co-methylation patterns induced by latent enzymatic regulators of the epitranscriptome." Molecular BioSystems 11, no. 1 (2015): 262-274.

}

\keyword{ Statistical Inference }

\examples{
# For usage, please check the main function with: 
?exomepeak
}
