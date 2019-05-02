\name{ctest}
\alias{ctest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ctest
}
\description{
c-test is used to compare two Poisson means, for peak calling or binding sites identification in exomePeak R-package
}
\usage{
ctest(IP, INPUT, TOTAL_IP, TOTAL_INPUT, FOLD = 1, minimal_counts_in_fdr = 10)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{IP}{ a vector of integers, each elemetnt represents the number of reads from a binding site in the IP sample
%%     ~~Describe \code{IP} here~~
}
  \item{INPUT}{ a vector of integers, each elements represents the number of reads from a binding site in the Input control sample
%%     ~~Describe \code{INPUT} here~~
}
  \item{TOTAL_IP}{an integer, which represents the total number of reads in IP sample
%%     ~~Describe \code{TOTAL_IP} here~~
}
  \item{TOTAL_INPUT}{ an integer, which represents the total number of reads in Input control sample
%%     ~~Describe \code{TOTAL_INPUT} here~~
}
  \item{FOLD}{ a decimal number, which indicates the ration of Possion mean to be tested, default: 1. Use a larger number for detection of highly enriched binding sites.
}
  \item{minimal_counts_in_fdr}{ an integer threshold, only the loci with reads more than this number are subjected for fdr calculation. default: 10
}
}
\details{
c-test is used to compare two Poisson means, for peak calling or binding sites identification in exomePeak R-package. The comparison of two
Possion distributions is converted into a binomial distribution based test. The number of reads at the same locus for the aligned reads are
 counted by other packages, such as Rsamtools or HTseq-count.
}
\value{
The function returns a list of length 3, which contains the log(p-value), log(fdr) and log(fold change), respectively.
}

\references{ Przyborowski, J. and Wilenski, H. (1940) Homogeneity of results in testing samples from Poisson series: with an application to testing clover seed for dodder. Biometrika, 31, 313-323
}

\examples{
result = ctest(c(20,10, 1), c(2,1,20), 100, 200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }
