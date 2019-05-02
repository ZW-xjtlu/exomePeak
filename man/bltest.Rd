\name{bltest}
\alias{bltest}
\title{
bltest
}
\description{
This is the default test for the differential post-transcriptional RNA modification sites. Differential from all existing tests the compare the absolute amount between two conditions, this test compares whether the percentage of modified molecules are the same.
}
\usage{
bltest(untreated_ip, untreated_input, 
       treated_ip, treated_input, 
       untreated_ip_total, untreated_input_total, 
       treated_ip_total, treated_input_total,
       minimal_count_fdr =10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{untreated_ip}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the IP sample under untreated condition
%%     ~~Describe \code{untreated_ip} here~~
}
  \item{untreated_input}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the Input control sample under untreated condition
%%     ~~Describe \code{untreated_input} here~~
}
  \item{treated_ip}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the IP sample under treated condition
%%     ~~Describe \code{treated_ip} here~~
}
  \item{treated_input}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the Input control sample under treated condition
%%     ~~Describe \code{treated_input} here~~
}
  \item{untreated_ip_total}{ an integer, total number of reads for the IP sample under untreated condition
%%     ~~Describe \code{untreated_ip_total} here~~
}
  \item{untreated_input_total}{ an integer, total number of reads for the Input control sample under untreated condition
%%     ~~Describe \code{untreated_input_total} here~~
}
  \item{treated_ip_total}{ an integer, total number of reads for the IP sample under treated condition
%%     ~~Describe \code{treated_ip_total} here~~
}
  \item{treated_input_total}{  an integer, total number of reads for the Input control sample under treated condition
%%     ~~Describe \code{treated_input_total} here~~
}
  \item{minimal_count_fdr}{ an integer threshold, only the loci with reads more than this number are subjected for fdr calculation. default: 10
}
}
\value{
The function returns a list of length 3, which contains the log(p-value), log(fdr) and log(fold change), respectively, from the test.
}

\details{
The comparison of 4 Poisson distributions are firstly collapsed into 2 Binomial distributions, and the function further tests whether the two binomial distributions have the same successful rate
 with a likelihood ratio test. The number of reads at the same locus for the aligned reads are
 counted by other packages, such as Rsamtools or HTseq-count.

}

\references{
Reference coming soon!
}
\author{ Lin Zhang, PhD  <laurenie.zhang@gmail.com>
%%  ~~who you are~~
}

\examples{
# input reads count of 3 binding sites
untreated_ip = c(10,20,30)
untreated_input = c(20,20,20)
treated_ip = c(30,10,20)
treated_input  = c(20,20,20)
# sequencing depths
untreated_ip_total = 10^7
untreated_input_total = 10^7
treated_ip_total = 10^7
treated_input_total = 10^7
# get the result
result = bltest(untreated_ip, untreated_input, 
	treated_ip, treated_input, 
	untreated_ip_total, untreated_input_total, 
	treated_ip_total, treated_input_total)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tatistical Inference }
