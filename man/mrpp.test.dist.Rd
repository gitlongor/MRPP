\name{mrpp.test}
\Rdversion{1.1}
\alias{mrpp.test.dist}
\alias{mrpp.test.default}
\alias{mrpp.test.formula}
\alias{mrpp.test}
\alias{mrpp}
\alias{p.value}
\alias{p.value.mrpp.test}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
MRPP test for one-way design
}
\description{
K-sample MRPP test
}
\usage{
mrpp.test(y, \dots)
\method{mrpp.test}{dist}(y, trt, B = as.integer(min(nparts(table(trt)), 1e4L)), 
		permutedTrt, weight.trt='df', eps=1e-8, \dots)
\method{mrpp.test}{formula}(y, data, \dots)
\method{mrpp.test}{default}(y, \dots)
\method{p.value}{mrpp.test}(x, type=c("raw","midp"), \dots)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{y}{A N*P data matrix, \code{\link{dist}} object of the data matrix, or a \code{\link{formula}}. 
}
  \item{trt}{A vector of treatment assignments. If missing, the first permutation in \code{permutedTrt} will be used as the original treatment assignment. As least one of \code{trt} or \code{permutedTrt} need to be given. 
}
  \item{B}{A positive integer number of permutations requested. This will only be used when \code{permutedTrt} is missing. This will be overwritten by the \code{ncol(permutedTrt)}.
}
  \item{permutedTrt}{An optional permutation index matrix. If missing, this will be computed as \code{permuteTrt(trt,B)}. See \code{\link{permuteTrt}}. 
}
  \item{weight.trt}{Either a character scalar being one of \code{c("df", "n", "equal", "terms")} or a positive numeric vector of length \code{nlevels(trt)} subject to sum-to-1 constraint. See details. 
}
  \item{eps}{A small non-negative number, differences below which among permuted test statistics are treated as equal.
}
  \item{data}{A data frame in which variables in the formula can be found.}
  \item{\dots}{Additional arguments passed to the methods. 
}
  \item{x}{An object of class \code{"mrpp.test"}.
}
  \item{type}{A character scalar being one of \code{c("raw","midp")}, specifying the type of the p-value to be extracted from \code{x}.
}
}
\details{
For the \code{weight.trt} argument, suppose that \code{n} denote the vector sample sizes for each treatment. Then, 
\itemize{
\item{\code{weight.trt="df"} uses group (treatment) weights proportional to \code{n-1};}
\item{ \code{weight.trt="df"} uses group (treatment) weights proportional to \code{n};}
\item{\code{weight.trt="equal"} uses equal group (treatment) weights; and }
\item{\code{weight.trt="terms"} uses group (treatment) weights proportional to \code{choose(n,2)}.}
} 

}

\value{
\code{mrpp.test} methods return an object of class \code{c("mrpp.test", "htest"}, which is a list with the following components:
\item{statistic }{The observed test statistic.}
\item{all.statistics }{\emph{All} permuted test statistics.}
\item{p.value }{The permutation p-value, with the \code{"midp"} attribute set to the mid-p-value.}
\item{parameter }{A vector of the number of permutations and the weighting method. The latter is coded as an integer, with levels given in the \code{"weight methods"} attribute. In addition, the \code{"weight.trt"} attribute gives a numeric vector of treatment group weights. }
\item{data.name }{The character name of the data.}
%\item{.Random.seed }{ The random number seed for reproducibility.} 
\item{method }{A string of test method.}
}
\references{
Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.
}
\author{
Long Qu
}
%%\note{
%%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{permuteTrt}}
}
\examples{
set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
nparts( table(trt))		## 92378 partitions =  choose(20,10)/2
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e2L)		## use 500 random permutations
\dontrun{
pmat=permuteTrt(trt, 1e6L)		## use all partitions, as 1e6L >= 92378 
}
### Distance matrix interface
(test.rslt = mrpp.test(dist(x), trt, permutedTrt=pmat, weight.trt="df" ) )
p.value(test.rslt)
p.value(test.rslt, type="midp")

### Data matrix interface
mrpp.test(x, trt, permutedTrt=pmat, weight.trt="df" ) 

### Formula interface
dat=data.frame(x, trt)
fmla=as.formula(sprintf("cbind(\%s)~trt",paste("X",1:5,sep="",collapse=",")))
mrpp.test(fmla, dat, permutedTrt=pmat, weight.trt="df" )

}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}