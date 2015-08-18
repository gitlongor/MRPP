\name{mrpp.test}
\Rdversion{1.1}
\alias{mrpp.test.mrpp}
\alias{mrpp.test.dist}
\alias{mrpp.test.default}
\alias{mrpp.test.formula}
\alias{mrpp.test}
\alias{p.value}
\alias{p.value.mrpp.test}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
MRPP Test for One-way Design
}
\description{
Multiresponse permutation procedure (MRPP) for testing the equality of (multivariate) distributions among \eqn{K} samples
}
\usage{
mrpp.test(y, \dots)

\method{mrpp.test}{mrpp}(y, method="pearson3gca", eps=1e-8, \dots)
\method{mrpp.test}{default}(y, \dots)
\method{mrpp.test}{formula}(y, \dots)
\method{mrpp.test}{dist}(y, \dots)

\method{p.value}{mrpp.test}(x, type=c("raw","midp"), \dots)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{y}{A N*P data matrix, a \code{\link{mrpp}} object, a \code{\link{dist}} object of the data matrix, or a \code{\link{formula}}. 
}
  \item{method}{A character scalar of the form \dQuote{primary.adjustment} where \dQuote{.adjustment} part is optional. Partial string matching is allowed for boths parts. See details. 
}
  \item{eps}{A small non-negative number, differences below which among permuted test statistics are treated as equal.
}
  \item{\dots}{Additional arguments that must be \emph{named}. Except for \code{mrpp.test.mrpp} and \code{p.value.mrpp.test}, typically at least one of \code{"trt"} and \code{"permutedTrt"} arguments (cf \code{\link{mrpp}}) should be given. See details. 
}
  \item{x}{An object of class \code{"mrpp.test"}.
}
  \item{type}{A character scalar being one of \code{c("raw","midp")}, specifying the type of the p-value to be extracted from \code{x}.
}
}
\details{
For the \code{method} argument, it is a character scalar of the form \dQuote{primary.adjustment} where the \dQuote{.adjustment} part is optional. Choices of the \dQuote{primary} part are \code{c("pearson3", "pearson3gca", "gammagca", "permutation")}. This determines the primary method to compute p-values. \code{"pearson3"} uses Pearson type III distribution (see \code{\link{pearson3}}) that matches the first three moments with the exact permutation distribution of the MRPP statistic. \code{"pearson3gca"} further tunes the result of \code{"pearson3"} using the fourth moment through a Gram-Charlier type A (GCA) expansion (see \code{\link[PDQutils]{dapx_gca}}) around the Pearson type III distribution. \code{"gammagca"} first uses a Gamma distribution that matches the first two moments to obtain an initial approximation, and then uses a GCA expansion around this Gamma distribution to adjust for the third and the fourth moment of the exact permutation distribution of the MRPP statistic. Occasional negative approximate density values from GCA approximation are trunced as zero. \code{"permutation"} method performs the actual permutation test leading to a discrete p-value, which can be either exact (when \code{B} is larger than the number of non-equivalent permutations) or a Monte Carlo approximation to the p-value. 

Except when \dQuote{primary} is \code{"permutation"}, the \dQuote{adjustment} part of the \code{method} argument is a kernel name acceptable to \code{\link{dkernel}}. Examples are \code{"gaussian"} and \code{"biweight"}. When \dQuote{.adjustment} is non-missing, the approximate density of the \dQuote{primary} is adjusted by kernel smoothing approximation as in a kernel density estimation (KDE; see \code{\link{density}}), with the bandwidth selected such that the resulting KDE approximation and \dQuote{primary} approximation, when both evaluated at the permutation MRPP statistics, has minimum sum of squared differences. The KDE approximation always results in a bona fide probability density, even if the initial GCA approximation is not a density. 

When \dQuote{primary} is \code{"permutation"}, then either \dQuote{.adjustment} is missing or \dQuote{adjustment} is \code{"midp"}.  The former produces the actual permutation p-value, whereas the latter returns the mid-p-value. 

For the "\code{\dots}" argument, allowed options depend on the actual S3 method as follows
\tabular{ll}{
\code{mrpp.test.mrpp} and \code{p.value.mrpp.test} \tab \code{\dots} not used\cr
\code{mrpp.test.formula} \tab \code{"method"} and \code{"eps"} arguments as in \code{mrpp.test.mrpp}, plus any formal arguments of \code{\link{model.frame.default}} (except for the \code{formula} argument), and any formal arguments of \code{\link{mrpp}} (except for the \code{y} argument) \cr
\code{mrpp.test.dist} and \code{mrpp.test.default} \tab  \code{"method"} and \code{"eps"} arguments as in \code{mrpp.test.mrpp}, plus any formal arguments of \code{\link{mrpp}} (except for the \code{y} argument) \cr
}

Please note that the "\code{\dots}" arguments need to \emph{named} pairs. Unnamed arguments will not be properly passed down. 
}

\note{
The \code{mrpp.test.mrpp} method is the actual workhorse, and other methods are wrappers of it. Thus, it is recommended to construct a \code{\link{mrpp}} object first before invoking \code{mrpp.test} functions. 
}

\value{
\code{mrpp.test} methods return an object of class \code{c("mrpp.test", "htest")}, which is a list with the following components:
\item{statistic }{The observed test statistic.}
\item{all.statistics }{All permuted test statistics that have actually been used. If permutations have not been performed, this will be the same as the observed test statistic. }
\item{p.value }{The p-value. When \code{method="permutation"}, the \code{"midp"} attribute is set to the mid-p-value; when \code{method="permutation.midp"}, the \code{"raw"} attributes is set to the raw permutation p-value; in all other cases, no attributes are set.}
\item{parameter }{The number of permutations performed. The \code{"weight.trt"} attribute is set to a numeric vector of treatment group weights being used. The \code{"pdfmethod"} attribute is set to the "primiary" part of the \code{method} argument. TheA vector of the number of permutations and the weighting method. The latter is coded as an integer, with levels given in the \code{"weight methods"} attribute. The \code{"kernel"} attribute is set to the "adjustment" part of the \code{method} argument when appropriate, or \code{character(0)} otherwise. The \code{"exact"} attribute is \code{TRUE} is the p-value is the so-called exact p-value, namely the raw permutation p-value when \emph{all} non-equivalent permutatinos are enumerated. }
\item{data.name }{The character name of the data.}
%\item{.Random.seed }{ The random number seed for reproducibility.} 
\item{method }{A string describing the test method.}
}
\references{
Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.
}
\author{
Long Qu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{permuteTrt}}, \code{\link{mrpp}} object.
}
\examples{
set.seed(2340)
urand.bigz(0,seed=1032940L) # init seed
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
nparts( table(trt))		## 92378 partitions =  choose(20,10)/2
mrpp.obj=mrpp(x, trt, B=5e2L)	## use 500 random permutations
\dontrun{
mrpp.obj=mrpp(trt, B=Inf)		## use all partitions, as 1e6L >= 92378 
}

### "mrpp" object interface
mrpp.test(mrpp.obj)
mrpp.test(mrpp.obj, method="permutation") ## raw permutation p-value
mrpp.test(mrpp.obj, method="pearson3") ## Pearson type III approximate p-value
mrpp.test(mrpp.obj, method="gammagca") ## GCA approx around Gamma 
	## GCA approx around Pearson type III, followed by Gaussian kernel smoothing
mrpp.test(mrpp.obj, method="pearson3gca.gaussian") 

(test.rslt = mrpp.test(mrpp.obj, method="permutation.midp")) ## permutation mid-p-value
p.value(test.rslt)  ## default to report raw p-value
p.value(test.rslt, type="midp")
p.value(test.rslt, type="raw")

### Distance matrix interface
mrpp.test(dist(x), trt=trt, B=1e3L, method="permutation")
pmat = mrpp.obj$permutedTrt.env$permutedTrt
mrpp.test(dist(x), permutedTrt=pmat, method="permutation") ## using precomputed permutedTrt

### Data matrix interface
mrpp.test(x, trt=trt, method="pearson3gca") 

### Formula interface
dat=data.frame(x, trt)
(fmla=as.formula(sprintf("cbind(\%s)~trt",paste("X",1:5,sep="",collapse=","))))
mrpp.test(fmla, data=dat, permutedTrt=pmat, weight.trt="df", method='permutation' )  


}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}