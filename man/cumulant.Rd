\name{moments}
\alias{mean}
\alias{var}
\alias{sd}
\alias{skewness}
\alias{kurtosis}
\alias{moment}
\alias{cumulant}
\alias{mean.mrpp}
\alias{var.mrpp}
\alias{sd.mrpp}
\alias{skewness.mrpp}
\alias{kurtosis.mrpp}
\alias{moment.mrpp}
\alias{cumulant.mrpp}
\alias{var.default}
\alias{sd.default}
\alias{skewness.default}
\alias{kurtosis.default}
\alias{moment.default}
\alias{cumulant.default}

\title{
Moments and Cumulants of MRPP Statistics
}
\description{
Generic functions and methods for computing moments and cumulants for class \code{mrpp} objects, as well as default methods for numeric vectors.
}
\usage{
\method{mean}{mrpp}(x, \dots)

var(x, \dots)
\method{var}{default}(x, y = NULL, na.rm = FALSE, use, \dots) 
\method{var}{mrpp}(x, \dots)

sd(x, \dots)
\method{sd}{default}(x, na.rm = FALSE, \dots) 
\method{sd}{mrpp}(x, \dots)

skewness(x, \dots)
\method{skewness}{default}(x, na.rm = FALSE, \dots) 
\method{skewness}{mrpp}(x, \dots)

kurtosis(x, excess=TRUE, \dots )
\method{kurtosis}{default}(x, excess = TRUE, na.rm = FALSE, \dots )
\method{kurtosis}{mrpp}(x, excess=TRUE, \dots)


moment(x, order, central, \dots )
\method{moment}{default}(x, order=1:4, central=FALSE,  \dots )
\method{moment}{mrpp}(x, order=1:4, central=FALSE, \dots)

cumulant(x, order,  \dots )
\method{cumulant}{default}(x, order=1:4,   \dots )
\method{cumulant}{mrpp}(x, order=1:4,  \dots)
}

\arguments{
  \item{x}{
	An object whose moments or cumulants to be computed. 
}
  \item{y, na.rm, use}{
	The same as in \code{\link[stats:var]{stats::var}} and/or \code{\link[stats:sd]{stats::sd}}. 
}
  \item{excess}{
	Logical, the same as in \code{\link[moments:kurtosis]{moments::kurtosis}}.
}
  \item{central}{
	Logical, the same as in \code{\link[moments:moment]{moments::moment}}.
}
  \item{order}{
  A \emph{vector} of non-negative integers, specifying the order of moments or cumulants.
}
  \item{\dots}{
	Placeholder for individual methods.
}
}

\details{
S3 methods \code{mean.mrpp}, \code{var.mrpp}, \code{sd.mrpp}, \code{skewness.mrpp} and \code{kurtosis.mrpp} compute the mean, variance, standard deviation, skewness and kurtosis for MRPP statistics of \code{"mrpp"} objects, respectively, with respective to the permutation distribution. The implementation does not require enumerating all MRPP statistics. 

\code{var} and \code{sd} are S3 generic functions, with \code{default} methods being \code{stats::var} and \code{stats::sd}, respectively. 

\code{skewness} and \code{kurtosis} are S3 generic functions, with \code{default} methods being \code{moments::skewness} and \code{moments::kurtosis}, respectively. 

\code{moment} and \code{cumulant} are S3 generic functions, with default methods computing the moments and cumulants on numeric vector \code{x}. 

}
\value{
\code{mean}, \code{var}, \code{sd}, \code{skewness}, and \code{kurtosis} returns a single numeric scalar. 

\code{moment} and \code{cumlant} return a vector of the same length as the \code{order} argument, with names being \code{as.character(order)}. 
}

\note{
\code{moment.mrpp} and \code{cumulant.mrpp} with \code{max(order)>=5L} have not been implemented yet. 

}
\references{
Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.

Jack Siemiatycki. (1978) Mantel's Space-Time Clustering Statistic: Computing Higher Moments and a Comparison of Various Data Transforms. \emph{J Statist Comput Simul}. \bold{7}, 13--31.

D. S. Tracy and I. H. Tajuddin. (1985) Extended Moment Results for Improving Inferences Based on MRPP. \emph{Commun Statist - Theor Meth}. \bold{14}(6), 1485--1496.
}
\author{
Long Qu
}

\seealso{
\code{\link[stats:var]{var}}, \code{\link[stats:sd]{sd}} in \code{stats} package. 

\code{\link[moments:skewness]{skewness}}, \code{\link[moments:kurtosis]{kurtosis}}, and \code{\link[moments:all.moments]{all.moments}} in the \code{moments} package.

\code{\link[PDQutils:moment2cumulant]{moment2cumulant}} and \code{\link[PDQutils:cumulant2moment]{cumulant2moment}} in the \code{PDQutils} package. 

\code{\link{mrpp}} object. 
}

\examples{
set.seed(2340)
n=5
x=matrix(rnorm(n*5*2),ncol=5)
apply(x, 2, moment)
apply(x, 2, cumulant)
apply(x, 2, mean)
apply(x, 2, sd)
apply(x, 2, skewness)
apply(x, 2, kurtosis)

trt=gl(2,n)
mrpp.obj = mrpp(x, trt, B=Inf)
mrpp.tst = mrpp.test(mrpp.obj, method="permutation")

mean(mrpp.obj) - mean(mrpp.tst$all.statistics)
var(mrpp.obj) - mean((mrpp.tst$all.statistics-mean(mrpp.tst$all.statistics))^2)
cumulant(mrpp.obj, order=3) - mean((mrpp.tst$all.statistics-mean(mrpp.tst$all.statistics))^3)
moment(mrpp.obj, order=4, central=TRUE) - moment(mrpp.tst$all.statistics, order=4, central=TRUE)

moment(mrpp.obj, order=1:4, central=FALSE)
moment(mrpp.obj, order=1:4, central=TRUE)
cumulant(mrpp.obj, order=1:4)

skewness(mrpp.obj)
kurtosis(mrpp.obj)

}
\keyword{ univar }

