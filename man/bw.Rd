\name{bandwidth}
\Rdversion{1.1}
%\alias{bw.mse.pdf.asym}
\alias{bw.smoothp}
\alias{bw.range}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
Bandwidth Selection for Kernel Approximation to Discrete Cumulative Distribution Functions
}
\description{
\code{bw.smoothp} chooses bandwidth based on the derivative of smoothed permutation p-values. 

\code{bw.range} returns a vector of candidate bandwidths that cover a large range of reasonable values. 
}
\usage{
%bw.mse.pdf.asym(x,iter.max=1L,eps=1e-6,start.bw=bw.nrd, verbose=FALSE)
bw.smoothp(y, permutedTrt, r=seq_len(NCOL(y)), bw = NULL, 
	distFunc=dist,  kernel='gaussian', %weight.trt="df", 
	method=c('drop1','keep1','dropkeep1','ss.gradp','kde.mse1'), verbose=FALSE, \dots)
bw.range(x, length=200, lower=.05, upper=.95, safety=100)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{y}{A data matrix with rows being observations and columns being variables. 
}
  \item{permutedTrt}{A \code{"permutedTrt"} object from calling \code{\link{permuteTrt}}.
}
  \item{r}{An integer vector, indicating the columns (variables) of \code{y} that needs to be considered. 
}
  \item{bw}{Either a vector of positive real bandwidths to be searched from, or \code{NULL}. If \code{NULL}, then \code{bw.range} will be called to obtain a vector of candidate of bandwidths.
}
  \item{distFunc}{A function that returns a \code{"dist"} object, as in the \code{\link{dist}} function. 
}
  \item{kernel}{A character scalar indicating the kernel to be used. Allowed values are given at \code{\link{pkernel}}.
}
%  \item{weight.trt}{As in \code{\link{mrpp.test.dist}}.}

  \item{method}{A character scalar specifying the bandwidth selection method. See details. 
}
  \item{verbose}{Logical. When \code{TRUE}, the bandwidth for each iteration will be printed for \code{method="kde.mse0"}, and diagnostic plots will be produced for other methods. 
}
  \item{\dots}{Additional arguments passed on to \code{\link{mrpp.test.dist}} or additional arguments for \code{method="kde.mse1"}. See details.
}
  \item{x}{A vector discrete statistics.
}
  \item{length}{A positive integer, specifying the number of candidate bandwidths to be searched from. 
}
  \item{lower, upper}{Numeric values between 0 and 1, indicate the quantiles. See details. 
}
  \item{safety}{A positive number, indicating the safety zone to be expanded. See details. 
}

}
\details{
For \code{bw.smoothp},  if \code{method="drop1"}, the bandwidth best approximating MRPP mid-p-values from dropping one variable at a time will be selected. 

If \code{method="keep1"}, the bandwidth best approximating univariate MRPP mid-p-values will be selected. 

If \code{method="dropkeep1"}, the geometric mean from bandwidths based on \code{method="drop1"} and \code{method="keep1"} will be selected. 

If \code{method="ss.gradp"}, the bandwidth that maximize the sum of squared derivatives of smoothed p-values will be selected.

If \code{method="kde.mse1"}, the bandwidth selection problem for approximating permutation statistics is treated as if finding the bandwidth in kernel \dQuote{density} estimation and an asymptotically optimal bandwidth minimizing mean squared error at the first data point will be used. In this case, the \code{\dots} argument can be named pairs with names being any subset of \code{c("iter.max","eps","start.bw")}. The procedure starts from the starting bandwidth specified by \code{start.bw}. The new bandwidth will be computed as 
\deqn{
\left[\frac{Rf(x_1)}{N\lbrace f''(x_1) \rbrace^2} \right]^{1/5}
}{
( {R f(x_1)} / [N {f''(x_1)}^2] )^{1/5}
}
where \code{R} is approximately 0.2821 if \code{kernel="gaussian"} (c.f. \code{\link{density}} with \code{give.Rkern=TRUE}), \code{N} is the number of permutation statistics and \code{f} is the density estimate using the current bandwidth. 

Details of the \code{\dots} argument are as follows: 
\itemize{
  \item{iter.max}{The max number of iterations. For an iterative procedure, set this to be  greater than 1L.
}  
  \item{eps}{A small positive number such that iterations are terminated when the difference between successive bandwidths is below \code{eps}.
}
  \item{start.bw}{A function, a character, or a positive number. If being a function or a character, the corresponding function will be called with \code{x} as the only argument to get the initial bandwidth. When it is a number, it will be used as the starting bandwidth directly.
}

}

\code{bw.range} is currently implemented as 
\preformatted{
bw.range=function(x, length=200, lower=.05, upper=.95, safety=100)
{
	uniqstat=sort(unique(x))
	lo = quantile(diff(uniqstat), lower)
	hi = diff(quantile(uniqstat, prob=c(lower, upper)))
	lf = log10(safety)
	10^seq(log10(lo)-lf, log10(hi)+lf, length=length)
}
}
}
%\note{
%	When \code{method="kde.mse1"}, only \code{"gaussian"} and \code{"biweight"} are allowed for the \code{kernel} argument.
%}

\value{
A numeric positive scalar, giving the final bandwidth. 
}
\references{
Parzen, E. (1962). On Estimation of a Probability Density Function and Mode. The Annals of Mathematical Statistics 33: 1065--1076.

Scott, D. W. (1992). Multivariate Density Estimation: Theory, Practice, and Visualization. Wiley, New York.

Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall, London. 
}
\author{Long Qu
}
\seealso{
\code{\link{density}}, \code{\link[ks]{kdde}}
}
\examples{
set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e3L)		## use 5000 random permutations
bw.smoothp(x, pmat, method='drop1', verbose=interactive())
bw.smoothp(x, pmat, method='keep1', verbose=interactive())
bw.smoothp(x, pmat, method='dropkeep1', verbose=interactive())
bw.smoothp(x, pmat, method='ss.gradp', verbose=interactive())
bw.smoothp(x, pmat, method='kde.mse1', verbose=interactive())

mrpp=mrpp.test(x~trt, permutedTrt=pmat, method='permutation')
summary( bw.range(mrpp$all.statistics) )


}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{distribution}
\keyword{smooth}
