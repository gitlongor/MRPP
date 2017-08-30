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
bw.smoothp(y, start.bw = NULL, kernel='triweight', 
	bw.method='mse(z[1])', scale='raw', verbose=TRUE, n.subset, 
	mrpp.stats=NULL, r=seq_len(y$R))
bw.range(x, length=50L, lower=.05, upper=.95, safety=2)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{y}{An \code{mrpp} object. 
}
  \item{start.bw}{An optional numeric value, being the starting value of the optimization for searching for an optimal bandwidth. 
}
  \item{kernel}{A character scalar indicating the kernel to be used. The allowed values are given in \code{\link{.smooth1.kernels}}. 
}
  \item{bw.method}{A character scalar specifying the bandwidth selection method. See details. 
}
  \item{scale}{Either \code{"raw"} or \code{"log"}. This specifies whether derivative approximation to discrete changes is minimized on the original p-value scale or it log scale. This is only effective when \code{bw.method} is one of \code{c('sym1','drop1','add1','keep1','dropadd1','dropaddsym1')}.
}
  \item{verbose}{Logical. When \code{TRUE}, diagnostic plots will be produced. 
}
  %\item{permutedTrt}{A \code{"permutedTrt"} object from calling \code{\link{permuteTrt}}.
%}
  \item{r}{An integer vector, indicating the columns (variables) of \code{y} that needs to be considered. 
}
  \item{n.subset}{An integer, only effective when \code{bw.method} is one of \code{c('sym1','drop1','add1','keep1','dropadd1','dropaddsym1','ss.gradp')}. When \code{n.subset} is at least \code{length(r)}, all variables will be used. When \code{n.subset} is less than \code{length(r)}, then the search for optimal bandwidth will only be performed on a subset of variables to save computational time.  If \code{n.subset} is missing, then the default behavior is to use all \code{r} variables when \code{length(r)<=100L}; otherwise, a subset of 100 plus 20\% additional variables will be used.
}

%  \item{distFunc}{A function that returns a \code{"dist"} object, as in the \code{\link{dist}} function. 
%}

%  \item{\dots}{Additional arguments passed on to \code{\link{mrpp.test.dist}} or additional arguments for \code{bw.method="mse(z[1])"}. See details.
%}
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
For \code{bw.method},  currently supported values are as follows: 

\itemize{
\item\code{bw.method="drop1"}: The bandwidth best approximating MRPP mid-p-values from dropping one variable at a time will be selected. This is equivalent to a backward finite difference approximation to the derivative. 

\item\code{bw.method="add1"}: The bandwidth best approximating MRPP mid-p-values from adding one more redundant variable will be selected. This is equivalent to a forward finite difference approximation to the derivative. 

\item\code{bw.method="sym1"}: The bandwidth selected is such that the resulting derivative will best approximate two-sided symmetric finite difference. 

\item\code{bw.method="dropadd1"}: The bandwidth selected is such that the resulting derivative will minimize the sum of squared approximation errors of both backward and forward finite difference. 

\item\code{bw.method="dropaddsym1"}: The bandwidth selected is such that the resulting derivative will minimize the sum of squared approximation errors of backward,  forward and symmetric finite difference. 

\item\code{bw.method="keep1"}: The bandwidth best approximating univariate MRPP mid-p-values will be selected. 

\item\code{bw.method='pearson3'}: The bandwidth selected is such that the resulting kernel density estimate best approximates an initial density. The initial density is a Pearson type III distribution, whose first 3 moments match those of MRPP statistics exactly. See also \code{\link{dpearson3}}.

\item\code{bw.method='pearson3gca'}: The bandwidth selected is such that the resulting kernel density estimate best approximates an initial density. The initial density is the Gram-Charlier type A approximation (GCA) to the Pearson type III distribution. The first 4 moments of the GCA are expected to match those of MRPP statistics. See also \code{\link{dpearson3gca}}. 

\item\code{bw.method="ss.gradp"}: The bandwidth that maximize the sum of squared derivatives of density of smoothed p-values will be selected.

\item\code{bw.method="mse(z[1])"}: The bandwidth selection problem for approximating permutation statistics is treated as if finding the bandwidth in kernel \dQuote{density} estimation and an asymptotically optimal bandwidth minimizing mean squared error at the first data point will be used. The procedure starts from plug-in estimates of bandwidths to estimate \eqn{f(x_1)} and \eqn{f''(x_1)} using \code{\link{hpi}}. The final optimal bandwidth will be computed as 
\deqn{
\left[\frac{R f(x_1)}{4N \lbrace k_2 f''(x_1) \rbrace^2} \right]^{1/5}
}{
( {R f(x_1)} / [4 N {k_2 f''(x_1)}^2] )^{1/5}
}
where \code{R} and \code{k_2} are constants determined by the \code{kernel} under use (c.f. \code{\link{krRkernel}}), \code{N} is the number of permutation statistics and \code{f} and \code{f''} are the estimates of the density and its second derivative. This formula corrects a typo in Parzen (1962). 
}

%Details of the \code{\dots} argument are as follows: 
%\describe{
%  \item{iter.max}{The max number of iterations. For an iterative procedure, set this to be  greater than 1L.
%}  
%  \item{eps}{A small positive number such that iterations are terminated when the difference between successive bandwidths is below \code{eps}.
%}
%  \item{start.bw}{A function, a character, or a positive number. If being a function or a character, the corresponding function will be called with \code{x} as the only argument to get the initial bandwidth. When it is a number, it will be used as the starting bandwidth directly.
%}

%}

\code{bw.range} is currently implemented as 
\preformatted{
bw.range=function(x, length=50L, lower=.05, upper=.95, safety=2)
{
	uniqstat=sort(unique(x))
	lo = quantile(diff(uniqstat), lower)
	hi = diff(quantile(uniqstat, prob=c(lower, upper)))
	lf = log10(safety)
	10^seq.int(from=log10(lo)-lf, to=log10(hi)+lf, length.out=length)
}
}
}

\value{
A numeric positive scalar, giving the final bandwidth. 
}
\references{
Parzen, E. (1962). On Estimation of a Probability Density Function and Mode. The Annals of Mathematical Statistics 33: 1065--1076.

Scott, D. W. (1992). Multivariate Density Estimation: Theory, Practice, and Visualization. Wiley, New York.

Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall, London. 

Wand, M.P. and Jones, M.C. (1994) Multivariate plugin bandwidth selection. Computational Statistics. 9, 97-116.
}
\author{Long Qu
}
\seealso{
\code{\link{dkde}}, \code{\link{dpearson3}}, \code{\link{mrpp}}, \code{\link{permuteTrt}}, \code{\link{hpi}}, \code{\link{krRkernel}}
}
\examples{
set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
mrpp.obj=mrpp(x, trt, B=5e3L)	## use 5000 random permutations
bw.smoothp(mrpp.obj, bw.method='drop1', verbose=interactive())
bw.smoothp(mrpp.obj, bw.method='add1', verbose=interactive())
bw.smoothp(mrpp.obj, bw.method='sym1', verbose=interactive())
bw.smoothp(mrpp.obj, bw.method='keep1', verbose=interactive())
bw.smoothp(mrpp.obj, bw.method='pearson3', verbose=interactive())
bw.smoothp(mrpp.obj, bw.method='pearson3gca', verbose=interactive())
bw.smoothp(mrpp.obj, bw.method='mse(z[1])', verbose=interactive())

mrpp.tst=mrpp.test(mrpp.obj, method='permutation')
summary( bw.range(mrpp.tst$all.statistics) )


}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{distribution}
\keyword{smooth}
