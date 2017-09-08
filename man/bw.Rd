\name{bandwidth}
\Rdversion{1.1}
%\alias{bw.amse.pdf}
\alias{bw.kdep}
\alias{bw.range}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
Bandwidth Selection for Kernel Approximation to Discrete Cumulative Distribution Functions
}
\description{
\code{bw.kdep} chooses bandwidth based on the derivative of smoothed permutation p-values. 

\code{bw.range} returns a vector of candidate bandwidths that cover a large range of reasonable values. 
}
\usage{
%bw.amse.pdf(x,iter.max=1L,eps=1e-6,start.bw=bw.nrd, verbose=FALSE)
bw.kdep(y, start.bw = NULL, kernel='triweight', 
	bw.method='amse(z[1]).pearson3', scale='raw', verbose=TRUE, n.subset, 
	mrpp.stats=NULL, r=seq_len(y$R))
bw.range(x, kernel, length=15L, var.infl=0.1)
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

%  \item{\dots}{Additional arguments passed on to \code{\link{mrpp.test.dist}} or additional arguments for \code{bw.method="amse(z[1])"}. See details.
%}
  \item{x}{A vector discrete statistics.
}
  \item{length}{A positive integer, specifying the number of candidate bandwidths to be searched from. 
}
%  \item{lower, upper}{Numeric values between 0 and 1, indicate the quantiles. See details. 
%}
%  \item{safety}{A positive number, indicating the safety zone to be expanded. See details. 
%}
  \item{var.infl}{A positive numeric value. The variance of the estimated distribution using the kernel method is always larger than the sample variance by a quantity proportional to squared bandwidth. \code{var.infl} specifies the maximal percentage of such inflation of variance. By providing an upper bound of such inflation of variance, a corresponding upper bound of sensible bandwidths can be found. With the default 10\% variance inflation, the difference between the estimated kernel density (oversmoothed) and a fine histogram of data is usually visually discernable and visually unappealling. 
}
}
\details{
For \code{bw.method},  currently supported values are as follows: 

\itemize{
\item\code{bw.method="drop1"}: The bandwidth best approximating MRPP mid-p-values from dropping one variable at a time will be selected. This is equivalent to a backward finite difference approximation to the derivative. 

\item\code{bw.method="add1"}: The bandwidth best approximating MRPP mid-p-values from adding one more redundant variable will be selected. This is equivalent to a forward finite difference approximation to the derivative. 

\item\code{bw.method="sym1"}: The bandwidth selected is such that the resulting derivative will best approximate two-sided symmetric finite difference. 

%\item\code{bw.method="dropadd1"}: The bandwidth selected is such that the resulting derivative will minimize the sum of squared approximation errors of both backward and forward finite difference. 

%\item\code{bw.method="dropaddsym1"}: The bandwidth selected is such that the resulting derivative will minimize the sum of squared approximation errors of backward,  forward and symmetric finite difference. 

\item\code{bw.method="keep1"}: The bandwidth best approximating univariate MRPP mid-p-values will be selected. 

\item\code{bw.method='pearson3'}: The bandwidth selected is such that the resulting kernel density estimate best approximates an initial density. The initial density is a Pearson type III distribution, whose first 3 moments match those of MRPP statistics exactly. See also \code{\link{dpearson3}}.

\item\code{bw.method='tlnorm'}: The bandwidth selected is such that the resulting kernel density estimate best approximates an initial density. The initial density is affine transformed log-normal distribution, whose first 3 moments match those of MRPP statistics exactly. See also \code{\link{dtlnorm}}.

\item\code{bw.method='pearson3gca'}: The bandwidth selected is such that the resulting kernel density estimate best approximates an initial density. The initial density is the Gram-Charlier type A approximation (GCA) to the Pearson type III distribution. The first 4 moments of the GCA are expected to match those of MRPP statistics. See also \code{\link{dpearson3gca}}. 

\item\code{bw.method='p3tlnormmix'}: The bandwidth selected is such that the resulting kernel density estimate best approximates an initial density. The initial density is a mixture of Pearson type III distribution and an affine transformed log-normal distribution. The first 4 moments of the mixture are expected to match those of MRPP statistics. See also \code{\link{dpearson3}} and \code{\link{dtlnorm}}. 

%\item\code{bw.method="ss.gradp"}: The bandwidth that maximize the sum of squared derivatives of density of smoothed p-values will be selected.

\item\code{bw.method="amise.WJ"}: The bandwith selection problem for approximating permutation statistics is treated as if finding the bandwidth in kernel \dQuote{density} estimation and an asymptotically optimal bandwidth minimizing asymptotic mean integrated squared error (AMISE) will be used. The procedure uses Wand and Jones' plug-in estimates of bandwidths computed \code{\link{hpi}}, which is based on the \code{"gaussian"} kernel. The result is then converted to the requested kernel.

\item\code{bw.method="amise.pearson3"}: The bandwith selection problem for approximating permutation statistics is treated as if finding the bandwidth in kernel \dQuote{density} estimation and an asymptotically optimal bandwidth minimizing asymptotic mean integrated squared error (AMISE) will be used. The procedure first assumes a reference Pearson type III distribution that matches the first three moments exactly to estimate \eqn{R(f'')}{R(f'')}, where \code{f''} is the second derivative of the assumed density of test statistics. If this value is non-finite, the Wand and Jones' plug-in estimates of bandwidths computed by \code{\link{hpi}} will be used. The result is then plugged into formula

\deqn{\lbrace \frac{R(K)}{R(f'') N V(K)} \rbrace^{1/5}}{( (R(K)/R(f'')/N/V(K))^(1/5)}

where \code{R(K)} is the \emph{roughness} of kernel \code{K} computed by \code{\link{Rkernel}}, \code{V(K)} is the variance of kernel \code{K} computed by  \code{\link{mkernel}}, and \code{N} is the number of permutation statistics.

\item\code{bw.method="amse(z[1]).WJ"}: The bandwidth selection problem for approximating permutation statistics is treated as if finding the bandwidth in kernel \dQuote{density} estimation and an asymptotically optimal bandwidth minimizing asymptotic mean squared error (AMSE) at the first data point will be used. The procedure starts from Wand and Jones' plug-in estimates of bandwidths to estimate \eqn{f(x_1)} and \eqn{f''(x_1)} using \code{\link{hpi}}, with \code{'epanechnikov'} and \code{'triweight'} kernels used, respectively. The final optimal bandwidth will be computed as 
\deqn{
\left[\frac{R(K) f(x_1)}{N \lbrace V(K) f''(x_1) \rbrace^2} \right]^{1/5}
}{
( {R(K) f(x_1)} / [ N {V(K) f''(x_1)}^2] )^{1/5}
}
where \code{R(K)} is the \emph{roughness} of kernel \code{K} computed by \code{\link{Rkernel}}, \code{V(K)} is the variance of kernel \code{K} computed by  \code{\link{mkernel}}, \code{N} is the number of permutation statistics and \code{f} and \code{f''} are the estimates of the density and its second derivative. This formula corrects a typo in Parzen (1962). 


\item\code{bw.method="amse(z[1]).pearson3"}: The bandwidth selection problem for approximating permutation statistics is treated as if finding the bandwidth in kernel \dQuote{density} estimation and an asymptotically optimal bandwidth minimizing asymptotic mean squared error (AMSE) at the first data point will be used. The procedure starts from Wand and Jones' plug-in estimates of bandwidths to estimate \eqn{f(x_1)} and \eqn{f''(x_1)} the method described for the \code{amise.pearson3} method, with \code{'epanechnikov'} and \code{'triweight'} kernels used, respectively. The final optimal bandwidth will be computed as 
\deqn{
\left[\frac{R(K) f(x_1)}{N \lbrace V(K) f''(x_1) \rbrace^2} \right]^{1/5}
}{
( {R(K) f(x_1)} / [ N {V(K) f''(x_1)}^2] )^{1/5}
}
where \code{R(K)} is the \emph{roughness} of kernel \code{K} computed by \code{\link{Rkernel}}, \code{V(K)} is the variance of kernel \code{K} computed by  \code{\link{mkernel}}, \code{N} is the number of permutation statistics and \code{f} and \code{f''} are the estimates of the density and its second derivative. This formula corrects a typo in Parzen (1962). 
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
bw.kdep(mrpp.obj, bw.method='drop1', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='add1', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='sym1', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='keep1', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='pearson3', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='tlnorm', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='p3tlnormmix', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='pearson3gca', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='amse(z[1]).WJ', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='amse(z[1]).pearson3', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='amise.WJ', verbose=interactive())
bw.kdep(mrpp.obj, bw.method='amise.pearson3', verbose=interactive())

mrpp.tst=mrpp.test(mrpp.obj, test.method='permutation')
summary( bw.range(mrpp.tst$all.statistics, kernel='biweight') )


}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{distribution}
\keyword{smooth}
