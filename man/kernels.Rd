\name{kernel}
\Rdversion{1.1}
\alias{skernel}
\alias{pkernel}
\alias{dkernel}
\alias{ddkernel}
\alias{d2dkernel}
\alias{fourier.kernel}
\alias{ckernel}
\alias{mkernel}
\alias{pkde}
\alias{dkde}
\alias{.kernels}
\alias{.smooth0.kernels}
\alias{.smooth1.kernels}
\alias{.smooth2.kernels}
\alias{Rkernel}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
Support Functions For Univariate Kernel Densities
}
\description{
\code{pkernel}, \code{dkernel} and \code{ckernel} return the cumulative distribution function, probability density function and characteristic function, respectively, for common univariate kernels specified by names. \code{fourier.kernel} is a synonym of \code{ckernel}. 

\code{skernel} returns the support of the requested kernel. \code{ddkernel} and \code{d2dkernel} return the 1st and the 2nd derivative of the density of the requested kernel, respectively.
\code{mkernel} returns the raw moments of the requested kernel distribution. 

\code{pkde} and \code{dkde} return the \emph{fitted} cumulative distribution function and probability density function. 

\code{Rkernel} computes roughness of (derivatives) of kernel function. 

\code{.kernels} is a vector containing the names of the kernel functions currently supported. 
\code{.smooth0.kernels} is a vector containing the names of the supported kernel functions that are continuous over the real line. 
\code{.smooth1.kernels} is a vector containing the names of the supported kernel functions of which the first derivatives are continuous over the real line. 
\code{.smooth2.kernels} is a vector containing the names of the supported kernel functions of which the second derivatives are continuous over the real line. 
}
\usage{
pkernel(kernel= .kernels)
dkernel(kernel= .kernels)
ddkernel(kernel=.kernels)
d2dkernel(kernel=.kernels)
ckernel(kernel= .kernels, root.2pi=FALSE, const=NULL)
fourier.kernel(kernel= .kernels, root.2pi=FALSE, const=NULL)

skernel(kernel=.kernels, eps=1e-3)
mkernel(kernel=.kernels, order=1:4, R=TRUE)

pkde(x, bw=bw.nrd, kernel=.kernels)
dkde(x, bw=bw.nrd, kernel=.kernels,	from, to, neval, eps=1e-3)
Rkernel(kernel=.kernels)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{kernel}{A character scalar of the name of the kernel. Partial string matching is allowed. For \code{pkde}, this can also be a valid univariate cumulative distribution function. 
}
  \item{root.2pi}{A logical scalar. When \code{TRUE}, a factor of \code{1/sqrt(2*base::pi)} will be multiplied to the integral transform; otherwise, the integral transform itself is computed. 
}
  \item{const}{Either \code{NULL} or a character vector with elements among \code{c('exponent','coefficient')}. This specifies whether the characteristic exponent and/or the characteristic coefficient of the transform will be added to the result as an attribute. See Parzen (1962) for their definitions. 
}
  \item{eps}{A small non-negative number. The effective support reported exclude regions with density values below \code{eps}. If \code{eps=0},  non-compact kernels like \code{gaussian} will have a reported support being the whole real line. 
}
  \item{order}{An integer vector, being a subvector of \code{1:4}. Orders larger than 4 are not yet supported. This could be a \code{0}-length vector, in which case \code{R=TRUE} will be assumed. 
}
  \item{R}{A logical scalar. If \code{TRUE}, then the integral of square of the density will be reported. Cf. the \code{give.Rkern} argument of \code{\link{density}}.
}
  \item{x}{A numeric vector of data values.
}
  \item{bw}{A positive vector indicating the bandwidth to be used in kernel density estimation. Alternatively, this can be a function that computes the bandwidth, or the character name of such a function. 
}
  \item{from,to}{Numeric scalars, the left and right boundaries within which the estimated density is to be evaluated, as in \code{\link{density}}. When missing, they are determined by extending the minimum and maximum of \code{x} by \code{s*bw}, where \code{s} is the effective support boundary reported by \code{skernel(kernel, eps)[2L]}.
}
  \item{neval}{An positive integer, the number of evaluation points, as in \code{\link{density}}. When missing, this is chosen as the minimum of 1024 and \code{nextn(2*lenght(x),2L)}.
}
}
\format{
Currently, \code{.smooth2.kernels} is the character vector

\code{c("gaussian", 'normal', 'triweight', 'tricube',"logistic",'sech')}.

\code{.smooth1.kernels} additionally includes \code{'biweight', 'quartic'}. 

\code{.smooth0.kernels} further includes \code{'epanechnikov','parabolic', 'triangular','cosine'}. 

And \code{.kernels} includes all above, plus \code{'uniform','rectangular'}. 


High-order kernels that are not bona fide distributions are not allowed. 

In the above \code{"uniform"} and \code{"rectangular"} are synonymous; \code{'gaussian'} and \code{'normal'} are synonymous; \code{'epanechnikov'} and \code{'parabolic'} are synonymous; \code{'biweight'} and \code{'quartic'} are synonymous. 

Whenever the kernel has compact support, the support will be \code{[-1,1]}; thus, the variance of each kernel is not necessarily one, which behaves differently than \code{\link{density}}.


}

\note{
The \code{pkde} implementation does \emph{not} make use of the fast Fourier transform (FFT); but the \code{dkde} uses FFT. 
}
\value{
For \code{pkernel}, \code{dkernel}, \code{ckernel}, \code{ddkernel} and \code{d2dkernel}, the result is a function (closure) that accepts a numeric vector as its argument, with formal arguments being \code{q}, \code{x}, \code{s}, \code{x} and \code{x}, respectively.

For \code{dkde}, if \code{length(bw)==1}, then the result is a function; otherwise, the result is a list of functions, corresponding to each bandwidth. In either case, the \code{attr(*, "environment")} attribute is set to the evaluation frame of \code{dkde}, allowing the access of some intermediate objects. 

For \code{mkernel}, the result is a numeric vector of the requested moments. 

For \code{skernel}, the result is a numeric vector of length 2, being the lower and upper bound of the (effective) support. 

For \code{ckernel} and \code{fourier.kernel}, when the \code{const} is set to non-\code{NULL}, an additional attribute \code{const} will be set to a numeric vector requested. 

\code{Rkernel} returns a named numeric vector of the same length as \code{deriv.order}, with names being \code{as.character(deriv.order)}. This is the integral of the square of the \code{deriv.order}th derivative of the density of the kernel.
}

\references{
\url{https://en.wikipedia.org/wiki/Kernel_(statistics)}

Section 3.5 of Bernard W. Silverman (1986) Density Estimation for Statistics and Data Analysis. London; Newyork: Chapman & Hall. 

Emanuel Parzen (1962) On Estimation of a Probability Density Function and Mode. Annals of Mathematical Statistics, 33(3), 1065-1076.

}
\author{Long Qu
}
\seealso{
\code{\link{density}}
}
\examples{
print(.kernels)
f=dkernel('biweight')
F=pkernel('biweight')
f(0) ## 15/16 = 0.9375
F(0) ## 0.5
f(-2:2)
F(-2:2)
ft=fourier.kernel('biweight')
ft(-2:2)
sapply(.kernels, mkernel)
sapply(.kernels, skernel, eps=0)
sapply(.kernels, skernel, eps=1e-5)

if(interactive()){
	curve(dkernel('gau')(x),-2,2, ylim=c(0,1.1))
	curve(dkernel('biw')(x),-2,2, add=TRUE, col=2)
	curve(dkernel('triw')(x),-2,2, add=TRUE, col=3)
	curve(dkernel('tric')(x),-2,2, add=TRUE, col=4)
	curve(dkernel('logi')(x),-2,2, add=TRUE, col=5)
}
if(interactive()){
	curve(pkernel('gau')(x),-2,2, ylim=c(0,1))
	curve(pkernel('biw')(x),-2,2, add=TRUE, col=2)
	curve(pkernel('triw')(x),-2,2, add=TRUE, col=3)
	curve(pkernel('tric')(x),-2,2, add=TRUE, col=4)
	curve(pkernel('logi')(x),-2,2, add=TRUE, col=5)
}
if(interactive()){
	curve(fourier.kernel('gau')(x),-20,20, ylim=c(-.05,1/sqrt(2*base::pi)))
	curve(fourier.kernel('biw')(x),-20,20, add=TRUE, col=2)
	curve(fourier.kernel('triw')(x),-20,20, add=TRUE, col=3)
	curve(fourier.kernel('tric')(x),-20,20, add=TRUE, col=4)
	curve(fourier.kernel('logi')(x),-20,20, add=TRUE, col=5)
}
set.seed(23425)
x=rnorm(100)
Fhat=pkde(x)
Fhat=pkde(x, bw=.8, kernel=function(x)pbeta(.5*x+.5, 3, 3))
if(interactive()){
	plot(ecdf(x))
	curve(pnorm(x), min(x), max(x), add=TRUE, col='red', lwd=2)
	curve(Fhat(x), min(x), max(x), add=TRUE, col='green', lwd=2)
}

fhat=dkde(x)
fhat=dkde(x, bw=.8, kernel='biw')
if(interactive()){
	hist(x, pr=TRUE, br=20)
	curve(dnorm(x), min(x), max(x), add=TRUE, col='red', lwd=2)
	curve(fhat(x), min(x), max(x), add=TRUE, col='green', lwd=2)
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{distribution}
\keyword{smooth}
