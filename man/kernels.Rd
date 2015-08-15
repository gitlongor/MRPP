\name{kernel}
\Rdversion{1.1}
\alias{pkernel}
\alias{dkernel}
\alias{fourier.kernel}
\alias{pkde}
\alias{dkde}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
CDF, PDF and Fourier Transform of Univariate Kernel Densities
}
\description{
\code{pkernel}, \code{dkernel} and \code{fourier.kernel} return the cumulative distribution function, probability density function and Fourier transform of the density, respectively, for common univariate kernels specified by names. 

\code{pkde} and \code{dkde} return the \emph{fitted} cumulative distribution function and probability density function. 
}
\usage{
pkernel(kernel= c("gaussian", "biweight", 'triweight', 'tricube',"logistic"))
dkernel(kernel= c("gaussian", "biweight", 'triweight', 'tricube',"logistic"))
fourier.kernel(kernel= c("gaussian", "biweight", 'triweight', 'tricube',"logistic"), 
	root.2pi=TRUE)

pkde(x, bw=bw.nrd, kernel=c("gaussian", "biweight", 'triweight', 'tricube',"logistic"))
dkde(x, bw=bw.nrd, kernel=c("gaussian", "biweight", 'triweight', 'tricube',"logistic"))
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{kernel}{A character scalar of the name of the kernel. Partial string matching is allowed. For \code{pkde}, this can also be a valid univariate cumulative distribution function. 
}
  \item{root.2pi}{A logical scalar. When \code{TRUE}, a factor of \code{1/sqrt(2*base::pi)} will be multiplied to the integral transform; otherwise, the integral transform itself is computed. 
}
  \item{x}{A numeric vector of data values.
}
  \item{bw}{A positive integer indicating the bandwidth to be used in kernel density estimation. Alternatively, this can be a function that computes the bandwidth, or the character name of such a function. 
}
}
\note{
Currently, kernels whose derivative of the density function is continuous are implemented. 

High-order kernels that are not bona fide distributions or non-smooth kernels are not allowed. 

The \code{pkde} implementation does \emph{not} make use of the fast Fourier transform (FFT); but the \code{dkde} uses FFT. 
}
\value{
A function (closure) that accepts a numeric vector as the argument. 
}
\references{
\url{https://en.wikipedia.org/wiki/Kernel_(statistics)}

Section 3.5 of Bernard W. Silverman (198) Density Estimation for Statistics and Data Analysis. London; Newyork: Chapman & Hall. 
}
\author{Long Qu
}
\seealso{
\code{\link{density}}
}
\examples{
f=dkernel('biweight')
F=pkernel('biweight')
f(0) ## 15/16 = 0.9375
F(0) ## 0.5
f(-2:2)
F(-2:2)
ft=fourier.kernel('biweight')
ft(-2:2)
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
