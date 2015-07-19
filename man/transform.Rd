\name{transform}
\alias{transform}
\alias{transform.matrix}
\alias{transform.numeric}
\alias{binom.symm}
\alias{binom.quadLL}
\alias{pois.vst}
\alias{pois.symm}
\alias{pois.quadLL}
\alias{nbinom.vst}
\alias{nbinom.symm}
\alias{nbinom.quadLL}
\alias{fisher.z}
\alias{chisq.symm}
\alias{binom.vst}
\alias{norm.vst}
\alias{norm.symm}
\alias{norm.quadLL}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
Transformation of Matrices or Numeric Vectors
}
\description{
The \code{transform} method for \code{matrix} and \code{numeric} vectors is a wrapper that perform common transformations for common exponential family distributions. 
}
\usage{
\method{transform}{matrix}(`_data`, method = "I", \dots)
\method{transform}{numeric}(`_data`, \dots)

binom.symm(y, ...)
binom.quadLL(y, ...)
pois.vst(y)
pois.symm(y)
pois.quadLL(y)
nbinom.vst(y, size, tau = 1/size)
nbinom.symm(y, size, tau = 1/size)
nbinom.quadLL(y, size, tau = 1/size, standardize = FALSE)
chisq.symm(y)
binom.vst(y)
norm.vst(x)
norm.symm(x)
norm.quadLL(x)
fisher.z(x)

}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{_data, x, y}{
	A numeric vector or matrix of data to be transformed. 
}
  \item{method}{
	A character name of the function that performs the actual transformation. The default is \code{base::I} function, i.e., no transformation at all. See details. 
}
  \item{size}{
	A positive numeric (not necessarily integer) scalar or vector as in \code{\link{dnbinom}}. In case of a vector when \code{_data} (\code{y}) is a matrix, then if  \code{length(size}==ncol(y) (or \code{ncol(`_data`)}, it is assumed that the \code{size} parameter is applicable to each whole column.  
}
  \item{tau}{
	A reparameterization that should always be \code{1/size}.
}
  \item{standardize}{
	Logical; if \code{TRUE}, the quadratic log-likelihood transformation is further multiplied by a factor depending on \code{y}. See details. 
}
  \item{\dots}{
	Optional arguments to be passed to the function determined by the \code{method} argument. 
}
}

\details{
\code{method} argument can be any (full) function name. 
If \code{method} is any one of 
\preformatted{c("nbinom.vst", "nbinom.quadLL", "nbinom.symm", 
"pois.vst", "pois.quadLL", "pois.symm", 
"binom.vst", "binom.quadLL", "binom.symm", 
"norm.vst", "norm.quadLL", "norm.symm", 
"chisq.vst", "chisq.quadLL", "chisq.symm", 
"fisher.z"),} 
partial string matching is allowed, if no other function in the search path has a name that matches \code{method} perfectly. Most of the above function names have the form of \dQuote{distribution.transform} where the \dQuote{distribution} part is typically the abbreviated distribution name as those used in the suffices in functions like \code{\link{rnorm}}. The \dQuote{transform} part is one of \code{"vst", "quadLL"} and \code{"symm"}, indicating the asymptotic variance stabilizing transformation, asymptotic quadratic log-likelihood transformation, and asymptotic zero skewness transformation. 

Only \code{transform.matrix} and \code{transform.numeric} are exported. 

FIXME: Give details of each transformation here... 
}
\value{
A numeric vector or a numeric matrix of transformed data. 
}
\references{
Francis J. Anscombe (1948) The transformation of Poisson, Binomial and Negative-binomial Data. \emph{Biometrika}. \bold{35}, 246--254.
Francis J. Anscombe (1964) Normal Likelihood Functions. \emph{Ann Inst Stat Math}. \bold{16}, 1--19.
Philip Hougaard (1982) Parametrization of Non-linear Models. \emph{J R Statist Soc B}. \bold{44}(2), 244--252.
Edwin B. Wilson and Margaret M. Hilferty (1931) The Distribution of Chi-Square. \emph{Proc Nat Acad Sci USA}. \bold{17}, 684--688.
}
\author{
Long Qu
}

\seealso{
\code{\link[base:transform]{base::transform}} for transformation of variables in \code{data.frame}. 
}
\examples{
y=rnbinom(1e3, size=10, prob=.5)  ## simulate example data
identical(y, transform(y))
y.log1p=transform(y, method="log1p")      ## base::log1p
y.vst=transform(y, method="nbinom.vst", tau=0.1)  ## towards variance stabilization
y.quadll=transform(y, method="nbinom.quadLL", tau=0.1)  ## towards quadratic log-likelihood
y.symm=transform(y, method="nbinom.symm", tau=0.1)  ## towards symmetry

if(interactive()){ ## visually assessing normality thru QQ plot
	qqnorm(y ); 	qqline(y)
	qqnorm(y.log1p);qqline(y.log1p);
	qqnorm(y.vst);qqline(y.vst);
	qqnorm(y.quadll);qqline(y.quadll);
	qqnorm(y.symm);qqline(y.symm);
}
}
\keyword{ manip }
\keyword{ distribution }% __ONLY ONE__ keyword per line
