\name{Transformed Log-Normal}
\Rdversion{1.1}
\alias{tlnorm}
\alias{dtlnorm}
\alias{ptlnorm}
\alias{qtlnorm}
\alias{rtlnorm}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
Affine Transformed Log-Normal Distribution
}
\description{
\code{dtlnorm}, \code{ptlnorm}, \code{qtlnorm} and \code{rtlnorm} return the density, cumulative probability, quantile, and pseudo-random numbers, respectively, from an affine transformed Log-Normal distribution with given mean, standard deviation, and skewness coefficients.

}
\usage{
dtlnorm(x, mean, sd, skew, log=FALSE)
ptlnorm(q, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
qtlnorm(p, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
rtlnorm(n, mean, sd, skew)

}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{x,p,q,n,  mean,sd,
  log,lower.tail,log.p}{As in \code{\link{dnorm}}.
}
  \item{skew}{Skewness coefficient, i.e., standardized 3rd cumulant.
}
}
\note{
When \code{skew=0}, results from the usual normal distribution with specified \code{mean} and \code{sd} will be returned. Other wise an Log-Normal distribution with given \code{skew} will be chosen, and will be affine transformed to match its mean and sd with given values.


}
\value{
A numeric vector of (possibly log-transformed) evaluated densities, cumulative probabilities, quantiles, or pseudo-random numbers.
}
%\references{
%Page 25--26 of Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.
%}
\author{Long Qu
}
\seealso{
\code{\link{dnorm}}, \code{\link{dpearson3}}
}
\examples{
dtlnorm(1,2,3,4)
ptlnorm(4,3,2,1)
qtlnorm(ptlnorm(4,3,2,1),3,2,1)
rtlnorm(4,3,2,1)

if(interactive()){
	hist(rtlnorm(1e5L,1,2,3), pr=TRUE, br=50)
	curve(dtlnorm(x,1,2,3), -1, 20, add=TRUE)
	
	qqplot(
		rtlnorm(1e4L, 2,2,2), 
		qtlnorm(runif(1e4),2,2,2)
	)
	abline(0,1)
	
	plot(ecdf(rtlnorm(1e4L,3,2,1)))
	curve(ptlnorm(x,3,2,1), -3,13,add=TRUE,col=2,lwd=3,lty=3)
	
	curve(dtlnorm(x,1,2,3,log=TRUE), -1, 20)
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{distribution}
