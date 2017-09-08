\name{PearsonIII}
\Rdversion{1.1}
\alias{pearson3}
\alias{pearsonIII}
\alias{dpearson3}
\alias{ppearson3}
\alias{qpearson3}
\alias{rpearson3}
\alias{Rdpearson3}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
Pearson type III distribution
}
\description{
\code{dpearson3}, \code{ppearson3}, \code{qpearson3} and \code{rpearson3} return the density, cumulative probability, quantile, and pseudo-random numbers, respectively, from a Pearson type III distribution with given mean, standard deviation, and skewness coefficients.

\code{Rdpearson3} returns the roughness of the density (derivative) of the Pearson type III distribution with given mean, standard deviation and skewness coefficients.

}
\usage{
dpearson3(x, mean, sd, skew, log=FALSE)
ppearson3(q, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
qpearson3(p, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
rpearson3(n, mean, sd, skew)
Rdpearson3(sd, skew, deriv.order=0, power=1)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{x,p,q,n,  mean,sd,
  log,lower.tail,log.p}{As in \code{\link{dnorm}}.
}
  \item{skew}{Skewness coefficient, i.e., standardized 3rd cumulant.
}
  \item{deriv.order}{A vector of non-negative integers. This is the order of derivatives of the density, with 0 being the density itself.
}
  \item{power}{A numeric scalar. The integral is raised to this power. 
}
}
\note{
Pearson type III distribution implemented here is simply an affine transformed Gamma distribution when \code{skew!=0}, and the normal distribution when \code{skew=0}.


}
\value{
\code{dpearson3}, \code{ppearson3}, \code{qpearson3} and \code{rpearson3} return a numeric vector of (possibly log-transformed) evaluated densities, cumulative probabilities, quantiles, or pseudo-random numbers. 

\code{Rdpearson3} returns a numeric vector, equal to 
 \deqn{\left[\int_{-\infy}^{\infy} \lbrace f^{(r)}(x)\rbrace^2 dx \right]^p}{[integral of {f^(r)}^2 ]^p}

where \emph{r} is \code{deriv.order} and \emph{p} is \code{power}. When \code{deriv.order} is not smaller than \code{4/skew-0.5}, the result will be \code{Inf}.
}

\references{
Page 25--26 of Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.
}
\author{Long Qu
}
\seealso{
\code{\link{dnorm}}, \code{\link{tlnorm}}
}
\examples{
dpearson3(1,2,3,4)
ppearson3(4,3,2,1)
qpearson3(ppearson3(4,3,2,1),3,2,1)
rpearson3(4,3,2,1)
Rdpearson3(.5,1.1,0:2,1)

if(interactive()){
	hist(rpearson3(1e5L,1,2,3), pr=TRUE, br=50)
	curve(dpearson3(x,1,2,3), -1, 20, add=TRUE)
	
	qqplot(
		rpearson3(1e4L, 2,2,2), 
		qpearson3(runif(1e4),2,2,2)
	)
	abline(0,1)
	
	plot(ecdf(rpearson3(1e4L,3,2,1)))
	curve(ppearson3(x,3,2,1), -3,13,add=TRUE,col=2,lwd=3,lty=3)
	
	curve(dpearson3(x,1,2,3,log=TRUE), -1, 20)
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{distribution}
