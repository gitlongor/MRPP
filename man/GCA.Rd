\name{Gram-Charlier}
\Rdversion{1.1}
\alias{dpearson3gca}
\alias{ppearson3gca}
\alias{dgammagca}
\alias{pgammagca}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
Gram-Charlier approximation up to 4th order moment around a Gamma distribution
}
\description{
\code{dpearson3gca} and \code{ppearson3gca} compute the approximate density and cumulative probability, respectively, by using the Gram-Charlier type A expansion up to the fourth order moment with basis distribution being the Pearson type III distribution that matches the first 3 moments. 

\code{dgammagca} and \code{pgammagca} are similar but the expansions are around a Gamma distribution that matches the first two moments.  
}
\usage{
dpearson3gca(x, mean, sd, skew, exkurt, log=FALSE)
ppearson3gca(q, mean, sd, skew, exkurt, lower.tail=TRUE, log.p=FALSE)

dgammagca(x, mean, sd, skew, exkurt, log=FALSE)
pgammagca(q, mean, sd, skew, exkurt, lower.tail=TRUE, log.p=FALSE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{x,q,  mean,sd,
  log,lower.tail,log.p}{As in \code{\link{dnorm}}.
}
  \item{skew}{Skewness coefficient, i.e., standardized 3rd cumulant.
}
  \item{exkurt}{Excess kurtosis coefficient, i.e., standardized 4th cumulant.
}
}
\note{
\code{dpearson3gca} and \code{ppearson3gca}, an affine translation is performed such that the first three moments matches a Gamma distribution, whereas the fourth moment is adjusted polynomially through Gram-Charlier expansion around this Gamma distribution. 

\code{dgammagca} and \code{pgammagca}, a simple scale transformation (without location shift) is performed such that the first two moments matches a Gamma distribution, whereas the third and the fourth moments are adjusted polynomially through Gram-Charlier expansion around this Gamma distribution. 

Negative density values are truncated as zero. 
}
\value{
A numeric vector of (possibly log-transformed) evaluated densities or cumulative probabilities.
}
\references{
Page 25--26 of Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.

Newton L. Bowers, Jr. (1966) Expansion of Probability Density Functions as a Sum of Gamma Densities with Applications in Risk Theory (with Discussions). \emph{Transactions of Society of Actuaries}, \bold{18}, pt. 1, 125--147.
}
\author{Long Qu
}
\seealso{
\code{\link{dpearson3}}, \code{\link[PDQutils]{dapx_gca}}
}
\examples{
dpearson3gca(1,2,3,4,5)
ppearson3gca(4,3,2,1,0)

dgammagca(1,2,3,4,5)
pgammagca(4,3,2,1,0)
if(interactive()){
	curve(dpearson3(x,1,2,3,log=FALSE), -1, 10)
	curve(dpearson3gca(x,1,2,3,4, log=FALSE), -1, 10, add=TRUE, col=4)
	curve(dgammagca(x,1,2,3,4, log=FALSE), -1, 10, add=TRUE, col=2)
	
	kurtosis(rpearson3(1e6,1,2,3), excess=TRUE) ## 1.5*3^2
	curve(dpearson3(x,1,2,3,log=FALSE), -1, 10)
	curve(dpearson3gca(x,1,2,3,1.5*3^2, log=FALSE), -1, 10, add=TRUE, col=4, lwd=3,lty=3)
	curve(dgammagca(x,1,2,3,1.5*3^2, log=FALSE), -1, 10, add=TRUE, col=2)
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{distribution}
