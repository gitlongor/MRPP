\name{Rising Factorial}
\alias{factorial.rising}
\title{
Rising Factorial (Pochhammer Symbol)
}
\description{
For positive integer \code{x} and {n}, \code{factorial.rising(x,n)} computes \deqn{x(x+1) \cdots x(x+n-1)}{x(x+1) \dots (x+n-1)}.
}
\usage{
factorial.rising(start, nterms)
}
\arguments{
	\item{start}{A vector of non-negative integers, the smallest factor in the rising factorial.
}
	\item{nterms}{A vector of non-negative integers, the number of terms in the rising factorial. 
}
}

\value{
The result is a vector of class \code{integer} when the maximum of the result is representable in the integer mode; otherwise the result will be a vector of class \code{bigz}. 
}
\author{
Long Qu
}

\seealso{
\code{\link{factorial}}, \code{\link{factorialZ}}, \code{\link{prod.bigz}}
}

\examples{
factorial.rising(1,5)  ## should be equal to 
factorial(5)
\dontshow{
stopifnot(identical(
	factorial.rising(1,5),
	as.integer(factorial(5))
))
}

factorial.rising(1,50)  ## shoulde be identical to 
factorialZ(50)			## or
prod.bigz(seq(from=1, by=1, length=50))
\dontshow{
stopifnot(identical(
	factorial.rising(1,50), 
	factorialZ(50)
))
stopifnot(identical(
	factorial.rising(1,50), 
	prod.bigz(seq(from=1, by=1, length=50))
))
}

factorial.rising(10,5)
factorial.rising(10,0)
factorial.rising(10,0:5)

factorial.rising(1:4, 10)
factorial.rising(1:5, 10)
factorial.rising(5:8, 10:13)

}

\keyword{ math }

