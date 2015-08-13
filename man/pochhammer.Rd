\name{Rising Factorial}
\alias{factorial.rising}
\title{
Rising Factorial (Pochhammer Symbol)
}
\description{
For non-negative integer \code{n} and a real \code{x}, \code{factorial.rising(x,n)} computes \deqn{x(x+1) \cdots x(x+n-1)}{x(x+1) \dots (x+n-1)}.

For negative integer \code{n} and a real \code{x}, \code{factorial.rising(x,n)} computes \deqn{(x+n)(x+n+1) \cdots (x-1)}{(x+n)(x+n+1) \dots (x-1)}.

}
\usage{
factorial.rising(start, nterms)
}
\arguments{
	\item{start}{A vector of starting values, i.e. the smallest factor in the rising factorial. These can be integer, double, bigz, or bigq vectors. 
}
	\item{nterms}{A vector of non-negative integers, the number of terms in the rising factorial. Negative values and non-integer values are not yet supported. 
}
}

\value{
When \code{nterms>=0}, the result is a vector with class being the same as that of \code{start}, except when \code{start} is an integer but the result is too large to be represented as an integer; in the latter case the result will be a \code{'bigz'} vector.  

When \code{nterms<0}, the result is either a double vector or a \code{'bigq'} vector, depending on whether \code{start} is numeric or \code{'big[zq]'}.
}
\author{
Long Qu
}

\seealso{
\code{\link{factorial}}, \code{\link{factorialZ}}, \code{\link{prod.bigz}}, \code{\link[orthopolynom]{pochhammer}}
}

\examples{
factorial.rising(1,5)  ## should be equal to 
factorial(5)

factorial.rising(10L,5)  ## integer
factorial.rising(10,5)   ## double
factorial.rising(as.bigz(10),5)  ## bigz
factorial.rising(as.bigq(10),5)  ## bigq

factorial.rising(10.5,5)  ## double
factorial.rising(as.bigq(10.5),5)  ## bigq

factorial.rising(10, -5)  ## double
factorial.rising(10L, -5)  ## double
factorial.rising(as.bigz(10),-5)  ## bigq
factorial.rising(as.bigq(10),-5)  ## bigq

\dontshow{
stopifnot(identical(
	factorial.rising(1L,5),
	as.integer(factorial(5))
))
}

factorial.rising(1L,50)  ## shoulde be identical to 
factorialZ(50)			## or
prod.bigz(seq(from=1, by=1, length=50))
\dontshow{
stopifnot(identical(
	factorial.rising(1L,50), 
	factorialZ(50)
))
stopifnot(identical(
	factorial.rising(1L,50), 
	prod.bigz(seq(from=1, by=1, length=50))
))
}

factorial.rising(10,5)
factorial.rising(10,0)
factorial.rising(10,0:5)

factorial.rising(1:4, 10)
factorial.rising(1:5, 10)
factorial.rising(5:8, 10:13)
factorial.rising(5:8, -3:0)
factorial.rising(5:8+0, -3:0)
factorial.rising(as.bigz(5:8), -3:0)

}

\keyword{ math }

