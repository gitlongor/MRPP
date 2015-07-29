\name{Extract.dist}
\alias{[[.dist}
\title{
Extracting Distances from a dist Object
}
\description{
For an object \code{x} of class \code{"dist"}, \code{x[[i,j]]} is equivalent to \code{as.matrix(x)[[i,j]]}, but without the intermediate \code{as.matrix} stage. 
}
\usage{
\method{`[[`}{dist}(x, i, j)
}
\arguments{
	\item{x}{An object of class \code{"dist"} (see \code{\link{dist}}). 
}
	\item{i,j}{Positive integer scalars. 
}
}

\value{
A single numeric scalar that is equal to \code{as.matrix(x)[[i,j]]}.
}
\author{
Long Qu
}

\seealso{
\code{\link{dist}}
}

\examples{
set.seed(2340)
x=matrix(rnorm(10*5),10)
d=dist(x)
dmat=as.matrix(d)
d[[1,1]]  ## dmat[[1,1]]
d[[5,3]]  ## dmat[[5,3]]
d[[2,4]]  ## dmat[[2,4]]
}

\keyword{ array }

