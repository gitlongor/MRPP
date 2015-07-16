\name{factoradic}
\Rdversion{1.1}
\alias{dec2FR}
\alias{dec2permvec}
\alias{FR2dec}
\alias{FR2permvec}
\alias{permvec2dec}
\alias{permvec2FR}
\title{Factoradic numbers
}
\description{
Utility functions to convert factoradic numbers among the vector , decimal  and permutation index representation  
}
\usage{
FR2dec(FR)
dec2FR(dec,N)
FR2permvec(FR,base=1L)
permvec2FR(permvec)
dec2permvec(dec,N,base=1L)
permvec2dec(permvec)
}
\arguments{
  \item{FR}{Vector representation of factoradic number. 
}
  \item{dec}{A \code{\link{bigz}} or \code{character} decimal representation of factoradic number. 
}
  \item{N}{An integer indicating the length of vector or permutation vector representation.}
  
  \item{base}{ 0 or 1, indicating the start index for a vector.}
  
  \item{permvec}{An integer vector of permutations of integers of size \code{N}.}
}
\details{
\code{FR2dec} / \code{dec2FR} converts factoradic numbers to/from decimal from/to vector representations. 

\code{FR2permvec} / \code{permvec2FR} converts factoradic numbers from/to vector representations to/from permutation indices. 

\code{dec2permvec} / \code{permvec2dec} are convenient wrappers that convert factoradic numbers from/to decimal representations to/from permutation indices. 
}
\value{
For \code{FR2dec} and \code{permvec2dec}, a \code{bigz} integer; 

for \code{dec2FR} and \code{permvec2FR}, an integer vector;

for \code{FR2permvec} and \code{dec2permvec}, an integer vector of permutations.
}
\references{
James McCaffrey (2012) .NET Test Automation Recipes: A Problem-Solution Approach. Apress [Chapter 10]

http://en.wikipedia.org/wiki/Factoradic

http://www.mathe2.uni-bayreuth.de/frib/KERBER/h00/node30.html
}
\seealso{
\code{\link{permuteTrt}}
}
\examples{
FR2dec(dec2FR(10,10))
set.seed(2340)
FR2permvec(permvec2FR(sample(10)))
permvec2dec(dec2permvec(100,10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{arith}
