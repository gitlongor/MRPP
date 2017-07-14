\name{sample1}
\Rdversion{1.1}
\alias{sample1}
\alias{rejsample}
\alias{HSEL.bigz}
\alias{sample.bigz}
\title{
Sample without replacement when the population can be huge. 
}
\description{General a small sample of size M out of a huge population of size N \emph{without} replacement. 
}
\usage{
HSEL.bigz(N, M, seed=0)
rejsample(N, M, seed=0)
sample1(N, seed=0)
sample.bigz(N, M)
}
\arguments{
  \item{N}{A \code{bigz} integer of the size of the population.
}
  \item{M}{The sample size. This needs to be a usual integer (no larger than \code{.Machine$integer.max}). 
}
  \item{seed}{A \code{bigz} integer of random number seed. See \code{\link{urand.bigz}}. 
}
}
\details{
\code{HSEL.bigz} generates a sample of size M from a population of size N, with space complexity linear in M, following the algorithm of Jarmo Ernvall and Olli Nevalainen (1982). For fixed M, the speed of \code{HSEL.bigz} does not deteriorate much as N grows. 

\code{rejsample} is similar to \code{HSEL.bigz}, but implemented using direct rejection through \code{\link{unique}}. 

\code{sample1} generates one number of 0 to N-1 using simple rejection and \code{\link{urand.bigz}}. 

\code{sample.bigz} is similar to \code{HSEL.bigz}, but using an approximate implementation as follows. This function breaks all \code{N} integers into chunks of size \code{.Machine$integer.max}, except for the last chunk which might be smaller. Then a multinomial vector of sample sizes is simulated for these chunks. Within each chunk, the base version of \code{\link{sample.int}} is called to simulate the required number of integers without replacement and the results are shifted by appropriate multiples of \code{.Machine$integer.max}. This is approximate in the sense that multinomial probabilities being used are only double precision float points rather than exact \code{bigq} objects.  If the difference between as.bigq(.Machine$integer.max, N) and as.numeric(as.bigq(.Machine$integer.max, N)) is sufficiently small, the approximateness can be ignored for practical purposes. Empirially, \code{sample.bigz} is typically faster than \code{HSEL.bigz} if \code{N} is no larger than \code{6e15}. 

}
\note{
\code{HSEL.bigz}, \code{rejsample} and \code{sample1} are exact algorithms and rely on \code{gmp}'s pseudo-random number generator. But \code{sample.bigz} only uses \code{R}'s default pseudo-random number generator. 
}
\section{warning}{
\code{rejsample} can be very slow when M is large. 

\code{sample.bigz} can consume a large amount of memory if N is large. 
}
\value{
For \code{HSEL.bigz}, \code{rejsample}, and \code{sample.bigz}, a \code{bigz} vector of length \code{N}, between 1 and N. 
For \code{sample1}, a \code{bigz} integer, between 0 and N-1. 
}
\references{
Jarmo Ernvall and Olli Nevalainen. 1982. An Algorithm for Unbiased Random Sampling. THE COMPUTER JOURNAL, VOL. 25, NO. 1, 45--47
}
\seealso{
\code{\link{sample}}, \code{\link{urand.bigz}}, \code{\link{RNGkind}}
}
\examples{
HSEL.bigz(factorialZ(15),5,12345)
rejsample(factorialZ(15),5,12345)
sample1(factorialZ(15), 12345)

set.seed(12345)
sample.bigz(factorialZ(15),5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{iteration}
\keyword{arith}