\name{MRPP object}
\Rdversion{1.1}
\alias{mrpp}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
MRPP Object
}
\description{
The \code{"mrpp"} object is a setup that allows the multiresponse permutation procedure (MRPP) method. 
}
\usage{
mrpp(y, trt, B=as.integer(min(nparts(table(trt)), 1e4L)),
	permutedTrt, weight.trt="df", distFunc=dist, idxOnly=FALSE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{y}{A numeric data matrix, with rows being samples and columns being response variables. 
}
  \item{trt}{A vector of treatment assignments. If missing, the first permutation in \code{permutedTrt} will be used as the original treatment assignment. As least one of \code{trt} or \code{permutedTrt} need to be given. 
}
  \item{B}{A positive integer number of permutations requested. This will only be used when \code{permutedTrt} is missing. This will be overwritten by the \code{ncol(permutedTrt)}.
}
  \item{permutedTrt}{An optional permutation index matrix. If missing, this will be computed as \code{permuteTrt(trt,B)}. See \code{\link{permuteTrt}}. 
}
  \item{weight.trt}{Either a character scalar being one of \code{c("df", "n", "equal", "pairs")} or a positive numeric vector of length \code{nlevels(trt)} subject to sum-to-1 constraint. See details. 
}
  \item{distFunc}{A function that computes a \code{\link{dist}} object given data \code{y}.
}
  \item{idxOnly}{As in \code{\link{permuteTrt}}.
}
}
\details{
For the \code{weight.trt} argument, suppose that \code{n} denote the vector sample sizes for each treatment. Then, 
\itemize{
\item{\code{weight.trt="df"} uses group (treatment) weights proportional to \code{n-1};}
\item{ \code{weight.trt="n"} uses group (treatment) weights proportional to \code{n};}
\item{\code{weight.trt="equal"} uses equal group (treatment) weights; and }
\item{\code{weight.trt="pairs"} uses group (treatment) weights proportional to \code{choose(n,2)}.}
} 
}

\value{
\code{mrpp} returns an object of class \code{"mrpp"}, which is a list with the following components:

\item{distObj }{The \code{"dist"} object of data.
}
\item{n }{An integer vector of sample sizes in each treatment. 
}
\item{ntrt }{A positive integer of the number of treatment groups. 
}
\item{nparts }{A \code{bigz} integer of the total number of non-equivalent permutations. 
}
\item{nobs }{A positive integer of the total sample size. 
}
\item{R }{A positive integer of the number of dimensions of the response variable \code{y}. 
}
\item{B.requested }{The number of permutations requested. This may not be the actual number of permutations being used, e.g., when it is larger than \code{nparts}. 
}
\item{trt }{A factor of original treatment assignments. This is identical to the \code{trt} argument when it is non-missing; otherwise, this is reconstructed from the \code{permutedTrt} argument. 
}
\item{permutedTrt.env}{An environment that contains two objects, \code{"B"} and \code{"permutedTrt"}. If the \code{permutedTrt} argument is missing, then both objects \code{permutedTrt.env} environment are assigned as \emph{promises} initially. There is no guarantee that they will be eventually evaluated, e.g., when moment-matching approximations are used to obtain approximate p-values in \code{mrpp.test}. 
}
\item{idxOnly}{A logical scalar, as in \code{\link{permuteTrt}}. 
}
\item{weight.trt}{A numeric vector of the treatment group weights, with the \code{"method"} attribute being the name of the treatment group weighting method as the \code{weight.trt} argument. 
}
\item{distFun}{The same as the \code{distFunc} argument. 
}
}
\note{
Currently, the following S3 methods are available for the \code{"mrpp"} class: 
\code{\link[MRPP:cumulant]{mean, sd, var, skewness, kurtosis,  moment, cumulant}} and \code{\link{mrpp.test}},  
}
\section{Warning}{
The \code{permutedTrt.env} component is an environemnt such that the duplication of a \code{"mrpp"} object will share the same environment. 
}
\references{
Paul W. Mielke, Kenneth J. Berry. (2007) Permutation Methods: A Distance Function Approach. 2nd ed. Springer.
}
\author{
Long Qu
}
%%\note{
%%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{permuteTrt}}, \code{\link{moment}}, \code{\link{mrpp.test}}
}
\examples{
set.seed(2340)
n=10
trt=gl(2,n)
x=matrix(rnorm(n*5*2), ncol=5)
mrpp.obj = mrpp(x, trt, B=Inf)
cumulant(mrpp.obj)
mrpp.test(mrpp.obj)
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}