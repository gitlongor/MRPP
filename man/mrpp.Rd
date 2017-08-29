\name{MRPP object}
\Rdversion{1.1}
\alias{mrpp}
\alias{mrpp.matrix}
\alias{mrpp.dist}
\alias{mrpp.formula}
\alias{$.mrpp}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
MRPP Object
}
\description{
The \code{"mrpp"} object is a setup that allows the multiresponse permutation procedure (MRPP) method. 
}
\usage{
mrpp(y, B, permutedTrt, weight.trt='df', distFunc=dist,...)

\method{mrpp}{matrix}(y, trt, B, permutedTrt, weight.trt='df', distFunc=dist,...)
\method{mrpp}{dist}(y, trt, B, permutedTrt, weight.trt='df', distFunc=dist, R, data, ...)
\method{mrpp}{formula}(y,  B, permutedTrt, weight.trt='df', distFunc=dist, ...)

\method{$}{mrpp}(x,name)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{y}{For \code{mrpp.matrix}, \code{y} is a numeric data matrix, with rows being samples and columns being response variables. For \code{mrpp.dist}, this is a \code{\link{dist}} object containing the pairwise inter-sample distances. For \code{mrpp.formula}, this is a two-sided formula. 
}
  \item{trt}{A vector of treatment assignments. If missing, the first permutation in \code{permutedTrt} will be used as the original treatment assignment. As least one of \code{trt} or \code{permutedTrt} need to be given. 
}
  \item{B}{A positive integer number of permutations requested. This will only be used when \code{permutedTrt} is missing. This will be overwritten by the \code{ncol(permutedTrt)}.
}
  \item{permutedTrt}{An optional \code{permutedTrt} object specifying the permutations of treatment assignments. If missing, this will be computed as \code{permuteTrt(trt,B)}. See \code{\link{permuteTrt}}. 
}
  \item{weight.trt}{Either a character scalar being one of \code{c("df", "n", "equal", "pairs")} or a positive numeric vector of length \code{nlevels(trt)} subject to sum-to-1 constraint. See details. 
}
  \item{distFunc}{A function that computes a \code{\link{dist}} object given data \code{y}.
}
%  \item{idxOnly}{As in \code{\link{permuteTrt}}.
%}
  \item{R}{A positive integer, equal to the number of response variables. When missing, it is determined from \code{data} (when available) or computed as \code{attr(y,'Size')-1L}.
}
 \item{data}{A numeric data matrix for \code{mrpp.matrix}, with rows being samples and columns being response variables. 
}
  \item{\dots}{For \code{mrpp.formula}, these are the additional arguments to be passed on to \code{\link{model.frame.default}} 
}
  \item{x}{An \code{mrpp} object.
}
  \item{name}{A symbol, specifying the name of the element to be extracted. 
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
\item{R }{A positive integer of the number of dimensions of the response variable. 
}
\item{B.requested }{The number of permutations requested. This may not be the actual number of permutations being used, e.g., when it is larger than \code{nparts}. 
}
\item{trt }{A factor of original treatment assignments. This is identical to the \code{trt} argument when it is non-missing; otherwise, this is reconstructed from the \code{permutedTrt} argument. 
}
\item{permutedTrt.env}{An environment that contains two objects, \code{"B"} and \code{"permutedTrt"}. If the \code{permutedTrt} argument is missing, then both objects \code{permutedTrt.env} environment are assigned as \emph{promises} initially. There is no guarantee that they will be eventually evaluated, e.g., when moment-matching approximations are used to obtain approximate p-values in \code{mrpp.test}. 
}
%\item{idxOnly}{A logical scalar, as in \code{\link{permuteTrt}}. 
%}
\item{weight.trt}{A numeric vector of the treatment group weights, with the \code{"method"} attribute being the name of the treatment group weighting method as the \code{weight.trt} argument. 
}
\item{distFun}{The same as the \code{distFunc} argument. 
}
\item{data.env}{An environment containing one object \code{"y"}, being the raw data set. This is created as a promise. Some methods (e.g., \code{\link{mrpp.test}}) do not require this and it will not be evaluated. For \code{mrpp.dist}, if \code{data} is not given, then the result is only a \emph{possible} data set reconstructed by \code{\link{cmdscale}}.
}

\code{x$name} returns the element named \code{
"name"} of an \code{"mrpp"} object \code{x}, except that \code{x$permutedTrt} returns \code{x[['permutedTrt.env']][['permuteTrt']]} and that \code{x$y} returns \code{x$data.env$y}.

}
\note{
Currently, the following S3 methods are available for the \code{"mrpp"} class: 
\code{\link[MRPP:cumulant]{mean, sd, var, skewness, kurtosis,  moment, cumulant}} and \code{\link{mrpp.test}}.
}
\section{Warning}{
The \code{permutedTrt.env} component and \code{data.env} are environemnts such that the duplication of a \code{"mrpp"} object will share the same environments. 
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
identical(x, mrpp.obj$y) # mrpp.obj$y extracts mrpp.obj$data.env$y

cumulant(mrpp.obj)
mrpp.test(mrpp.obj, method='pearson3gca') # fast approximation w/o permutation
str(mrpp.obj$permutedTrt) #now takes time to force the promise
ls(mrpp.obj$permutedTrt.env)
mrpp.test(mrpp.obj, method='permutation') 

dist.x=dist(x)
mrpp.obj2 = mrpp(dist.x, trt, B=Inf)
cumulant(mrpp.obj2) # same as above
mrpp.test(mrpp.obj2) # same as above
mrpp.obj2$y # recovered data, typically not the same as x
max(abs(dist(mrpp.obj2$y)-dist.x)) # max error of recovery of distances

mrpp.obj3 = mrpp(x~trt, N=Inf) # equivalent formula interface
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}