\name{MRPP importance}
\Rdversion{1.1}
\alias{grad.kdep}
\alias{p.value.grad.kdep}
\alias{hessian.smoothp}
\alias{get.p.dd.dw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation-based variable importance measure
}
\description{
\code{grad.kdep} computes variable importance measure of partial derivative of continuous approximation to MRPP p-values with respect to hypothetical weights. 

\code{get.p.dd.dw} computes variable importance measure by using partial derivative of distance function with respect to hypothetical weights as the new distance function in MRPP

\code{p.value} method for class \code{grad.kdep} computes the approximate mid-p-values from smoothing. 
}
\usage{
grad.kdep(y,  bw, kernel='triweight', adjust=NULL,
	mrpp.stats=NULL, r=seq_len(y$R), test=FALSE)
get.p.dd.dw(y, permutedTrt, r = seq_len(ncol(y)), weight.trt="df", eps=1e-8)
\method{p.value}{grad.kdep}(x, type=c("keep1","drop1","add1"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ An \code{\link{mrpp}} object.
}
  %\item{permutedTrt}{A list matricies that list all random treatment assignments to be used in the MRPP test. 
%}
  \item{r}{A positive integer vector, specifying variables for which variable importance measure is requested. 
}
  \item{test}{\code{TRUE} or \code{FALSE}. If \code{FALSE}, the first permutation in \code{y$permutedTrt} will be used as the original treatment assignment. If \code{TRUE}, each permutation in \code{y$permutedTrt} will be treated as if it is the original treatment assignment, useful for doing permutation tests using the computed variable importance measure as the test statistic.
}
  %\item{distObj}{The numeric vector computed from \code{dist(y)}, i.e., the lower triangle of the Euclidean distance matrix. 
%}
  \item{mrpp.stats}{An optional vector of all observed MRPP test statistics. If this is \code{NULL}, then 
  
 \code{mrpp.test(y,test.method='permutation')$all.statistics}

  will be called to obtain it.  
}
  \item{kernel}{A character scalar taking values in \code{\link{.smooth1.kernels}}.  
}
  \item{bw}{A numeric or character scalar. If being numeric, it gives the bandwidth used for kernel density estimation of \code{mrpp.stats}. This might be positive infinity. If \code{bw} is a character, it takes values as the \code{bw.method} argument of \code{\link{bw.smoothp}}. 
}
%  \item{cperm.mat}{}
  \item{weight.trt}{As in \code{\link{mrpp.test.dist}}.
}
  \item{adjust}{Either \code{NULL} or a character scalar being one of \code{c('none','weighted.mean','scale','log scale')}. See details. 
}
  \item{eps}{A small positive tolerance, the same as described in \code{\link{mrpp.test}}. }
  \item{x}{An \code{"grad.kdep"} object returned from calling \code{grad.smooothp}.
}
  \item{type}{A character scalar, being one of \code{c("keep1","drop1","add1")}, specifying the type of approximated mid-p-values using gradients. 
}
  \item{\dots}{Not used. 
}
}

\details{
\code{grad.kdep} first uses kernel density estimator to approximate the distribution of \code{mrpp.stats} and obtain a continuous approximate p-value for the MRPP test. Then it computes the partial derivative of the approximate p-value with respect to the hypothetical weight of each variable, evaluated at weight 1. 
That is, if the partial derivative is negatively large, then the corresponding variable will be more important, because increasing its weight will quickly decreasing the approximate MRPP p-value. Vice versa. 

When \code{adjust='none'}, the partial derivatives of p-value with respect to weights will be returned. When \code{adjust='weighted.mean'}, the partial derivatives will be treated as a weighted average, i.e., will be divided by total weights. When \code{adjust='scale'}, they are further divided standard deviations. When \code{adjust='log scale'}, the partial derivatives are first computed on the log p-value scale and the results are then exponentiated. When \code{adjust=NULL}, then as long as \code{bw} being used is finite, \code{adjust='none'} will be the default setting; otherwise, \code{adjust='weighted.mean'} will be the default setting. 

\code{get.p.dd.dw} first computes the partial derivative of the Euclidian distance measure with respect to the hypothetical weights for each variable. This is then used as the distance function for the variable to be used in MRPP. The resulting MRPP p-value is returned as the variable importance. The smaller, the more important.
}

\value{
\code{grad.kdep} returns an object of class \code{"grad.kdep"}. This is an \code{length(r)}-vector of importance measure when \code{test} is \code{FALSE}.  When \code{test} is \code{TRUE}, this will be an \code{B}-by-\code{length(r)} matrix of importance measure, where \code{B} is the number of permutations in \code{permutedTrt}. The (\code{i,j})th element is the importance measure for the \code{r[j]}th variable when the \code{i}th column of each component of \code{permutedTrt} is treated as the original treatment assignment. 

In either case, the \code{"parameters"} attribute will be set, which is a list with the following components: 
\itemize{
\item{\code{kernel}: a string of the name of the kernel, or \code{NULL} when \code{bw=Inf};}
\item{\code{weight.trt}: a numeric vector of treatment weights;}
\item{\code{adjust}: a string of the name of adjustment.}
}

In addition, the \code{grad.kdep} object also has the \code{midp} attribute set to the mid-p-value across the \code{r} dimensions.

For \code{get.p.dd.dw}, the result is an \code{length(r)}-vector of importance measure. 

For \code{p.value.grad.kdep}, the result is a numeric vector of approximate mid-p-values when \code{x} is produced by a call to \code{grad.kdep} with \code{adjust="none"}. Otherwise, a vector of \code{NA_real_}'s. 
}

\note{
When \code{test} is \code{TRUE} and when the number of permutations is large, \code{grad.kdep} will consume a large amount of memory. 

The result from \code{p.value.grad.kdep} is only approximate. There is no guarantee that the approximate mid-p-values lie between 0 and 1. 
}

\references{
Manuscript under preparation 
}
\author{
Long Qu
}

\seealso{\code{\link{mrpp.test}}, \code{\link{permuteTrt}}, \code{\link{bw.smoothp}}, \code{\link{mrppBVS.test}}, \code{\link{mrppBVS}}
}
\examples{

set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
mrpp.obj = mrpp(x, trt=trt, B=5e3L)
(ir1=grad.kdep(mrpp.obj))

ir2=grad.kdep(mrpp.obj, adjust='weighted.mean')
ir3=grad.kdep(mrpp.obj, adjust='scale')
ir4=get.p.dd.dw(x, mrpp.obj$permutedTrt)
ir5=grad.kdep(mrpp.obj, adjust='log scale')
ir6=grad.kdep(mrpp.obj, adjust='logit scale')
ir7=grad.kdep(mrpp.obj, bw=Inf)
cor(cbind(ir1, ir2, ir3, ir4, ir5, ir6, ir7), method="s")

p.value(ir1, type="drop1")
p.value(ir1, type="add1")
p.value(ir1, type="keep1")

p.value(ir5, type="drop1")
p.value(ir5, type="add1")
p.value(ir5, type="keep1")

p.value(ir6, type="drop1")
p.value(ir6, type="add1")
p.value(ir6, type="keep1")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{multivariate}
\keyword{nonparametric}