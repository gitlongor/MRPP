\name{MRPP importance}
\Rdversion{1.1}
\alias{grad.smoothp}
\alias{hessian.smoothp}
\alias{get.p.dd.dw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation-based variable importance measure
}
\description{
\code{grad.smoothp} computes variable importance measure of partial derivative of continuous approximation to MRPP p-values with respect to hypothetical weights. 
\code{get.p.dd.dw} computes variable importance measure by using partial derivative of distance function with respect to hypothetical weights as the new distance function in MRPP
}
\usage{
grad.smoothp(y, permutedTrt, bw, r=seq_len(NCOL(y)), test=FALSE, 
        distObj=dist(y), mrpp.stats=NULL, 
        kernel='biweight', weight.trt="df", adjust='none')
get.p.dd.dw(y, permutedTrt, r = seq_len(ncol(y)), weight.trt="df", eps=1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ The response data matrix with each row being a data point and each column being a variable. 
}
  \item{permutedTrt}{A list matricies that list all random treatment assignments to be used in the MRPP test. 
}
  \item{r}{A positive integer vector, specifying variables for which variable importance measure is requested. 
}
  \item{test}{\code{TRUE} or \code{FALSE}. If \code{FALSE}, the first permutation in \code{permutedTrt} will be used as the original treatment assignment. If \code{TRUE}, each permutation in \code{permutedTrt} will be treated as if it is the original treatment assignment, useful for doing permutation tests using the computed variable importance measure as the test statistic.
}
  \item{distObj}{The numeric vector computed from \code{dist(y)}, i.e., the lower triangle of the Euclidean distance matrix. 
}
  \item{mrpp.stats}{A vector of all observed MRPP test statistics. If this is \code{NULL}, then 
  
 \code{mrpp.test.dist(distObj,permutedTrt=permutedTrt,weight.trt=weight.trt)$all.statistics}

  will be called to obtain it.  
}
  \item{kernel}{A character scalar among \code{c("gaussian", "biweight", 'triweight', 'tricube',"logistic")}
}
  \item{bw}{A numeric scalar giving the bandwidth used for kernel density estimation of \code{mrpp.stats}. This might be positive infinity. 
}
%  \item{cperm.mat}{}
  \item{weight.trt}{As in \code{\link{mrpp.test.dist}}.
}
  \item{adjust}{A character scalar being one of \code{c('none','weighted.mean','scale')}. See details. 
}
%  \item{min.wts}{}
  \item{eps}{A small positive tolerance, the same as described in \code{\link{mrpp.test}}. }
}

\details{
\code{grad.smoothp} first uses kernel density estimator to approximate the distribution of \code{mrpp.stats} and obtain a continuous approximate p-value for the MRPP test. Then it computes the partial derivative of the approximate p-value with respect to the hypothetical weight of each variable, evaluated at weight 1. 
That is, if the partial derivative is negatively large, then the corresponding variable will be more important, because increasing its weight will quickly decreasing the approximate MRPP p-value. Vice versa. 

When \code{adjust='none'}, the partial derivatives will be returned. When \code{adjust='weighted.mean'}, the partial derivatives will be treated as a weighted average, i.e., will be divided by total weights. When \code{adjust='scale'}, they are further divided standard deviations. 

\code{get.p.dd.dw} first computes the partial derivative of the Euclidian distance measure with respect to the hypothetical weights for each variable. This is then used as the distance function for the variable to be used in MRPP. The resulting MRPP p-value is returned as the variable importance. The smaller, the more important.
}
\value{
For \code{grad.smoothp}, when \code{test} is \code{FALSE}, the result is an \code{length(r)}-vector of importance measure. 
When \code{test} is \code{TRUE}, the result is an \code{B}-by-\code{length(r)} matrix of importance measure, where \code{B} is the number of permutations in \code{permutedTrt}. The (\code{i,j})th element is the importance measure for the \code{r[j]}th variable when the \code{i}th column of each component of \code{permutedTrt} is treated as the original treatment assignment. 

In either case, the \code{'call'} attribute is set to the call of \code{grad.smoothp}. 

For \code{get.p.dd.dw}, the result is an \code{length(r)}-vector of importance measure. 
}
\note{
When \code{test} is \code{TRUE} and when the number of permutations is large, \code{grad.smoothp} will consume a large amount of memory. 
}

\references{
Long Qu, Dan Nettleton, Jack C. M. Dekkers. Relative Variable Importance and Backward Variable Selection for the Multiresponse Permutation Procedure, with Applications to High Dimensional Genomic Data. (Manuscript under review). 
}
\author{
Long Qu
}

\seealso{\code{\link{mrpp.test}}, \code{\link{permuteTrt}}, \code{\link{bw.grad.smoothp}}, \code{\link{mrppBVS.test}}, \code{\link{mrppBVS}}
}
\examples{

set.seed(2340)
x=matrix(rnorm(20*5),20)
trt=gl(2,10)
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 5e3L)		## use 5000 random permutations
(ir1=grad.smoothp(x, pmat))
grad.smoothp(x, pmat, adjust='weighted.mean')
grad.smoothp(x, pmat, adjust='scale')
ir2=get.p.dd.dw(x, pmat)
cor(ir1, ir2, method="s")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{multivariate}
\keyword{nonparametric}