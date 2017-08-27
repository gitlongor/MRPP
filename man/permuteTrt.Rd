\name{permuteTrt}
\Rdversion{1.1}
\alias{permuteTrt}
\alias{nparts}
\alias{nperms.permutedTrt}
\alias{ntrt.permutedTrt}
\alias{trt.permutedTrt}
\alias{[.permutedTrt}
\title{
Compute permutations of treatment assignments in a one-way design. 
}
\description{
Enumerate or sample permutations of treatment assignments in a one-way design. 

\code{permuteTrt} computes (a sample of) permutations of given a vector of treatment labels, represented as a \code{permutedTrt} object. 

\code{nparts} computes the total number of non-equivalent permutations of treatment labels. 

\code{nperms.permutedTrt} returns the number of permutations generated in a \code{permutedTrt} object. 

\code{ntrt.permutedTrt} returns the sample size in each treatment group of a \code{permutedTrt} object. 

\code{trt.permutedTrt} returns the treatment label vector of a \code{permutedTrt} object. 

\code{[.permutedTrt} returns a subset of the permutations of a \code{permutedTrt} object. 
}
\usage{
permuteTrt(trt, B = 100)
nparts(n)
nperms.permutedTrt(permutedTrt)
ntrt.permutedTrt(permutedTrt)
trt.permutedTrt(permutedTrt)
\method{[}{permuteTrt}(x,i,j)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{trt}{A vector or factor of treatment labels
}
  \item{B}{A positive integer of the number of permutations requested. Currently, a \code{bigz} object is not accepted. 
}
%  \item{idxOnly}{A logical scalar, indicating whether only factoradic indices are computed. This has no effect when \code{B} is not smaller than \code{nparts(table(trt))}. See the Value section for details.}
  
  \item{n}{An integer vector of the sample sizes in each treatment group.}
  
  \item{permutedTrt, x}{A \code{"permutedTrt"} object returned from a \code{permutedTrt} call.}
  
  \item{i}{An integer, character or logical vector of indices of the list to be subset. This is often unspecified. }
  
  \item{j}{An integer vector, specifying the columns of matrices in \code{x} to be extracted. }
  
  %\item{\dots}{Unused.}
}
\details{
The \code{nparts} functions computes the total number of partitions, with sizes being specified by \code{n}. 

Denote the totality by \code{N}, i.e., \code{length(trt)} or \code{sum(n)}. When \code{B} is no smaller than \code{nparts(table(trt))}, all partitions are returned by \code{permuteTrt} and \code{B} is set to \code{nparts(table(trt))}. 
Otherwise, a sample \code{B} random permutations will be selected from a total of \code{nparts(table(trt))} possible non-equivalent permutations without replacement. 
}
\value{
For \code{permuteTrt}, it returns an S3 object of class \code{"permutedTrt"}, which is a named list \code{ANS} of \code{length(unique(trt))} treatment assignment matrices. The names of the list correspond to the levels of \code{trt}. The ordering of the treatment is first determined by the sample size in each treatment (in \emph{decreasing} order) and then by lexical order of the treatment labels ( in \emph{increasing} order). 

%When \code{idxOnly} is \code{FALSE}, 
Each matrix has \code{B} columns, where \code{B} may be smaller than requested, such that \code{ANS[[i]][,b]} contains sorted individual indices that are allocated to treatment \code{i} under the \code{b}th random permutation. Note that, here \code{b} should not be interpreted as the factoradic number. That is, the column locations of the matrices have no special meaning. However, an exception is that the first column (i.e., the first permutation) is always equivalent to the original treatment assignment. In other words, it is guaranteed that \code{setequal( split(seq(length(trt)),trt), lapply(ANS, "[", , 1L))} should always be \code{TRUE}. %In addition, the \code{"idx"} attribute of the result is set to \code{NA_character_}. 

%When \code{idxOnly} is \code{TRUE}, each matrix has only one column, corresponding to the original treatment assignment. The permutations being sampled are represented by the correspoding decimal factoradic numbers (stored as character) returned in the \code{"idx"} attribute. In this case, the length of the \code{"idx"} attribute is \code{B} and the first element is always \code{"0"}, corresponding to no permutations.

\code{nparts} returns a \code{"\link{bigz}"} object which is the total number of non-equivalent treatment assignments. 
\code{nperms.permutedTrt} returns the total number of permutations. 
\code{ntrt.permutedTrt} returns a named vector of sample sizes. 
\code{trt.permutedTrt} returns a vector of treatment labels that is equivalent  to (i.e., might be a permutation of) the one passed to \code{permuteTrt}.

For \code{[.permutedTrt} the result is another \code{permutedTrt} object. If \code{i} is missing, the result only contains a subset of columns of \code{x} specified by \code{j}. If \code{i} is non-missing but \code{j} is missing, the result only contains a subset of treatments, equivalent to subsetting a usual \code{list}. If both \code{i} and \code{j} are specified, the result will subset both the treatments (determined by \code{i}) and the permutations (i.e., columns specified by \code{j}).
}
%\note{
%When \code{B} is less than \code{nparts(table(trt))}, ideally one could randomly sample \code{B} partitions from a total of \code{nparts(table(trt))} without replacement, subject to the first partition being equivalent to the original treatment assignment. For the sake of simplicity, this is \emph{not} implemented here. The random sample of size \code{B} is generated from a total of \code{N!} permutations without replacement. Therefore, it is possible for some of the permutations sampled corresponding to equivalent partitions. However, for the purpose of permutation tests, the type I error rate is still under control. 
%}
\author{
Long Qu
}
\seealso{
\code{\link{combn}}, \code{\link{urand.bigz}}, \code{\link{HSEL.bigz}}, \code{\link{sample1}}, \code{\link{dec2permvec}}
}
\examples{
set.seed(2340)
trt=sample(gl(2,10))[-1] ## two-group comparison with slight imbalance
nparts(table(trt))		## 92378 partitions =  choose(19,10)
pmat=permuteTrt(trt, 1e3L)		## use 1000 random permutations
ntrt.permutedTrt(pmat)  ## 2 treatments
nperms.permutedTrt(pmat) ## 1e3L permutations

pmat[,1L] # extract the 1st permutation of all treatments, i.e., the original treatment assignment
pmat['1', 1:10] # observations in treatment '1' under the first 10 permutations
str(pmat['1']) # all permutations in treatment '1'
stopifnot(identical(pmat, pmat[]))

\dontrun{
pmat=permuteTrt(trt, 1e6L)		## use all partitions, as 1e6L >= 92378 
stopifnot(identical(
  permuteTrt(trt.permutedTrt(pmat), 1e6L), 
  pmat
))
}

}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{utilities}
\keyword{iteration}
\keyword{arith}