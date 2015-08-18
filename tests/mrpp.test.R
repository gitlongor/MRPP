(R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != '')

library(combinat)
mrpp.stats.R=function(y, trt, wt='df', round=8)
{
	trt=as.factor(trt)
	tabtrt=table(trt)
	ntrt=length(tabtrt)
	N=sum(tabtrt)
	levs=levels(trt)
	
	wt=switch(wt, 
		df = (tabtrt-1), 
		n =  (tabtrt-0), 
		equal= rep(1,ntrt), 
		pairs=choose(tabtrt, 2)
	)
	wt = wt/sum(wt)
	distMat=as.matrix(dist(y))
	
	sum1trt=function(idx)sum(distMat[idx,idx])
	stat=function(trt){
		idxes = outer(seq(ntrt), trt, '==')
		sums = apply(idxes, 1L, sum1trt)
		sum(sums * wt / tabtrt / pmax(1,tabtrt - 1) )
	}
	all.stats = unlist(permn(as.integer(trt), stat))
	sort(unique(round(all.stats, 8L)))
}

library(MRPP)



do.test=function(n=c(3,5), wts=c('df','n','equal','pairs'), R=5)
{
	trt=ordered(sample(rep(seq_along(n), n)), levels=as.character(seq_along(n)))
	y=matrix(rnorm(length(trt)*R), nc=R)
	pmat=permuteTrt(trt, Inf)		## use 1000 random permutations
	
	for(wt in wts){
		mrpp.rslt=mrpp.test(y, trt, permutedTrt=pmat, weight.trt=wt, method='permutation')
		mrpp.stat.rslt = mrpp.stats.R(y, trt, round=8L, wt=wt)
		stopifnot(identical(
			sort(unique(round(mrpp.rslt$all.statistics, 8L))), 
			mrpp.stat.rslt
		))
		cat('weight = ', wt, ':\tOK\n')
	}
	invisible(TRUE)
}

set.seed(2340)
if(!R_CHECK_TIMINGS_){
	for(i in seq(1e1L)){
		ntrt = sample(2:4, 1)
		n=sample(rpois(ntrt, 2)+1)
		R = sample(rpois(1, 10)+1)
		if(max(nparts(n), factorial(sum(n)))<1e7L)	{
			cat("sample sizes = ", n, "\n")
			do.test(n)
		}
	}
}
