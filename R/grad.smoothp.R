grad.kdep.mrpp <-
function(y,  bw, kernel='triweight', adjust=NULL, mrpp.stats=NULL, r=seq_len(y$R), test=FALSE)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
	on.exit({ # recover original data.env
			if(!is.null(y$data.env$y.bakGradKdeP)){
				y$data.env$y=y$data.env$y.bakGradKdeP
				rm(y.bakGradKdeP, envir=y$data.env)
			}
	})
	if(missing(bw)) bw='sym1'
	if(length(r)!=y$R){
		y[['data.env']]$y.bakGradKdeP=y$y
		y[['data.env']]$y=y[['data.env']]$y[,r,drop=FALSE]
		y$R=length(r)
		y$distObj=y$distFunc(y$y)
		r=seq_len(y$R)

	}	
	if(is.numeric(bw)&&is.infinite(bw)&&bw>0){
		this.call=as.list(match.call())
		if('bw'%in%names(this.call)) this.call[['bw']]=NULL
		if('adjust'%in%names(this.call)) this.call[['adjust']]=NULL
		if('kernel'%in%names(this.call)) this.call[['kernel']]=NULL
		if('r'%in%names(this.call)) this.call[['r']]=NULL
		this.call[[1L]]=NULL
		return(do.call(grad.kdep.Inf.mrpp, this.call, envir=parent.frame()))
	}
	if(!isTRUE(test)){
		this.call=as.list(match.call())
		if('test'%in%names(this.call))this.call[['test']]=NULL
		if('r'%in%names(this.call))this.call[['r']]=NULL
		this.call[[1L]]=NULL
		return(do.call(grad.kdep.notest.mrpp, this.call,envir=parent.frame()))
	}

	if(is.null(mrpp.stats)) {
		mrppt=mrpp.test(y, test.method='permutation'); 
		mrpp.stats=mrppt$all.statistics
	}
	
	weight.trt = y$weight.trt
    B=length(mrpp.stats)
    b=if(isTRUE(test)) 1:B else 1L
    ans=matrix(NA_real_, length(b), y$R)
    N=y$nobs

	if(is.null(adjust[1L])) adjust='none'
	if(is.numeric(bw) && is.infinite(bw)) adjust='weighted.mean'
	adjust=match.arg(adjust, c('none','weighted.mean','scale', 'log scale', 'logit scale'))

	if(is.character(bw))
		bw=bw.kdep(y,kernel=kernel, bw.method=bw, verbose=FALSE, 
			scale=switch(adjust, none=,weighted.mean=,scale='raw',`log scale`='log',`logit scale`='logit'),
			mrpp.stats=mrpp.stats
		)


	pval0=p.empirical(mrpp.stats,midp=FALSE)-as.numeric(0.5/y$nparts)

	
	pars=list(kernel=kernel, weight.trt=weight.trt, adjust=adjust, bw=bw)
	adjust0=adjust; if(adjust=='log scale' || adjust=='logit scale') adjust='none'
#    weight=matrix(NA_real_, B, length(b))   ## this may require large memory when test=TRUE
#    for(b.i in 1:length(b))
#      weight[,b.i]=pmax(min.wts,dnorm((mrpp.stats[b[b.i]]-mrpp.stats),0,bw))
	weight = dkernel(kernel)( (mrpp.stats-rep(mrpp.stats[b],each=B) )/ bw)
	#dim(weight)=c(B, length(b))
	w.idx=which(weight>0)
	if(adjust=='none'){
		  weight=weight/bw/B
	}else weight=weight/rep(.colSums(weight, B, length(b)), each=B)
	
#    contrast.mat=matrix(0,choose(N,2),N); k=1
#    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
#    #all.ddelta.dw=abs(contrast.mat%*%y)^2/pmax(1e-8,distObj)/2 ## avoiding division by zero
#    all.ddelta.dw=(contrast.mat%*%y)^2/distObj/2 ## when denom is zero, the numerator is also zero. 
#        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
    all.ddelta.dw=apply(y$y,2L,y$distFunc)^2/y$distObj*0.5   ## these 2 lines replace the above 5 lines
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
	
	if(adjust=='scale'){
			w2s=sqrt(.colSums(weight^2, B, length(b)))
			expr = expression(ans[, r.i] <- .colSums(weight * (rep(dz.dw[b], each=B)-dz.dw), B, length(b))/sd(dz.dw)/w2s)
	}else	expr = expression(ans[, r.i] <- .colSums(weight * (rep(dz.dw[b], each=B)-dz.dw), B, length(b)))
	for(r.i in r){ # r=seq_len(R)
		#dz.dw=.C('mrppstats',all.ddelta.dw[,r.i],permutedTrt,cperm.mat,nrow(permutedTrt),B,N,as.integer(weight.trt),ans=double(B),PACKAGE='MRPP',DUP=FALSE)$ans
		dz.dw=.Call(mrppstats, all.ddelta.dw[,r.i], y$permutedTrt, weight.trt, PACKAGE='MRPP')
#        for(b.i in 1:length(b)){
#            dd.dw=dz.dw[b[b.i]]-dz.dw
#            ans[r.i, b[b.i]]=sum(weight[,b.i]*dd.dw)/B  #length(b)
#        }
		eval(expr)    ## this lines replace the above 3 lines
	}
	if(adjust0=='log scale') {
		ans=exp(ans/pval0) 
	}else if(adjust0=='logit scale') ans=exp(ans/pval0/(1-pval0))
	#eval(recover.data.env)
    structure(drop(ans), parameters=pars, midp=pval0, class=c('grad.kdep','grad.mrpp'))
}


## Liuhua's very fast version with bw=Inf, because no permutations are needed: 
grad.kdep.Inf.mrpp <-
function(y, test=FALSE, mrpp.stats=NULL)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{	bw=Inf; adjust='weighted.mean'; kernel=NULL

	if(is.null(mrpp.stats)) {
		mrppt=mrpp.test(y, test.method='permutation'); 
		mrpp.stats=mrppt$all.statistics
		pval0=p.value(mrppt, 'midp')
	}else pval0=p.empirical(mrpp.stats,midp=FALSE)-as.numeric(0.5/y$nparts)
	weight.trt = y$weight.trt
    B=length(mrpp.stats)
    b=if(isTRUE(test)) 1:B else 1L
    ans=matrix(NA_real_, length(b), y$R)
    N=y$nobs

	pars=list(kernel=kernel, weight.trt=weight.trt, adjust=adjust, bw=bw)
	
	all.ddelta.dw=apply(y$y,2L,y$distFunc)^2/y$distObj*0.5
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
	
	for(r.i in seq_len(y$R)){
		ans[, r.i]=.Call(mrppstats_subset, all.ddelta.dw[,r.i], y$permutedTrt, weight.trt, b, PACKAGE='MRPP') - mean(all.ddelta.dw[,r.i])
	}
	
    structure(drop(ans), parameters=pars, midp=pval0, class=c('grad.kdep','grad.mrpp'))
}



grad.kdep.notest.mrpp <-
function(y, bw, mrpp.stats=NULL, kernel='triweight', adjust=NULL)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{	if(missing(bw)) bw='sym1'
	if(is.null(mrpp.stats)) {
		mrppt=mrpp.test(y, test.method='permutation'); 
		mrpp.stats=mrppt$all.statistics
		pval0=p.value(mrppt, 'midp')
	}else pval0=p.empirical(mrpp.stats,midp=FALSE)-as.numeric(0.5/y$nparts)
	weight.trt = y$weight.trt
    B=length(mrpp.stats)
    ans=matrix(NA_real_, 1L, y$R)
    N=y$nobs

	if(is.null(adjust[1L])) adjust='none'
	if(is.numeric(bw) && is.infinite(bw)) adjust='weighted.mean'
	adjust=match.arg(adjust, c('none','weighted.mean','scale', 'log scale', 'logit scale'))

	if(is.character(bw)) # this slow when B is large
		bw=bw.kdep(y, kernel=kernel,  bw.method=bw, verbose=FALSE, 
			scale=switch(adjust, none=,weighted.mean=,scale='raw',`log scale`='log',`logit scale`='logit'),
			mrpp.stats=mrpp.stats
		)


	pars=list(kernel=kernel, weight.trt=weight.trt, adjust=adjust, bw=bw)
	adjust0=adjust; if(adjust=='log scale' || adjust=='logit scale') adjust='none'
	weight = dkernel(kernel)( (mrpp.stats-mrpp.stats[1L] )/ bw)
	w.idx=which(weight>0)
	w.pos=weight[w.idx]
	if(adjust=='none'){
		  w.pos=w.pos/bw/B
	}else w.pos=w.pos/sum(w.pos)
	
    all.ddelta.dw=apply(y$y,2L,y$distFunc)^2/y$distObj*0.5   
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
	
	if(adjust=='scale'){
			w2s=sqrt(sum(w.pos^2))
			expr = expression(ans[, r.i] <- sum(w.pos * (dz.dw[1L]-dz.dw) )/sd(dz.dw)/w2s)
	}else	expr = expression(ans[, r.i] <- sum(w.pos * (dz.dw[1L]-dz.dw) ) )

	#the loop below is sufficiently fast, when compact kernel is used (i.e., w.idx is short)
	for(r.i in seq_len(y$R)){
		dz.dw=.Call(mrppstats_subset, all.ddelta.dw[,r.i], y$permutedTrt, weight.trt, w.idx, PACKAGE='MRPP')
		eval(expr)    
	}

	if(adjust0=='log scale') {
		ans=exp(ans/pval0)
	}else if(adjust0=='logit scale') ans=exp(ans/pval0/(1-pval0))
    structure(drop(ans), parameters=pars, midp=pval0, class=c('grad.kdep','grad.mrpp'))
}

grad.kdep = grad.kdep.mrpp
