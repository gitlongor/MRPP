
.importance <-
function(y, permutedTrt, bw, r=seq_len(NCOL(y)), test=FALSE, distFunc=dist, 
	distExpr.1d = expression(distFunc(y[,index,drop=FALSE])), mrpp.stats=NULL, 
    kernel='triweight', weight.trt="df", measure=c('grad','weighted.mean','scale','p'))
{
	distObj = distFunc(y)
	if(is.null(mrpp.stats)) mrpp.stats=mrpp.test.dist(distObj,permutedTrt=permutedTrt,weight.trt=weight.trt, method='permutation')$all.statistics
	weight.trt = mrpp.weight.trt(weight.trt, trt.permutedTrt(permutedTrt))$weight.trt[names(permutedTrt)]
    B=length(mrpp.stats)
    b=if(isTRUE(test)) 1:B else 1L
    if(!is.matrix(y) && !is.data.frame(y)) y = as.matrix(y)
    ans=matrix(NA_real_, length(b), length(r))
    N=as.integer(nrow(y))
	measure=match.arg(measure)

	pars=list(kernel=kernel, weight.trt=weight.trt, measure=measure)
	if(is.finite(bw)){
		weight = dkernel(kernel)( (mrpp.stats-rep(mrpp.stats[b],each=B) )/ bw)
		dim(weight)=c(B, length(b))
		if(measure=='grad'){
		  weight=weight/bw/B
		}else weight=weight/rep(.colSums(weight, B, length(b)), each=B)
	}
	
	distFunc.1d=function(index)eval(distExpr.1d)
    all.ddelta.dw=sapply(r, distFunc.1d)   
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
	
	if(measure=='scale'){
			w2s=sqrt(.colSums(weight^2, B, length(b)))
			expr = expression(ans[, r.i] <- .colSums(weight * (rep(dz.dw[b], each=B)-dz.dw), B, length(b))/sd(dz.dw)/w2s)
	}else	expr = expression(ans[, r.i] <- .colSums(weight * (rep(dz.dw[b], each=B)-dz.dw), B, length(b)))
	
	for(r.i in seq_along(r)){
		dz.dw=.Call(mrppstats, all.ddelta.dw[,r.i], permutedTrt, weight.trt, PACKAGE='MRPP')
		eval(expr)    ## this lines replace the above 3 lines
	}
    structure(drop(ans), parameters=pars, midp=midp.empirical(mrpp.stats), class='.importance')
}

grad.smoothp <-
function(y, permutedTrt, bw, r=seq_len(NCOL(y)), test=FALSE, 
        distObj=dist(y), mrpp.stats=NULL, 
        kernel='triweight', weight.trt="df", adjust=NULL)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
    ## min.wts=1e-8  ### this was handled by bw.safety()
	if(is.null(mrpp.stats)) mrpp.stats=mrpp.test.dist(distObj,permutedTrt=permutedTrt,weight.trt=weight.trt, method='permutation')$all.statistics
	pval0=midp.empirical(mrpp.stats)
	weight.trt = mrpp.weight.trt(weight.trt, trt.permutedTrt(permutedTrt))$weight.trt[names(permutedTrt)]
    B=length(mrpp.stats)
    b=if(isTRUE(test)) 1:B else 1L
    if(!is.matrix(y) && !is.data.frame(y)) y = as.matrix(y)
    ans=matrix(NA_real_, length(b), length(r))
    N=as.integer(nrow(y))
    #if(missing(cperm.mat)) cperm.mat=apply(permutedTrt,2,function(kk)(1:N)[-kk])

	if(missing(bw)) bw='sym1'
	if(is.character(bw))
		bw=bw.smoothp(y,permutedTrt=permutedTrt,r=r, kernel=kernel, weight.trt=weight.trt, method=bw, verbose=FALSE)

	if(is.null(adjust[1L])) adjust='none'
	if(is.numeric(bw) && is.infinite(bw)) adjust='weighted.mean'
	adjust=match.arg(adjust, c('none','weighted.mean','scale', 'log scale'))

	pars=list(kernel=kernel, weight.trt=weight.trt, adjust=adjust, bw=bw)
	adjust0=adjust; if(adjust=='log scale') adjust='none'
#    weight=matrix(NA_real_, B, length(b))   ## this may require large memory when test=TRUE
#    for(b.i in 1:length(b))
#      weight[,b.i]=pmax(min.wts,dnorm((mrpp.stats[b[b.i]]-mrpp.stats),0,bw))
	weight = dkernel(kernel)( (mrpp.stats-rep(mrpp.stats[b],each=B) )/ bw)
	dim(weight)=c(B, length(b))
	if(adjust=='none'){
		  weight=weight/bw/B
	}else weight=weight/rep(.colSums(weight, B, length(b)), each=B)
	
#    contrast.mat=matrix(0,choose(N,2),N); k=1
#    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
#    #all.ddelta.dw=abs(contrast.mat%*%y)^2/pmax(1e-8,distObj)/2 ## avoiding division by zero
#    all.ddelta.dw=(contrast.mat%*%y)^2/distObj/2 ## when denom is zero, the numerator is also zero. 
#        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
    all.ddelta.dw=apply(y[,r,drop=FALSE],2L,dist)^2/distObj*0.5   ## these 2 lines replace the above 5 lines
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
	
	if(adjust=='scale'){
			w2s=sqrt(.colSums(weight^2, B, length(b)))
			expr = expression(ans[, r.i] <- .colSums(weight * (rep(dz.dw[b], each=B)-dz.dw), B, length(b))/sd(dz.dw)/w2s)
	}else	expr = expression(ans[, r.i] <- .colSums(weight * (rep(dz.dw[b], each=B)-dz.dw), B, length(b)))
	for(r.i in seq_along(r)){
		#dz.dw=.C('mrppstats',all.ddelta.dw[,r.i],permutedTrt,cperm.mat,nrow(permutedTrt),B,N,as.integer(weight.trt),ans=double(B),PACKAGE='MRPP',DUP=FALSE)$ans
		dz.dw=.Call(mrppstats, all.ddelta.dw[,r.i], permutedTrt, weight.trt, PACKAGE='MRPP')
#        for(b.i in 1:length(b)){
#            dd.dw=dz.dw[b[b.i]]-dz.dw
#            ans[r.i, b[b.i]]=sum(weight[,b.i]*dd.dw)/B  #length(b)
#        }
		eval(expr)    ## this lines replace the above 3 lines
	}
	if(adjust0=='log scale') ans=exp(ans/pval0)
    structure(drop(ans), parameters=pars, midp=pval0, class='grad.smoothp')
}

p.value.grad.smoothp = function(x, type=c('keep1','drop1','add1'),...)
{
	adj=attr(x, 'parameters')$adjust
	if(!(adj %in% c('none','log scale'))) return(rep(NA_real_, length(x)))
	type=match.arg(type)
	x0=x; attributes(x0)=NULL
	if(adj=='none'){
		switch(type, 
			drop1 = attr(x, 'midp') - x0, 
			keep1 = attr(x, 'midp') - sum(x0) + x0, 
			add1 = attr(x, 'midp') + x0
		)
	}else if(adj=='log scale'){
		switch(type, 
			drop1 = attr(x, 'midp') / x0, 
			keep1 = attr(x, 'midp') * exp( -sum(log(x0)) + log(x0) ) , 
			add1 = attr(x, 'midp') * x0
		)
	}else stop('"adjust" unsupported')
}

grad.smoothp.bw <-
expression( ## simplified from grad.smoothp; allowing a vector of bw's; only used in bw.smoothp
{
    B=length(mrpp.stats)
    b=1L
	n.bw=length(bw)
    ans=matrix(NA_real_, n.bw, length(r))
    N=as.integer(nrow(y))

	all.ddelta.dw=apply(y[,r,drop=FALSE],2L,dist)^2/distObj*0.5   
	all.ddelta.dw[is.nan(all.ddelta.dw)]=0
	
	bw0=rep(bw, each=B)
	weights=dkernel(kernel)((mrpp.stats-mrpp.stats[b])/bw0)/bw0
	dim(weights) = c(B, n.bw)
	for(r.i in seq_along(r)){
		dz.dw=.Call(mrppstats, all.ddelta.dw[,r.i], permutedTrt, as.numeric(weight.trt), PACKAGE='MRPP')  ## character weight.trt not allowed
		ans[, r.i] = scale/B* .colSums(weights* (dz.dw[b]-dz.dw), B, n.bw)
	}
})


hessian.smoothp <-
function(y, permutedTrt, r=seq_len(NCOL(y)), test=FALSE, 
        distObj=dist(y), 
        mrpp.stats=mrpp.test.dist(distObj,permutedTrt=permutedTrt,weight.trt=weight.trt)$all.statistics,
        kernel='triweight', bw=bw.mse.pdf.asym(mrpp.stats), #cperm.mat, 
        weight.trt="df", scale=1, standardized=FALSE)
## y=N-by-p data matrix; r=dimension index; 
{
    B=length(mrpp.stats)
    b=if(isTRUE(test)) 1:B else 1L
    if(!is.matrix(y) && !is.data.frame(y)) y = as.matrix(y)
    ans=array(NA_real_, dim=c(length(b), length(r), length(r)))
    N=as.integer(NROW(y))
	
	out0 = outer(-mrpp.stats, mrpp.stats[b], '+')
    if(is.finite(bw)) weight = dkernel(kernel)(out0, 0, bw)

    all.ddelta.dw=apply(y[,r,drop=FALSE],2L,dist)^2/distObj*0.5   
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0

	diffs = apply(y[,r, drop=FALSE], 2L, function(z)vech(outer(z,z,'-'))[-(1+cumsum(c(0, safeseq(length(z),2,by=-1))))])
	all.d2delta.dw2=array(NA_real_, dim=c(length(distObj), length(r), length(r)))
	for(r.i in seq_along(r)){
		all.d2delta.dw2[,,r.i]= diffs * diffs[,r.i] / distObj^3 * -0.25
	}
	all.d2delta.dw2[is.nan(all.d2delta.dw2)]=0
	
    if(is.finite(bw)){
		out1s = list()
        for(r.i in seq_along(r)){
            dz.dw=.Call(mrppstats, all.ddelta.dw[,r.i], permutedTrt, as.numeric(weight.trt), PACKAGE='MRPP')
			out1s[[r.i]] = outer(-dz.dw, dz.dw[b], '+')
		}
		for(r.i in seq_along(r)){
			for(s.i in safeseq(r.i, length(r))) {
				dz2=.Call(mrppstats, all.d2delta.dw2[,s.i,r.i], permutedTrt, as.numeric(weight.trt), PACKAGE='MRPP')
				ans[, r.i, s.i] = ans[, s.i, r.i] = 
					scale/B* colSums(weight * 
						outer(-dz2, dz2[b], '+') - 1/bw^2 * out0 * out1s[[r.i]] * out1s[[s.i]]
					)
			}
        }
    }else{  ## infinite bw (the same as finite bw, except no weights(i.e. weight=1)
		.NotYetImplemented()
        if(isTRUE(standardized)){
            for(r.i in seq_along(r)){
                dz.dw=.Call(mrppstats, all.ddelta.dw[,r.i], permutedTrt, as.numeric(weight.trt), PACKAGE='MRPP')
                ans[, r.i] = scale/B* (B*dz.dw[b]-sum(dz.dw))/sd(dz.dw) 
            }
        }else{
            for(r.i in seq_along(r)){
                dz.dw=.Call(mrppstats, all.ddelta.dw[,r.i], permutedTrt, as.numeric(weight.trt), PACKAGE='MRPP')
                ans[, r.i] = scale/B* (B*dz.dw[b]-sum(dz.dw)) 
            }
        }
    }
    drop(ans)
}

if(FALSE){
    weight=matrix(NA_real_, B, length(b))   ## this may require large memory when test=TRUE
    for(b.i in 1:length(b))
      weight[,b.i]=pmax(min.wts,dnorm((mrpp.stats[b[b.i]]-mrpp.stats),0,bw))
#    weight = dnorm(outer(mrpp.stats, mrpp.stats[b], '-'), 0, bw)    ## this replaces the above 3 lines


    contrast.mat=matrix(0,choose(N,2),N); k=1
    for(i in 1:(N-1))for(j in (i+1):N){contrast.mat[k,i]=1;contrast.mat[k,j]=-1;k=k+1}
    #all.ddelta.dw=abs(contrast.mat%*%y)^2/pmax(1e-8,distObj)/2 ## avoiding division by zero
    all.ddelta.dw=(contrast.mat%*%y)^2/distObj/2 ## when denom is zero, the numerator is also zero. 
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
#    all.ddelta.dw=apply(y,2L,dist)^2/distObj*0.5   ## these 2 lines replace the above 5 lines
#        all.ddelta.dw[is.nan(all.ddelta.dw)]=0


        for(b.i in 1:length(b)){
            dd.dw=dz.dw[b[b.i]]-dz.dw
            ans[r.i, b[b.i]]=sum(weight[,b.i]*dd.dw)/B  #length(b)
        }
#        ans[, r.i] = colMeans(weight * outer(-dz.dw, dz.dw[b], '+'))    ## this lines replace the above 3 lines
}