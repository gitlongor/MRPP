## this did not take account of the dependence of 
## Pearson3 distribution on the hypothetical weights 
## in distance measure through the permutational 
## moments. 
## Temporarily give up this idea...

grad.pear3p.mrpp <-
function(y, adjust=NULL, r=seq_len(y$R), test=FALSE)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
	on.exit({ # recover original data.env
			if(!is.null(y$data.env$y.bakGradPear3P)){
				y$data.env$y=y$data.env$y.bakGradPear3P
				rm(y.bakGradPear3P, envir=y$data.env)
			}
	})
	if(length(r)!=y$R){
		y[['data.env']]$y.bakGradPear3P=y$y
		y[['data.env']]$y=y[['data.env']]$y[,r,drop=FALSE]
		y$R=length(r)
		y$distObj=y$distFunc(y$y)
		r=seq_len(y$R)

	}	
	if(!isTRUE(test)){
		this.call=as.list(match.call())
		if('test'%in%names(this.call))this.call[['test']]=NULL
		if('r'%in%names(this.call))this.call[['r']]=NULL
		this.call[[1L]]=NULL
		return(do.call(grad.kdep.notest.mrpp, this.call,envir=parent.frame()))
	}

	mrppt=mrpp.test(y, test.method='pearson3'); 
	mrpp.z=mrppt$statistic
	pval0=p.value(mrppt)
	
	weight.trt = y$weight.trt
    B=length(mrpp.stats)
    b=if(isTRUE(test)) 1:B else 1L
    ans=matrix(NA_real_, length(b), y$R)
    N=y$nobs

	if(is.null(adjust[1L])) adjust='none'
	adjust=match.arg(adjust, c('none', 'log scale', 'logit scale'))


	pars=list(weight.trt=weight.trt, adjust=adjust)
	
	
    all.ddelta.dw=apply(y$y,2L,y$distFunc)^2/y$distObj*0.5
        all.ddelta.dw[is.nan(all.ddelta.dw)]=0
	dz.dw=numeric(y$R)
	permtrt1=permutedTrt1(y$trt)
	for(r.i in r){
		dz.dw[r.i]=.Call(mrppstats, all.ddelta.dw[,r.i], permtrt1, weight.trt, PACKAGE='MRPP')
	}
	
	cums=cumulant.mrpp(y,order=1:3)
	cums[2L]=sqrt(cums[2L]); cums[3L]=cums[3L]/cums[2L]^3
	shape = 4/cums[3]^2; scale = abs(cums[3L]*.5)
	mrpp.gam=sign(cums[3L])*((mrpp.z-cums[1L])/cums[2L]+2/cums[3L])
	# the derivative below did not take account of the dependence of pearson3 distribution on the introduced weights in the distance measure
	ans = dgamma(mrpp.gam, shape=shape, scale=scale)/cums[2L]*dz.dw
	
	if(adjust=='log scale') {
		ans=exp(ans/pval0) 
	}else if(adjust0=='logit scale') ans=exp(ans/pval0/(1-pval0))

    structure(drop(ans), parameters=pars, midp=pval0, class=c('grad.kdep','grad.mrpp'))
}
