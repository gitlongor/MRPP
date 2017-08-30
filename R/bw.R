# starting from Scott's/SJ's rule as prelim bw
# then using kernels with this bw to estimate f(x[x]) and f''(x[1]) 
bw.mse.pdf.asym.SJ=function(x,iter.max=1L,eps=1e-6,start.bw, kernel='triweight', verbose=FALSE)
{
	if(missing(start.bw)) {
		start.bw=bw.SJ(x)*sqrt(mkernel(kernel,2L,R=FALSE)) 
		start.bw=function(x)bw.nrd(x)*sqrt(mkernel(kernel,order=2L,R=FALSE))
	}
    if(is.function(start.bw)) bw0=start.bw(x)
    if(is.numeric(start.bw)) bw0=start.bw
    if(is.character(start.bw)) bw0=call(start.bw, x)

	kr.r.R=krRkernel(kernel);
		kr=kr.r.R['kr'];  Rk=kr.r.R['R']
		stopifnot(kr.r.R['r']==2)
	dkern=dkernel(kernel)
	d2dkern=d2dkernel(kernel)
    
    n.iter=1L
    xdiff=x[1L]-x
    repeat{
            tmp=xdiff/bw0
            f.1=mean(dkern(tmp))/bw0
            ddf.1=mean(d2dkern(tmp))/bw0^3
        bw=(f.1*Rk/(length(x)*4*(kr*ddf.1)^2))^.2
        if(abs(bw-bw0)<eps || n.iter>=iter.max) return(bw[[1L]])
        if(verbose && isTRUE(n.iter%%verbose==0))cat(bw,fill=TRUE)
        bw0=bw
        n.iter=n.iter+1
    }
}

# starting from a bw that matches Pearson Type III dist'n
# then using kernels with this bw to estimate f(x[x]) and f''(x[1]) 
bw.mse.pdf.asym.pearson3gca=function(x, mrpp,iter.max=1L,eps=1e-6,kernel='triweight', verbose=FALSE)
{
	cums=cumulant(mrpp, order=1:4)
	cums[2L]=sqrt(cums[2L])
	cums[3L]=cums[3L]/cums[2L]^3
	cums[4L]=cums[4L]/cums[2L]^4

	pdf0=	dpearson3gca(x, cums[1L], cums[2L], cums[3L], cums[4L])
	bw0=bw.matchpdf(x, kernel=kernel, pdf=pdf0)
	
	kr.r.R=krRkernel(kernel);
		kr=kr.r.R['kr'];  Rk=kr.r.R['R']
		stopifnot(kr.r.R['r']==2)
	dkern=dkernel(kernel)
	d2dkern=d2dkernel(kernel)
    
    n.iter=1L
    xdiff=x[1L]-x
	#bw0=Inf 
	#f.1=dpearson3(x[1], cums[1L], cums[2L], cums[3L])
	#ddf.1=d2dpearson3(x[1], cums[1L], cums[2L], cums[3L])
    repeat{
		tmp=xdiff/bw0
		f.1=mean(dkern(tmp))/bw0
		ddf.1=mean(d2dkern(tmp))/bw0^3

		bw=(f.1*Rk/(length(x)*4*(kr*ddf.1)^2))^.2
        if(abs(bw-bw0)<eps || n.iter>=iter.max) return(bw[[1L]])
        if(verbose && isTRUE(n.iter%%verbose==0))cat(bw,fill=TRUE)
        bw0=bw
        n.iter=n.iter+1
    }
}


# Using Wand & Jone's plug-in method to obtain bw for f() and for f''() separately
# Then using kernels with such estimated bw's to obtain f(x[1]) and f''(x[1]) 
# Plug-in final results to Parzen's formula
# This is often similar to the unbiased CV below, but much faster.
bw.mse.pdf.asym.WJpi=function(x, kernel='triweight', verbose=FALSE)
{
	h0=ks::hpi(x, deriv.order=0) #gaussian kernel bw
	h2=ks::hpi(x, deriv.order=2) #gaussian kernel bw
	
	kr.r.R=krRkernel(kernel);
		kr=kr.r.R['kr'];  Rk=kr.r.R['R']
		stopifnot(kr.r.R['r']==2)
	dkern.gau=dkernel('gaussian')
	d2dkern.gau=d2dkernel('gaussian')
    
    xdiff=x[1L]-x
            tmp=xdiff/h0
            f.1=mean(dkern.gau(tmp))/h0
			tmp=xdiff/h2
            ddf.1=mean(d2dkern.gau(tmp))/h2^3
			
			bw=(f.1*Rk/(length(x)*4*(kr*ddf.1)^2))^.2
	if(isTRUE(verbose)){
		attr(bw, 'bw.f(x[1])')=h0
		attr(bw, "bw.f''(x[1])")=h2
	}
	bw
}



# Using unbiased CV to obtain bw for f() and for f''() separately, with a max of 500 data values
# Then using kernels with such estimated bw's to obtain f(x[x]) and f''(x[1]) 
bw.mse.pdf.asym.ucv=function(x, kernel='triweight', verbose=FALSE)
{
	nmax=512L;	nx=length(x); iter.max=1L; eps=1e-6
	x.idx=if(nx<nmax) seq_len(nmax) else unique(as.integer(round(seq.int(from=1L, to=nx, length.out=nmax))))
	bw00=kedd::h.ucv(x[x.idx], deriv.order=0L, kernel=kernel)$h*(nmax/nx)^.2
	bw02=kedd::h.ucv(x[x.idx], deriv.order=2L, kernel=kernel)$h*(nmax/nx)^(1/9)
	
	kr.r.R=krRkernel(kernel);
		kr=kr.r.R['kr'];  Rk=kr.r.R['R']
		stopifnot(kr.r.R['r']==2)
	dkern=dkernel(kernel)
	d2dkern=d2dkernel(kernel)
    
    xdiff=x[1L]-x
	f.1=mean(dkern(xdiff/bw00))/bw00
	ddf.1=mean(d2dkern(xdiff/bw02))/bw02^3

	bw=(f.1*Rk/(length(x)*4*(kr*ddf.1)^2))^.2
	
	if(isTRUE(verbose)){
		attr(bw, 'bw.f(x[1])')=bw00
		attr(bw, "bw.f''(x[1])")=bw02
	}
	bw
}

bw.mse.pdf.asym =bw.mse.pdf.asym.WJpi

bw.range=function(x, length=50L, lower=.05, upper=.95, safety=2)
{
	uniqstat=sort(unique(x))
	lo = quantile(diff(uniqstat), lower)
	hi = diff(quantile(uniqstat, prob=c(lower, upper)))
	lf = log10(safety)
	10^seq.int(from=log10(lo)-lf, to=log10(hi)+lf, length.out=length)
}

bw.safety=function(x, kernel, nNonzero=3L, pdf.cut=1e-3)
{# avoiding bw that is too small to be useful
	supp=skernel(kernel, pdf.cut)[2L]
	
	diffs=x[1L]-x
	u.a.d=unique(abs(diffs))
	ord=order(u.a.d)
	n.u.a.d=length(u.a.d)
	ans0 = if(nNonzero<n.u.a.d) u.a.d[ord[nNonzero+1L]] else tail(u.a.d,1L)
	ans0 / supp
}

bw.smoothp.grid <-
function(y, permutedTrt, r=seq_len(NCOL(y)), bw = NULL, 
	distFunc=dist,  kernel='triweight', 
	bw.method='sym1', verbose=TRUE, subset, adjust='none',...)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
	bw.method=match.arg(bw.method, choices=c('sym1','drop1','add1','keep1','mse(z[1])','pearson3','pearson3gca')) # 'dropadd1','dropaddsym1', 'ss.gradp'
	adjust=match.arg(adjust, choices=c('none','log scale'))
    if(!is.matrix(y) && !is.data.frame(y)) y = as.matrix(y)
	R=NCOL(y)
	if(R==1L ) bw.method='ss.gradp'

	lst=list(...)
	lst$y=y[,r,drop=FALSE]
	lst$permutedTrt=permutedTrt
	lst$distFunc=distFunc 
	nms = setdiff(names(formals(mrpp)), '...')
	idx=names(lst)%in%nms
	mrpp.obj = do.call('mrpp', lst[idx])
	mrpp=mrpp.test(mrpp.obj, method='permutation')
	distObj=mrpp.obj$distObj
	
#	distObj = distFunc(y[,r,drop=FALSE])
#	lst=list(...)
#	lst$y = distObj
#	lst$permutedTrt=permutedTrt
#	nms = setdiff(c(names(formals(mrpp.test.dist)), names(formals(mrpp))), '...')
#	idx=names(lst)%in%nms
#	mrpp = do.call('mrpp.test.dist', c(lst[idx],list(method='permutation')))
	
	lower.bound = bw.safety(mrpp$all.statistics, kernel)
	
	plot.expr=expression({
		hist(mrpp$all.statistics, breaks='FD', freq=FALSE, main='Density of MRPP statistics', ylab='Density', xlab='MRPP z-statistics')
		pdf.final=dkde(mrpp$all.statistics, ans, kernel)
		marg = diff(range(mrpp$all.statistics))/5
		curve(pdf.final(x), min(mrpp$all.statistics)-marg, max(mrpp$all.statistics)+marg, add=TRUE)
		rug(mrpp$all.statistics)
		abline(v=mrpp$all.statistics[1L], col='blue', lty=3L)
		axis(1L, at=mrpp$all.statistics[1L], labels=expression(z[1L]), col='blue', line=1, col.ticks='blue', col.axis='blue')
		box()
		title(sub=sprintf('bandwidth = %.3g', ans))
	})
	if(bw.method=='mse(z[1])'){
		lst=list(...)
		lst$x=mrpp$all.statistics
		#lst$mrpp=mrpp.obj
		lst$kernel=kernel
		#lst$verbose=verbose
		nms=names(lst)
		idx= nms=='' | nms%in%names(formals(bw.mse.pdf.asym))
		ans = do.call('bw.mse.pdf.asym', lst[idx])
		if(ans<lower.bound){
			ans=lower.bound
			if(verbose)warning('bandwidth selected is too small; replaced by a safer lower bound.')
		}
		if(verbose) {
			dev.new(noRStudioGD=TRUE)
			eval(plot.expr)
			title(main='\nMethod: KDE MSE formula')
		}
		return(ans)
	}

	weight.trt = attr(mrpp$parameter, 'weight.trt')
	
	if(missing(bw) || is.null(bw)) bw = bw.range(mrpp$all.statistics)	
	bw=bw[bw>=lower.bound]
	n.bw=length(bw)
	if(n.bw==0L) stop('bandwidth range is not suitable')
	
	if(bw.method=='pearson3' || bw.method=='pearson3gca') {
		cums=cumulant(mrpp(y[,r,drop=FALSE], permutedTrt=permutedTrt, weight.trt=weight.trt, distFunc=distFunc),order=1:4); 
		cums[2L]=sqrt(cums[2L]); cums[3L]=cums[3L]/cums[2L]^3; cums[4L]=cums[4L]/cums[2L]^4
		pdf0x = if(bw.method=='pearson3') dpearson3(mrpp$all.statistics, cums[1L], cums[2L], cums[3L]) else 
		dpearson3gca(mrpp$all.statistics, cums[1L], cums[2L], cums[3L], cums[4L]) 
		ans=bw.matchpdf(mrpp$all.statistics, kernel=kernel, pdf=pdf0x, bw=start.bw, verbose=verbose, title=switch(bw.method, pearson3="Match Pearson III Dist'n", pearson3gca="Match GCA Adjust. of Pearson III Dist'n"))
		if(ans<lower.bound){
			ans=lower.bound
			if(verbose) warning('bandwidth selected is too small; replaced by a safer lower bound.')
		}
		return(ans)

	}
	
	## pre-computing all.ddelta.dw in grad.smoothp
	r.bak=r
	if(bw.method%in%c('sym1','drop1','add1','keep1','dropadd1','dropaddsym1','ss.gradp')) {
		if(missing(subset)){
			n.subset=round(100+.2*(length(r)-100))
			if(length(r)<=100L) {
				subset=seq_along(r)
			}else{
				this.call=as.list(match.call())
				this.call[[1L]]=
				this.call$bw.method=
				this.call$verbose=
				this.call$subset=
				this.call$kernel=
				this.call$adjust=NULL
				this.call$bw=Inf
				this.call$permutedTrt=lapply(permutedTrt, '[',,1L,drop=FALSE)
				imp.idx=names(this.call)%in%names(formals(grad.smoothp))
				imp0=do.call('grad.smoothp.Inf', this.call[imp.idx], envir=parent.frame())
				ord0=order(imp0)
				subset=ord0[round(seq.int(from=1L, to=length(r), length.out=n.subset))]
			}
		}else if(is.null(subset)){
			subset=seq_along(r)
		}
	    r=r[subset]
	}
	
	gradEnv=new.env()
	gradEnv$mrpp.stats=mrpp$all.statistics
	gradEnv$scale=1
	eval(grad.smoothp.bw, envir=gradEnv)
	
	if(bw.method=='ss.gradp'){
		ss = rowSums(gradEnv$ans^2)
		idx = which.max(ss)
		ans=bw[idx]
		if(verbose){
			dev.new(width=10, height=5, noRStudioGD=TRUE); par(mfrow=1:2)
			plot(log10(bw), (ss), xlab='bandwidth', ylab='SS of p-value gradients', type='o', axes=FALSE, main='\nMethod: Max SS of Derivative of Density' )
			axis(2)
			axis(1, at = log10(ans), labels=sprintf('%.3g',ans), col='blue',col.ticks='blue', col.axis='blue', line=1)
			axis(2, at = ss[idx], labels=sprintf('%.3g',ss[idx]), col='red',col.ticks='red', col.axis='red', line=1)
			ats=axTicks(1)
			axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
			abline(v=log10(ans), col='blue', lty=3)
			abline(h=(ss[idx]), col='red', lty=3)
			box()
			
			eval(plot.expr)
		}
		if(idx==1L || idx==length(bw))warning('optimal bandwidth occurs at the boundary')
		ans =  exp(mean(log(ans)))  
		if(ans<lower.bound){
			ans=lower.bound
			if(verbose)warning('bandwidth selected is too small; replaced by a safer lower bound.')
		}
		return(ans)
		
	}


	mrpp.obj1=mrpp.obj
	if(bw.method%in%c('drop1','sym1','dropadd1', 'dropaddsym1')){
		drop1p = function(r.i){
			#lst$y=distFunc(y[,r[-r.i],drop=FALSE])
			#p.value(do.call('mrpp.test.dist', lst),type="midp")
			mrpp.obj1$distObj=distFunc(y[,r[-r.i],drop=FALSE])
			mrpp.obj1$R=mrpp.obj$R-1L
			p.value(mrpp.test(mrpp.obj1), type='midp')
		}
		drop1pval= sapply(seq_along(r), drop1p)
	}
	if(bw.method%in%c('add1','sym1','dropadd1','dropaddsym1')){
		add1p = function(r.i){
			#lst$y=distFunc(y[,c(r,r[r.i]),drop=FALSE])
			#p.value(do.call('mrpp.test.dist', lst),type="midp")
			mrpp.obj1$distObj=distFunc(y[,c(r,r[r.i]),drop=FALSE])
			mrpp.obj1$R=mrpp.obj$R+1L
			p.value(mrpp.test(mrpp.obj1), type='midp')
			
		}
		add1pval=sapply(seq_along(r), add1p)
	}
	if(bw.method%in%c('keep1')){
		keep1 = function(r.i){
			#lst$y=distFunc(y[,r[r.i],drop=FALSE])
			#p.value(do.call('mrpp.test.dist', lst),type="midp")
			mrpp.obj1$distObj=distFunc(y[,r[r.i],drop=FALSE])
			mrpp.obj1$R=1L
			p.value(mrpp.test(mrpp.obj1), type='midp')
		}
		keep1pval=sapply(seq_along(r), keep1)
		#t(smps - gradEnv$ans%*%(1-diag(1, length(r), length(r))))
	}
	pval0=p.value(mrpp, type='midp')
	if(adjust=='none'){
		sses = switch(bw.method, 
			drop1= .rowSums((gradEnv$ans - (pval0 - drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) ,
			add1 = .rowSums((gradEnv$ans - (add1pval - pval0)[col(gradEnv$ans)])^2, n.bw, length(r)) , 
			sym1 = .rowSums((gradEnv$ans - .5*(add1pval- drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) , 
			keep1= .rowSums((gradEnv$ans%*%(1-diag(1, length(r), length(r)))-(pval0-keep1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) , 
			dropadd1 = .rowSums((gradEnv$ans - (pval0 - drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) 
				+ .rowSums((gradEnv$ans - (add1pval - pval0)[col(gradEnv$ans)])^2, n.bw, length(r)) ,
			dropaddsym1 = .rowSums((gradEnv$ans - (pval0 - drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) 
				+ .rowSums((gradEnv$ans - (add1pval - pval0)[col(gradEnv$ans)])^2, n.bw, length(r))
				+ .rowSums((gradEnv$ans - .5*(add1pval- drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) 
		)
	}else if(adjust=='log scale'){
		log.p0=log(pval0)
		gradEnv$ans = gradEnv$ans/pval0
		sses = switch(bw.method, 
			drop1= .rowSums((gradEnv$ans - (log.p0 - log(drop1pval))[col(gradEnv$ans)])^2, n.bw, length(r)) ,
			add1 = .rowSums((gradEnv$ans - (log(add1pval) - log.p0)[col(gradEnv$ans)])^2, n.bw, length(r)) , 
			sym1 = .rowSums((gradEnv$ans - .5*(log(add1pval)- log(drop1pval))[col(gradEnv$ans)])^2, n.bw, length(r)) , 
			keep1= .rowSums((gradEnv$ans%*%(1-diag(1, length(r), length(r)))-(log.p0-log(keep1pval))[col(gradEnv$ans)])^2, n.bw, length(r)) , 
			dropadd1 = .rowSums((gradEnv$ans - (log.p0 - log(drop1pval))[col(gradEnv$ans)])^2, n.bw, length(r)) 
				+ .rowSums((gradEnv$ans - (log(add1pval) - log.p0)[col(gradEnv$ans)])^2, n.bw, length(r)) ,
			dropaddsym1 = .rowSums((gradEnv$ans - (log.p0 - log(drop1pval))[col(gradEnv$ans)])^2, n.bw, length(r)) 
				+ .rowSums((gradEnv$ans - (log(add1pval) - log.p0)[col(gradEnv$ans)])^2, n.bw, length(r))
				+ .rowSums((gradEnv$ans - .5*(log(add1pval)- log(drop1pval))[col(gradEnv$ans)])^2, n.bw, length(r)) 
		)
	}
	idx=which.min(sses)
	ans=bw[idx]
	min.sse=sses[idx]

	if(ans<lower.bound){
		ans=lower.bound
		if(verbose)warning('bandwidth selected is too small; replaced by a safer lower bound.')
	}
	
	if(verbose){
		dev.new(width=10, height=5, noRStudioGD=TRUE)
		par(mfrow=1:2)

		plot(log10(bw), sses, xlab='bandwidth', ylab='SS of approx errors', type='o', main='', axes=FALSE)
		title(main=paste0('\nMethod: ',switch(bw.method, drop1='Backward Difference', add1='Forward Difference', keep1='Keep 1 Variable', sym1='Central Difference', dropadd1='Backward + Forward Difference', dropaddsym1='Backward + Forward + Central Difference')))
		#axis(3, at = log10(ans), labels=sprintf('%.3g',ans), col='red',col.ticks='red',col.axis='red')
		axis(2)
		ats=axTicks(1)
		axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
		abline(h=min.sse, col='red', lty=3); box()
		
		axis(1, at = log10(ans), labels=sprintf('%.3g',ans), col='blue',col.ticks='blue', col.axis='blue', line=1)
		axis(2, at = min.sse, labels=sprintf('%.3g',min.sse), col='red',col.ticks='red', col.axis='red', line=1)
		abline(v=log10(ans), col='blue', lty=3)

		eval(plot.expr)
	}
	
	ans
	
}

bw.smoothp.optim.mrpp <-
function(y, start.bw = NULL, kernel='triweight', 
	bw.method='mse(z[1])', scale='raw', verbose=TRUE, n.subset, mrpp.stats=NULL, r=seq_len(y$R))
## y mrpp object; r=dimension index; 
{
	on.exit({ # recover original data.env
			if(!is.null(y$data.env$y.bak)){
				y$data.env$y=y$data.env$y.bak
				rm(y.bak, envir=y$data.env)
			}
	})
	if(length(r)!=y$R){
		y$data.env$y.bak = y$y
		y[['data.env']]$y=y[['data.env']]$y[,r,drop=FALSE]
		y$R=length(r)
		y$distObj=y$distFunc(y$y)
		r=seq_len(y$R)
	}

	if(is.null(mrpp.stats)) {
		mrppt=mrpp.test(y, method='permutation'); 
		mrpp.stats=mrppt$all.statistics
	}	
	bw.method=match.arg(bw.method, choices=c('sym1','drop1','add1','keep1','mse(z[1])','pearson3','pearson3gca')) # 'dropadd1','dropaddsym1','ss.gradp',
	scale=match.arg(scale, choices=c('raw','log'))
    
	R=y$R
	if(R==1L ) bw.method='ss.gradp'

#	lst=list(...)
#	lst$y=y[,r,drop=FALSE]
#	lst$permutedTrt=permutedTrt
#	lst$distFunc=distFunc 
#	nms = setdiff(names(formals(mrpp)), '...')
#	idx=names(lst)%in%nms
#	mrpp.obj = do.call('mrpp', lst[idx])
#	mrpp=mrpp.test(mrpp.obj, method='permutation')
#	distObj=mrpp.obj$distObj
	
	
	lower.bound = bw.safety(mrpp.stats, kernel)
	
	plot.expr=expression({
		hist(mrpp.stats, breaks='FD', freq=FALSE, main='Density of MRPP statistics', ylab='Density', xlab='MRPP z-statistics')
		pdf.final=dkde(mrpp.stats, bw=ans, kernel=kernel)
		marg = diff(range(mrpp.stats))/5
		curve(pdf.final(x), min(mrpp.stats)-marg, max(mrpp.stats)+marg, add=TRUE)
		rug(mrpp.stats)
		abline(v=mrpp.stats[1L], col='blue', lty=3L)
		axis(1L, at=mrpp.stats[1L], labels=expression(z[1L]), col='blue', line=1, col.ticks='blue', col.axis='blue')
		box()
		title(sub=sprintf('bandwidth = %.3g', ans))
	})
	if(bw.method=='mse(z[1])'){
		ans=bw.mse.pdf.asym(mrpp.stats, kernel=kernel,verbose=verbose)
		if(ans<lower.bound){
			ans=lower.bound
			if(verbose)warning('bandwidth selected is too small; replaced by a safer lower bound.')
		}
		if(verbose) {
			dev.new(noRStudioGD=TRUE)
			eval(plot.expr)
			tmp=density(mrpp.stats, bw=attr(ans,'bw.f(x[1])'))
			lines(tmp$x, tmp$y, col=2)
			tmp=density(mrpp.stats, bw=attr(ans,"bw.f''(x[1])"))
			lines(tmp$x, tmp$y, col=5)
			legend('topleft', lty=1L, col=c(1:2,5), text.col=c(1:2,5), legend=c('final density est.', 'initial density', 'initial 2nd deriv. density'))
			title(main='\n\nMethod: KDE MSE formula',cex.main=.8)
		}
#		eval(recover.data.env)
		return(ans)
	}

	weight.trt = y$weight.trt
	
	if(bw.method=='pearson3' || bw.method=='pearson3gca') {
		cums=cumulant(y,order=1:4); 
		cums[2L]=sqrt(cums[2L]); cums[3L]=cums[3L]/cums[2L]^3; cums[4L]=cums[4L]/cums[2L]^4
		pdf0x = if(bw.method=='pearson3') dpearson3(mrpp.stats, cums[1L], cums[2L], cums[3L]) else 
		dpearson3gca(mrpp.stats, cums[1L], cums[2L], cums[3L], cums[4L]) 
		ans=bw.matchpdf(mrpp.stats, kernel=kernel, pdf=pdf0x, bw=NULL, verbose=verbose, title=switch(bw.method, pearson3="Match Pearson III Dist'n", pearson3gca="Match GCA Adjust. of Pearson III Dist'n"))
		if(ans<lower.bound){
			ans=lower.bound
			if(verbose) warning('bandwidth selected is too small; replaced by a safer lower bound.')
		}
#		eval(recover.data.env)
		return(ans)
	}
	
	r.bak=r
	if(bw.method%in%c('sym1','drop1','add1','keep1','dropadd1','dropaddsym1','ss.gradp')) {
		if(missing(n.subset))
			n.subset= if(R>100L) round(100+.2*(R-100)) else R
		if(n.subset>R) n.subset=R
		if(n.subset<R) {
			imp0=grad.smoothp.Inf.mrpp(y=y,test=FALSE,mrpp.stats=mrpp.stats)
			ord0=order(imp0)
			r=unique(sort(ord0[round(seq.int(from=1L, to=R, length.out=n.subset))]))

			if(is.null(y$data.env$y.bak)) y$data.env$y.bak = y$y
			y[['data.env']]$y=y[['data.env']]$y[,r,drop=FALSE]
			y$R=length(r)
			y$distObj=y$distFunc(y$y)
			r=seq_len(y$R)			
		}
	}
	
	gradEnv=new.env()
	gradEnv$mrpp.stats=mrpp.stats
	evalq({
		B=length(mrpp.stats)
		N=as.integer(y$nobs)
	
		all.uni.dist2=apply(y$y,2L,y$distFunc)^2
		all.ddelta.dw=all.uni.dist2/y$distObj*0.5   
		all.ddelta.dw[is.nan(all.ddelta.dw)]=0
		
		dz.dw=matrix(NA_real_, B, y$R)
		worked.idx=rep(FALSE, B)
	}, envir=gradEnv)
	
	thisEnv=environment()
	locs=numeric(0L)
	niter=0L
	mrpp.stats.diff1=mrpp.stats-mrpp.stats[1L]
	
	
	# starting values
	if(is.null(start.bw)){
		start.bw = bw.smoothp.optim.mrpp(y, kernel=kernel, bw.method='mse(z[1])', mrpp.stats=mrpp.stats, verbose=FALSE)
	}else start.bw=start.bw[1L]

	dkern=dkernel(kernel)
	ddkern=ddkernel(kernel)
	obj.common.expr=quote({
		niter<<-niter+1L
		thisEnv$locs[[niter]]=this.bw=exp(logbw)
		weights=dkern(mrpp.stats.diff1/this.bw)
		w.tf=weights>0; w.pos=weights[w.tf]
		perm.idx=which(!gradEnv$worked.idx & w.tf)
		if(length(perm.idx)>0L){
			gradEnv$worked.idx[perm.idx]=TRUE
			for(r.i in seq_along(r))
				gradEnv$dz.dw[perm.idx,r.i]=.Call(mrppstats_subset, gradEnv$all.ddelta.dw[,r.i], y$permutedTrt, weight.trt, perm.idx, PACKAGE='MRPP')
		}
		
		ans=1/this.bw/gradEnv$B*(sum(w.pos)*gradEnv$dz.dw[1,]-colSums(w.pos*gradEnv$dz.dw[w.tf,,drop=FALSE]))
	})

	if(bw.method=='ss.gradp'){ 
# 		# this heuristic does not seem to work at all
#		objfunc=eval(bquote(function(logbw)
#		{
#			niter<<-niter+1L
#			thisEnv$locs[[niter]]=this.bw=exp(logbw)
#			ddkde.fit=ddkde(mrpp.stats, bw=this.bw, kernel=kernel)
#			ans=sum(ddkde.fit(mrpp.stats)^2)
#			thisEnv$ss[[niter]]=ans
#			-ans
#		}))
		objfunc=eval(bquote(function(logbw){
			.(obj.common.expr)
			ans = sum(ans^2) # sum(abs(ans)^4)
			thisEnv$ss[[niter]]=ans
			-ans
		}))
		attr(objfunc, 'srcref')=NULL
		opt.rslt=suppressWarnings(optim(log(start.bw), objfunc, method='Nelder-Mead'))
		ord=order(locs)
		bw=locs[ord]
		ss=ss[ord]
		
		idx = which.max(ss)
		ans=bw[idx]
		if(verbose){
			dev.new(width=10, height=5, noRStudioGD=TRUE); par(mfrow=1:2)
			plot(log10(bw), (ss), xlab='bandwidth', ylab='SS of p-value gradients', type='o', axes=FALSE, main='\nMethod: Max SS of Derivative of Density' )
			axis(2)
			axis(1, at = log10(ans), labels=sprintf('%.3g',ans), col='blue',col.ticks='blue', col.axis='blue', line=1)
			axis(2, at = ss[idx], labels=sprintf('%.3g',ss[idx]), col='red',col.ticks='red', col.axis='red', line=1)
			ats=axTicks(1)
			axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
			abline(v=log10(ans), col='blue', lty=3)
			abline(h=(ss[idx]), col='red', lty=3)
			box()
			
			eval(plot.expr)
		}
		if(idx==1L || idx==length(bw))warning('optimal bandwidth occurs at the boundary')
		ans =  exp(mean(log(ans)))  
		if(ans<lower.bound){
			ans=lower.bound
			if(verbose)warning('bandwidth selected is too small; replaced by a safer lower bound.')
		}
#		eval(recover.data.env)
		return(ans)
		
	}


	mrpp.obj1=y
	distObj2=y$distObj^2
	if(any(bw.method==c('drop1','sym1'))){
		mrpp.obj1$R=y$R-1L
		drop1p = function(r.i){
			#mrpp.obj1$distObj=y$distFunc(y$data.env$y[,-r.i,drop=FALSE])
			#mrpp.obj1$distObj=sqrt(distObj2-y$distFunc(y$data.env$y[,r.i,drop=FALSE])^2)
			mrpp.obj1$distObj=sqrt(distObj2-gradEnv$all.uni.dist2[,r.i])
			mrpp.test.mrpp(mrpp.obj1, method='pearson3')$p.value
		} # >97% of time is spent in mrpp.test.mrpp
		drop1pval= sapply(r, drop1p)
	}
	if(any(bw.method==c('add1','sym1'))){
		mrpp.obj1$R=y$R+1L
		add1p = function(r.i){
			#mrpp.obj1$distObj=y$distFunc(y$data.env$y[,c(r,r.i),drop=FALSE])
			#mrpp.obj1$distObj=sqrt(distObj2+y$distFunc(y$data.env$y[,r.i,drop=FALSE])^2)
			mrpp.obj1$distObj=sqrt(distObj2+gradEnv$all.uni.dist2[,r.i])
			mrpp.test.mrpp(mrpp.obj1,method='pearson3')$p.value
		} 
		add1pval=sapply(r, add1p)
	}
	if(bw.method=='keep1'){
		mrpp.obj1$R=1L
		keep1 = function(r.i){
			mrpp.obj1$distObj[]=sqrt(gradEnv$all.uni.dist2[,r.i])
			mrpp.test.mrpp(mrpp.obj1,method='pearson3')$p.value
		}
		keep1pval=sapply(r, keep1)
		#t(smps - gradEnv$ans%*%(1-diag(1, length(r), length(r))))
	}
	
	pval0=p.empirical(mrpp.stats, midp=FALSE) - as.numeric(0.5/y$nparts)
	if(scale=='raw'){
		sses.expr=switch(bw.method, 
			drop1 = quote(sum((ans - (pval0 - drop1pval))^2)), 
			add1  = quote(sum((ans - (add1pval - pval0))^2)), 
			sym1  = quote(sum((ans - .5*(add1pval - drop1pval))^2)), 
			keep1 = quote(sum((ans%*%(1-diag(1, length(r), length(r)))-(pval0-keep1pval))^2)) #, 
#			dropadd1 = quote(sum((ans - (pval0 - drop1pval))^2)+
#							 sum((ans - (add1pval - pval0))^2)), 
#			dropaddsym1 = quote(sum((ans - (pval0 - drop1pval))^2)+
#							    sum((ans - (add1pval - pval0))^2)+
#								sum((ans - .5*(add1pval - drop1pval))^2))
		)
		objfunc=eval(bquote(function(logbw){
			.(obj.common.expr)
			thisEnv$sses[[niter]]=.(sses.expr)
		}))
	}else if(scale=='log'){
		log.p0=log(pval0)
		sses.expr=switch(bw.method, 
			drop1 = quote(sum((ans - (log.p0 - log(drop1pval)))^2)), 
			add1  = quote(sum((ans - (log(add1pval) - log.p0))^2)), 
			sym1  = quote(sum((ans - .5*(log(add1pval)- log(drop1pval)))^2)), 
			keep1 = quote(sum((ans%*%(1-diag(1, length(r), length(r)))-(log.p0-log(keep1pval)))^2)) #, 
#			dropadd1 = quote(sum((ans - (log.p0 - log(drop1pval)))^2)+
#							 sum((ans - (log(add1pval) - log.p0))^2)), 
#			dropaddsym1 = quote(sum((ans - (log.p0 - log(drop1pval)))^2)+
#							    sum((ans - (log(add1pval) - log.p0))^2)+
#								sum((ans - .5*(log(add1pval)- log(drop1pval)))^2))
		)
		objfunc=eval(bquote(function(logbw){
			.(obj.common.expr)
			ans=ans/pval0
			thisEnv$sses[[niter]]=.(sses.expr)
		}))
	}else stop('unknown "scale"')
		attr(objfunc, 'srcref')=NULL
		if(FALSE){
			opt.rslt=suppressWarnings(optim(log(start.bw), objfunc, method='Nelder-Mead'))
		}else {
			opt.rg.const=c(.5*(1-sqrt(5)),1)*log(5)
			opt.range=log(start.bw)+opt.rg.const
			repeat{
				opt.rslt=optimize(objfunc, opt.range)
				if(any(abs(opt.rslt$minimum-opt.range)< 1e-3)){
					opt.range=opt.rslt$minimum+opt.rg.const
					next
				}
				break
			}
		}
		ord=order(locs)
		bw=locs[ord]
		sses=sses[ord]
	idx=which.min(sses)
	ans=bw[idx]
	min.sse=sses[idx]

	if(ans<lower.bound){
		ans=lower.bound
		if(verbose)warning('bandwidth selected is too small; replaced by a safer lower bound.')
	}
	
	if(verbose){
		dev.new(width=10, height=5, noRStudioGD=TRUE)
		par(mfrow=1:2)

		plot(log10(bw), sses, xlab='bandwidth', ylab='SS of approx errors', type='o', main='', axes=FALSE)
		title(main=paste0('\nMethod: ',switch(bw.method, drop1='Backward Difference', add1='Forward Difference', keep1='Keep 1 Variable', sym1='Central Difference', dropadd1='Backward + Forward Difference', dropaddsym1='Backward + Forward + Central Difference')))
		#axis(3, at = log10(ans), labels=sprintf('%.3g',ans), col='red',col.ticks='red',col.axis='red')
		axis(2)
		ats=axTicks(1)
		axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
		abline(h=min.sse, col='red', lty=3); box()
		title(sub=paste('#iterations:', niter))
		
		axis(1, at = log10(ans), labels=sprintf('%.3g',ans), col='blue',col.ticks='blue', col.axis='blue', line=1)
		axis(2, at = min.sse, labels=sprintf('%.3g',min.sse), col='red',col.ticks='red', col.axis='red', line=1)
		abline(v=log10(ans), col='blue', lty=3)

		eval(plot.expr)
	}
	
#	eval(recover.data.env)
	ans
	
}

bw.smoothp=bw.smoothp.optim.mrpp

bw.matchpdf=function(x, kernel=.kernels, pdf, bw = NULL, verbose=FALSE, title=NULL)
{## TODO: Minimize Hellinger Distance
	if(missing(bw) || is.null(bw)) bw = bw.range(x)
	if(is.function(pdf)) pdf=pdf(x)
	
	#sse=function(bw) sum((dkde(x, bw, kernel)(x)-pdf)^2)
	#ss = sapply(bw, sse)
	if(FALSE){ ## grid search
		pdfs=dkde(x, bw, kernel)
		dens=sapply(pdfs, do.call, list(v=x))
		ss=.colSums((dens-pdf)^2, length(x), length(bw))
		idx=which.min(ss)
		ans=bw[idx]
		min.ss=ss[idx]
	}else{# golden sectioning
		pdf0=dkde(x, median(bw), kernel)
		env=attr(pdf0, 'environment')
		this.env=environment()
		niter=0L
		locs=numeric(0L)
		ss=numeric(0L)
		func=function(logbw){
			niter <<- niter+1L
			this.env$locs[[niter]] = thisbw = exp(logbw)
			zetalstar = env$ftkern(thisbw*env$sl) * env$Yl
			zetak=fft(zetalstar)*c(1,-1)
			kde.est=approx(env$tk0, pmax.int(0,Re(zetak)), x, ties='ordered')$y
			this.env$ss[[niter]] = sum((kde.est-pdf)^2)
		}
		opt.rg.const=c(.5*(1-sqrt(5)),1)*log(10)
		if(length(bw)>1L){
			log.range=log(range(bw))
		}else{
			log.range=log(bw)+opt.rg.const
		}
		repeat{
			opt.rslt=optimize(func, log.range)
			if(any(abs(opt.rslt$minimum-log.range)<1e-3)) {
				log.range=opt.rslt$minimum+opt.rg.const
				next
			}
			ans=exp(opt.rslt$minimum)
			min.ss=opt.rslt$objective
			ord=order(locs)
			bw=locs[ord]
			ss=ss[ord]
			break
		}
	}
	if(verbose){
		dev.new(width=10, height=5, noRStudioGD=TRUE)
		par(mfrow=1:2)
		plot(log10(bw), log10(ss), xlab='bandwidth', ylab='log10(SSE of p-value density)', type='o', axes=FALSE )
		if(!is.null(title))graphics::title(main=paste0('\nMethod: ', title)) else graphics::title(main='Matching Initial PDF')
		axis(2L)
		axis(1L, at = log10(ans), labels=sprintf('%.3g',ans), col='blue',col.ticks='blue', col.axis='blue', line=1)
		ats=axTicks(1)
		axis(1L, at=ats, labels=parse(text=paste0('10^',ats))  )
		abline(v=log10(ans), h=log10(min.ss), col='blue', lty=3L)
		title(sub=paste("#iterations:", niter))
		box()
		
		hist(x, breaks='FD', freq=FALSE, main='Density of MRPP statistics', ylab='Density', xlab='MRPP z-statistics')
		pdf.final=dkde(x, ans, kernel)
		marg = diff(range(x))/5
		curve(pdf.final(x), min(x)-marg, max(x)+marg, add=TRUE)
		rug(x)
		abline(v=x[1L], col='blue', lty=3L)
		axis(1L, at=x[1L], labels=expression(z[1L]), col='blue', line=1, col.ticks='blue', col.axis='blue')
		ord.x=order(x)
		lines(x[ord.x], pdf[ord.x], col=2L)
		legend('topleft', lty=1L, col=1:2, text.col=1:2, legend=c('kernel density', 'initial density'))
		box()
		title(sub=sprintf('bandwidth = %.3g', ans))
	}
	if(any( ans==range(bw) ) )warning('optimal bandwidth occurs at the boundary')
	ans
}


if(FALSE){

smoothp=function(x, bw=bw.mse.pdf.asym, kernel = c("gaussian", "biweight", 'triweight', 'tricube',"logistic"), eps=1e-8)
{
    if(is.character(bw)) bw=get(bw, mode='function')
    if(is.function(bw))  bw=bw(x)
    stopifnot(is.numeric(bw))
	if(is.character(kernel)) {
		Kernel = pkernel(kernel)
	}else if(is.function(kernel)) Kernel = kernel
	stopifnot(Kernel(-Inf)==0 && Kernel(Inf)==1)

	midp = p.empirical(x, eps)
	knots=x; n=length(knots); rm(x, eps)
	f0=function(x)
	{	nx=length(x)
		xdiff=x-rep(knots, each=nx)
		dim(xdiff) = c(nx, n)
		.rowMeans(Kernel(xdiff/bw), nx, n)
	}
	structure(list(p.value = f0(knots[1L]), statistic=knots[1L], cdf=f0, midp=midp, kernel=kernel), class = 'smoothp')
}



bw.constr.midp=function(x, eps = 1e-8, verbose=FALSE)
{## this won't work when x[1]==min(x)
	midp = mean(x<=x[1L]-eps) + .5* mean(x>x[1L]-eps & x<x[1L]+eps)
	B=length(x)
	
	xdiff=x[1]-x
	
	eqn =function(logbw){
		bw=exp(logbw)
		mean(pnorm(xdiff/bw)) - midp
	}
	
	lans = nlsolve(log(bw.nrd(x)), eqn, 
		control=list(verbose = FALSE, maxit=25L, tol=sqrt(.Machine$double.eps) )
	)
	exp(lans)
}

    f.1=kdde(x, bw0, deriv.order=0, eval.points=x[1])$estimate
    ddf.1=kdde(x, bw0, deriv.order=2, eval.points=x[1])$estimate

    locpoly(x,  drv = 0L, degree=0, kernel = "normal", gridsize=2L, range.x=0:1+x[1], bandwidth=bw0)$y

    f.1=mean(dnorm((x[1]-x)/bw0))/bw0


        f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
        ddf=drvkde(x,2,bw0,se=FALSE)
        ddf.1=approx(ddf$x.grid[[1]],ddf$est,xout=x[1])$y
        f.1;ddf.1
		

	set.seed(23490)
	x=rnorm(1000)
	bw.mse.pdf.asym(x)
	bw.constr.midp(x)
	if(interactive()) {
	plot(density(x, bw.mse.pdf.asym(x)), ylim=c(0,dnorm(0)))
	lines(density(x, bw.constr.midp(x)), ylim=c(0,dnorm(0)))
	abline(v=x[1L])
	curve(dnorm, lty=3, col="grey", add=TRUE)
	}
	
	xx=seq.int(from=min(x), to=max(x), length.out=1e5)
	midpKde=kde(x, ,bw.constr.midp(x), eval.points=xx)
	
}
