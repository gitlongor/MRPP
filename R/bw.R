if(FALSE){
# starting from Scott's/SJ's rule as prelim bw
# then using kernels with this bw to estimate f(x[x]) and f''(x[1]) 
bw.amse.pdf.SJ=function(x,iter.max=1L,eps=1e-6,start.bw, kernel='triweight', verbose=FALSE)
{
	if(missing(start.bw)) {
		start.bw=bw.SJ(x)*sqrt(mkernel(kernel,2L))
		start.bw=function(x)bw.nrd(x)*sqrt(mkernel(kernel,order=2L))
	}
    if(is.function(start.bw)) bw0=start.bw(x)
    if(is.numeric(start.bw)) bw0=start.bw
    if(is.character(start.bw)) bw0=call(start.bw, x)

	kr=mkernel(kernel,2L)/2; Rk=Rkernel(kernel,0L)
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
bw.amse.pdf.pearson3gca=function(x, mrpp,iter.max=1L,eps=1e-6,kernel='triweight', verbose=FALSE)
{
	cums=cumulant(mrpp, order=1:4)
	cums[2L]=sqrt(cums[2L])
	cums[3L]=cums[3L]/cums[2L]^3
	cums[4L]=cums[4L]/cums[2L]^4

	pdf0=	dpearson3gca(x, cums[1L], cums[2L], cums[3L], cums[4L])
	bw0=bw.matchpdf(x, kernel=kernel, pdf=pdf0)
	
	kr=mkernel(kernel,2L)/2; Rk=Rkernel(kernel,0L)
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
}

# Using Wand & Jone's plug-in method to obtain bw for f() and for f''() separately
# Then using kernels with such estimated bw's to obtain f(x[1]) and f''(x[1]) 
# Plug-in final results to Parzen's formula
# This is often similar to the unbiased CV below, but much faster.
bw.amse.pdf.WJpi=function(x, kernel='triweight', verbose=FALSE, method='WJ',sd,skew)
{
	method0=method=match.arg(method, c('WJ','pearson3'))
	if(method0=='pearson3'){
		if(missing(sd)) sd=stats::sd(x)
		if(missing(skew)) skew=moments::skewness(x)
	}
	
  ## f''(x[1]) esitmation bandwidth
	mu2k.3w=mkernel('triweight',2L) # optimal for 2nd pdf derivative estimation
	Rk.3w=Rkernel('triweight',deriv.order=2L)
	d2dkern=d2dkernel('triweight')
	if(method0=='pearson3'){
		if(4/skew^2-4.5 <=0 ) method='WJ'
	}
	if(method=='WJ'){
		mu2k.gau=mkernel('gaussian',2L)
		Rk.gau=Rkernel('gaussian',deriv.order=c(0L,2L))
		h2=hpi(x, deriv.order=2) #gaussian kernel bw
			h2=h2*(  (Rk.3w/mu2k.3w^2)/(Rk.gau['2']/mu2k.gau^2) )^(1/9)
	}else if(method=='pearson3'){
		h2=bw.amise.pdf(x,kernel='triweight',deriv.order=2L, verbose=FALSE, method='pearson3', sd=sd, skew=skew)
	}else .NotYetImplemented()
	
  ## f(x[1]) esitmation bandwidth
	mu2k.ep=mkernel('epanechnikov',2L) # optimal for density estimation
	Rk.ep=Rkernel('epanechnikov',deriv.order=0L)
	dkern=dkernel('epanechnikov')
	method=method0
	if(method0=='pearson3'){
		if(4/skew^2-2.5 <=0 ) method='WJ'
	}
	if(method=='WJ'){
		h0=hpi(x, deriv.order=0) #gaussian kernel bw
			h0=h0*(  (Rk.ep/mu2k.ep^2)/(Rk.gau['0']/mu2k.gau^2) )^.2
	}else if(method=='pearson3'){
		h0=bw.amise.pdf(x,kernel='epanechnikov',deriv.order=0L, verbose=FALSE, method='pearson3', sd=sd, skew=skew)
	}else .NotYetImplemented()
	

  ## final bandwidth
    Rk.kern=Rkernel(kernel, deriv.order=0L)
	mu2k.kern=mkernel(kernel, 2L)
    xdiff=x[1L]-x
	f.1=mean(dkern(xdiff/h0))/h0
	ddf.1=mean(d2dkern(xdiff/h2))/h2^3
	
	bw=(f.1*Rk.kern/(length(x)*(mu2k.kern*ddf.1)^2))^.2
  ## output
	if(isTRUE(verbose)){
		attr(bw, 'bw.f(x[1])')=h0
		attr(bw, 'kernel.f(x[1])')='epanechnikov'
		attr(bw, "bw.f''(x[1])")=h2
		attr(bw, "kernel.f''(x[1])")='triweight'
	}
	bw
}



# Using unbiased CV to obtain bw for f() and for f''() separately, with a max of 500 data values
# Then using kernels with such estimated bw's to obtain f(x[x]) and f''(x[1]) 
bw.amse.pdf.ucv=function(x, kernel='triweight', verbose=FALSE)
{
	nmax=512L;	nx=length(x)
	x.idx=if(nx<nmax) seq_len(nmax) else unique(as.integer(round(seq.int(from=1L, to=nx, length.out=nmax))))
	bw00=kedd::h.ucv(x[x.idx], deriv.order=0L, kernel='epanechnikov')$h*(nmax/nx)^.2
	bw02=kedd::h.ucv(x[x.idx], deriv.order=2L, kernel='triweight')$h*(nmax/nx)^(1/9)
	
	kr=mkernel(kernel,2L)/2; Rk=Rkernel(kernel,0L)
	dkern=dkernel('epanechnikov')
	d2dkern=d2dkernel('triweight')
    
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

bw.amse.pdf =bw.amse.pdf.WJpi
# Using Wand & Jone's plug-in method to obtain bw for optimal asym. mean integrated squared error assuming gaussian kernel
# Then transform the result to other kernels
bw.amise.pdf=function(x, kernel='triweight', deriv.order=0L, verbose=FALSE, method='WJ', sd,skew)
{
	method=match.arg(method, c('WJ','pearson3'))
	r=deriv.order
	if (method=='pearson3'){
		if(missing(sd)) sd=stats::sd(x)
		if(missing(skew)) skew=moment::skewness(x)
		beta.arg1 = (2*r+5)/2
		if((beta.arg2 = 4/skew^2-beta.arg1 ) <=0 ) method='WJ'
	}
	if(method=='WJ'){
		mu2k=mkernel(kernel,2L)
		mu2k.gau=mkernel('gaussian',2L)
		Rk=Rkernel(kernel,deriv.order=r)
		Rk.gau=Rkernel('gaussian',deriv.order=r)
		h0=hpi(x, deriv.order=r) #gaussian kernel bw
		ans=h0*(  (Rk/mu2k^2)/(Rk.gau/mu2k.gau^2) )^(1/(2*r+5))
		if(isTRUE(verbose)) attr(ans, 'bw.gaussian')=h0
		ans
	}else if (method=='pearson3'){
		.5*sd*abs(skew)*exp((1/(2*r+5))*(
			log(2*base::pi*(1+2*r))+log(Rkernel(kernel,r))-
			log(length(x))-lbeta(beta.arg1 ,beta.arg2)-
			2*log(mkernel(kernel, 2L))
		))
	}else .NotYetImplemented()
}

bw.range=function(x, length=50L, lower=.05, upper=.95, safety=2)
{
	uniqstat=sort(unique(x))
	lo = quantile(diff(uniqstat), lower)
	hi = diff(quantile(uniqstat, prob=c(lower, upper)))
	lf = log10(safety)
	10^seq.int(from=log10(lo)-lf, to=log10(hi)+lf, length.out=length)
}

# upper limit is an upper bound of optimal amise bandwidth
bw.range=function(x, kernel, length=15L, safety=1.2)
{
	lo=bw.safety(x,kernel)
	mu2k=mkernel(kernel,2L)
	Rk=Rkernel(kernel,0L)
	N=length(x)
	#hi=(243/35*Rk/mu2k^2/N)^.2*sd(x) # original formula
	sdx.ucl=sqrt(var(x)*(N-1)/qchisq( .005, N-1)) # 1% upp.conf.lim.
	hi=(243/35*Rk/mu2k^2/N)^.2* sdx.ucl
	10^seq.int(from=log10(lo), to=log10(hi*safety), length.out=length)
}

# upper limit is determined by percent variance inflation
bw.range=function(x, kernel, length=15L, var.infl=0.1)
{
	lo=bw.safety(x,kernel)
	mu2k=mkernel(kernel,2L)
	N=length(x)
	var.n=var(x)*(N-1)/N
	hi=sqrt(var.infl*var.n/mu2k)
	10^seq.int(from=log10(lo), to=log10(hi), length.out=length)
}

bw.safety=function(x, kernel, nNonzero=3L, pdf.cut=1e-3)
{# avoiding bw that is too small to be useful
	supp=skernel(kernel, pdf.cut)[2L]
	
	diffs=x[1L]-x
	u.a.d=unique(abs(diffs))
	ord=order(u.a.d)
	n.u.a.d=length(u.a.d)
	if(nNonzero<n.u.a.d) {
		ans0 = u.a.d[ord[nNonzero+1L]] 
		strict = u.a.d[ord[nNonzero]] 
	}else {
		strict = ans0 = tail(u.a.d,1L)
	}
	ans = ans0 / supp
	attr(ans, 'strict')=strict / supp
	ans
}


bw.kdep.optim.mrpp <-
function(y, start.bw = NULL, kernel='triweight', 
	bw.method='amse(z[1]).pearson3', scale='raw', verbose=TRUE, n.subset, mrpp.stats=NULL, r=seq_len(y$R))
## y mrpp object; r=dimension index; 
{
	on.exit({ # recover original data.env
			if(!is.null(y$data.env$y.bakBw)){
				y$data.env$y=y$data.env$y.bakBw
				rm(y.bakBw, envir=y$data.env)
			}
	})
	if(length(r)!=y$R){
		y$data.env$y.bakBw = y$y
		y[['data.env']]$y=y[['data.env']]$y[,r,drop=FALSE]
		y$R=length(r)
		y$distObj=y$distFunc(y$y)
		r=seq_len(y$R)
	}

	if(is.null(mrpp.stats)) {
		mrppt=mrpp.test(y, test.method='permutation'); 
		mrpp.stats=mrppt$all.statistics
	}	
	bw.method=match.arg(bw.method, choices=c('sym1','drop1','add1','keep1','amse(z[1]).WJ','amse(z[1]).pearson3','amise.WJ','amise.pearson3','pearson3','tlnorm', 'p3tlnormmix','pearson3gca')) # 'dropadd1','dropaddsym1','ss.gradp',
	scale=match.arg(scale, choices=c('raw','log','logit'))
    kernel=match.arg(kernel, .kernels)
	R=y$R
	if(R==1L ) bw.method='amse(z[1]).pearson3'

#	lst=list(...)
#	lst$y=y[,r,drop=FALSE]
#	lst$permutedTrt=permutedTrt
#	lst$distFunc=distFunc 
#	nms = setdiff(names(formals(mrpp)), '...')
#	idx=names(lst)%in%nms
#	mrpp.obj = do.call('mrpp', lst[idx])
#	mrpp=mrpp.test(mrpp.obj, test.method='permutation')
#	distObj=mrpp.obj$distObj
	
	bw.range=bw.range(mrpp.stats,kernel=kernel,length=10L,var.infl=0.1)
	log.bw.range=log(range(bw.range))
	
	lower.bound = bw.safety(mrpp.stats, kernel=kernel, nNonzero=3L, pdf.cut=1e-3)
	strict.lb = attr(lower.bound, 'strict')
	upper.bound = tail(bw.range, 1L)
	overshoot.msg='Bandwidth selected seems large; please check diagnostic plots.'
	undershoot.msg='Bandwidth selected is too small, replaced by a safer lower bound.'
	
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
	if(substr(bw.method,1L,5L)=='amise'){
		if( (tmp=substr(bw.method,7L,nchar(bw.method))) =='pearson3') {
			cums=cumulant(y, 2:3); cums[1L]=sqrt(cums[1L]); cums[2L]=cums[2L]/cums[1L]^3
			ans=bw.amise.pdf(mrpp.stats, kernel=kernel,verbose=verbose,deriv.order=0L, method=tmp, sd=cums[1L], skew=cums[2L])
		}else ans=bw.amise.pdf(mrpp.stats, kernel=kernel,verbose=verbose,deriv.order=0L, method=tmp)
		if(ans<=strict.lb ){
			if(verbose) warning(undershoot.msg)
			ans=lower.bound
		}else if(verbose && ans >= upper.bound){
			warning(overshoot.msg)
		}
		if(verbose) {
			dev.new(noRStudioGD=TRUE)
			eval(plot.expr)
			title(main='\n\nMethod: KDE AMISE formula',cex.main=.8)
		}
		return(ans)
	}
	if(substr(bw.method,1L,10L)=='amse(z[1])'){
		if( (tmp=substr(bw.method,12L,nchar(bw.method))) =='pearson3') {
			cums=cumulant(y, 2:3); cums[1L]=sqrt(cums[1L]); cums[2L]=cums[2L]/cums[1L]^3
			ans=bw.amse.pdf(mrpp.stats, kernel=kernel,verbose=verbose,method=tmp, sd=cums[1L], skew=cums[2L])
		}else ans=bw.amse.pdf(mrpp.stats, kernel=kernel,verbose=verbose,method=tmp)
		if(ans<=strict.lb){
			if(verbose) warning(undershoot.msg)
			ans=lower.bound
		}
		if(verbose) {
			dev.new(noRStudioGD=TRUE)
			eval(plot.expr)
			tmp=dkde(mrpp.stats, kernel=kernel,bw=attr(ans,'bw.f(x[1])'))
			curve(tmp(x), min(mrpp.stats)-marg, max(mrpp.stats)+marg, col=2,add=TRUE)
			tmp=dkde(mrpp.stats, kernel=kernel, bw=attr(ans,"bw.f''(x[1])"))
			curve(tmp(x), min(mrpp.stats)-marg, max(mrpp.stats)+marg, col=5,add=TRUE)
			legend('topleft', lty=1L, col=c(1:2,5), text.col=c(1:2,5), legend=c('final (AMSE) density', 'AMISE density', 'AMISE density for 2nd deriv.'))
			title(main='\n\nMethod: KDE AMSE formula',cex.main=.8)
		}
		return(ans)
	}

	weight.trt = y$weight.trt

	if(bw.method=='pearson3' || bw.method=='pearson3gca') {
		cums=cumulant(y,order=seq_len(.pdfnmoment[bw.method])); 
		cums[2L]=sqrt(cums[2L]); cums[3L]=cums[3L]/cums[2L]^3; 
		pdf0x = if(bw.method=='pearson3') dpearson3(mrpp.stats, cums[1L], cums[2L], cums[3L]) else 
		dpearson3gca(mrpp.stats, cums[1L], cums[2L], cums[3L], cums[4L]/cums[2L]^4) 
		ans=bw.matchpdf(mrpp.stats, kernel=kernel, pdf=pdf0x, bw=NULL, verbose=verbose, title=switch(bw.method, pearson3="Match Pearson III Dist'n", pearson3gca="Match GCA Adjust. of Pearson III Dist'n"))
		if(ans<=strict.lb){
			if(verbose) warning(undershoot.msg)
			ans=lower.bound
		}
		return(ans)
	}
	if(bw.method=='tlnorm') {
		cums=cumulant(y,order=seq_len(.pdfnmoment[bw.method])); 
		cums[2L]=sqrt(cums[2L]); cums[3L]=cums[3L]/cums[2L]^3; 
		pdf0x = dtlnorm(mrpp.stats, cums[1L], cums[2L], cums[3L]) 
		ans=bw.matchpdf(mrpp.stats, kernel=kernel, pdf=pdf0x, bw=NULL, verbose=verbose, title="Match Transformed Log-Normal Dist'n")
		if(ans<=strict.lb){
			if(verbose) warning(undershoot.msg)
			ans=lower.bound
		}
		return(ans)
	}
	if(bw.method=='p3tlnormmix') {
		cums=cumulant(y,order=seq_len(.pdfnmoment[bw.method])); 
		cums[2L]=sqrt(cums[2L]); cums[3L]=cums[3L]/cums[2L]^3; cums[4L]=cums[4L]/cums[2L]^4
		pdf0x = dMixP3Tln(mrpp.stats, cums[1L], cums[2L], cums[3L], cums[4L]) 
		ans=bw.matchpdf(mrpp.stats, kernel=kernel, pdf=pdf0x, bw=NULL, verbose=verbose, title="Match Mixture of Pearson III and Transformed Log-Normal Dist'n")
		if(ans<=strict.lb){
			if(verbose) warning(undershoot.msg)
			ans=lower.bound
		}
		return(ans)
	}

	r.bak=r
	if(bw.method%in%c('sym1','drop1','add1','keep1')) {#,'dropadd1','dropaddsym1','ss.gradp')) {
		if(missing(n.subset))
			n.subset= if(R>100L) round(100+.2*(R-100)) else R
		if(n.subset>R) n.subset=R
		if(n.subset<R) {
			imp0=grad.kdep.Inf.mrpp(y=y,test=FALSE,mrpp.stats=mrpp.stats)
			ord0=order(imp0)
			r=unique(sort(ord0[round(seq.int(from=1L, to=R, length.out=n.subset))]))

			if(is.null(y$data.env$y.bakBw)) y$data.env$y.bakBw = y$y
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
	
	

	dkern=dkernel(kernel)
	ddkern=ddkernel(kernel)
	adjust.string=switch(scale, raw='none', log='log scale', logit='logit scale')
	pval0=p.empirical(mrpp.stats, midp=FALSE) - as.numeric(0.5/y$nparts)
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
		ans=ans/switch(adjust.string, 'none'=1, 'log scale'=pval0, 'logit scale'=pval0*(1-pval0))
		ans=switch(adjust.string, 'none'=ans, 'log scale'=, 'logit scale'=exp(ans))
		attr(ans, 'parameters')=list(adjust=adjust.string)
		attr(ans, 'midp')=pval0
		class(ans)='grad.kdep'
	})

	if(FALSE && bw.method=='ss.gradp'){ 
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
		if(verbose && (idx==1L || idx==length(bw))) warning('optimal bandwidth occurs at the boundary')
		ans =  exp(mean(log(ans)))  
		if(ans<=strict.lb){
			if(verbose) warning(undershoot.msg)
			ans=lower.bound
		}
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
			mrpp.test.mrpp(mrpp.obj1, test.method='pearson3')$p.value
		} # >97% of time is spent in mrpp.test.mrpp
		drop1pval= sapply(r, drop1p)
	}
	if(any(bw.method==c('add1','sym1'))){
		mrpp.obj1$R=y$R+1L
		add1p = function(r.i){
			#mrpp.obj1$distObj=y$distFunc(y$data.env$y[,c(r,r.i),drop=FALSE])
			#mrpp.obj1$distObj=sqrt(distObj2+y$distFunc(y$data.env$y[,r.i,drop=FALSE])^2)
			mrpp.obj1$distObj=sqrt(distObj2+gradEnv$all.uni.dist2[,r.i])
			mrpp.test.mrpp(mrpp.obj1,test.method='pearson3')$p.value
		} 
		add1pval=sapply(r, add1p)
	}
	if(bw.method=='keep1'){
		mrpp.obj1$R=1L
		keep1 = function(r.i){
			mrpp.obj1$distObj[]=sqrt(gradEnv$all.uni.dist2[,r.i])
			mrpp.test.mrpp(mrpp.obj1,test.method='pearson3')$p.value
		}
		keep1pval=sapply(r, keep1)
		#t(smps - gradEnv$ans%*%(1-diag(1, length(r), length(r))))
	}
	
	sses.expr=switch(bw.method,
		drop1=quote(sum((p.value.grad.kdep(ans, type='drop1')-drop1pval)^2)),
		add1 =quote(sum((p.value.grad.kdep(ans, type='add1' )-add1pval)^2)),
		keep1=quote(sum((p.value.grad.kdep(ans, type='keep1')-keep1pval)^2)),
		sym1 =quote(sum((0.5*(p.value.grad.kdep(ans, type='add1' )
		                     -p.value.grad.kdep(ans, type='drop1'))
						 -.5*(add1pval - drop1pval)
						)^2))
	)
	objfunc=eval(bquote(function(logbw){
		.(obj.common.expr)
		thisEnv$sses[[niter]]=.(sses.expr)
	}))
	if(FALSE && scale=='raw'){
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
	}else if(FALSE && scale=='log'){
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
	} 
	

		

	# starting values when bw.method is not am[i]se, pearson3[gca]
	fixup.range=function(rg, mid){
		golden=(sqrt(5)-1)*.5; goldeni=1+golden
		if(missing(mid) || mid<=rg[1L]){
			mid=golden*(rg[1L]+rg[2L]*golden)
		}else if(mid>=rg[2L]){
			tmp.lo=mean(rg)
			rg=c(tmp.lo, mid*(2+golden)-goldeni*tmp.lo)
		}else if( (tmp.hi=mid*(2+golden)-goldeni*rg[1L]) > rg[2L] ){ #mid is close to upper limit
			rg=c(goldeni*mid - golden*rg[2L],  rg[2L])
		}else{ #mid is close to lower limit
			rg=c(rg[1L], tmp.hi)
		}
		return(list(range=rg, mid=mid))
	}
	if(is.null(start.bw)){
		#start.bw = bw.kdep.optim.mrpp(y, kernel=kernel, bw.method='amse(z[1])', mrpp.stats=mrpp.stats, verbose=FALSE)
		#start.bw = bw.kdep.optim.mrpp(y, kernel=kernel, bw.method='amise', mrpp.stats=mrpp.stats, verbose=FALSE)
		#range.start=fixup.range(log.bw.range, log(start.bw))
		start.bw=bw.range
	}
	n.bw=length(start.bw)
	if(n.bw==2L){
		range.start=fixup.range(log(sort(start.bw)))
	}else if(n.bw>2L){
		locs=sort(start.bw)
		# initial grid search
		gradEnv.grid=new.env(parent=gradEnv)
		evalq({
			ans=matrix(NA_real_, n.bw, length(r))
	
			bw0=rep(locs, each=B)
			weights=dkern(mrpp.stats.diff1/bw0)/bw0
			dim(weights) = c(B, n.bw)
			w.idx=which(apply(weights>0, 1L, any))
			w.pos=weights[w.idx,,drop=FALSE]
			B.sub=length(w.idx)
			for(r.i in r){
				dz.dw=.Call(mrppstats_subset, all.ddelta.dw[,r.i], y$permutedTrt, as.numeric(weight.trt), w.idx, PACKAGE='MRPP')
				ans[, r.i] = 1/B* .colSums(w.pos* (dz.dw[1L]-dz.dw), B.sub, n.bw)
			}
			ans=ans/switch(adjust.string, 'none'=1, 'log scale'=pval0, 'logit scale'=pval0*(1-pval0))
			ans=switch(adjust.string, 'none'=ans, 'log scale'=, 'logit scale'=exp(ans))
			attrs=list(parameters=list(adjust=adjust.string), midp=pval0,class='grad.kdep')
			ans=as.data.frame(t(ans))
			ans=lapply(ans, 'attributes<-', attrs)
			sses=switch(bw.method,
				drop1=.colSums((sapply(ans, p.value.grad.kdep, type='drop1')-drop1pval)^2,y$R,n.bw),
				 add1=.colSums((sapply(ans, p.value.grad.kdep, type= 'add1')- add1pval)^2,y$R,n.bw),
				keep1=.colSums((sapply(ans, p.value.grad.kdep, type='keep1')-keep1pval)^2,y$R,n.bw),
				 sym1=.colSums((.5*(sapply(ans, p.value.grad.kdep, type= 'add1')-
									sapply(ans, p.value.grad.kdep, type='drop1'))
							   -.5*(add1pval - drop1pval))^2, y$R, n.bw)
			)
			idx=which.min(sses)
			log.start=log(locs[idx])
			log.brack=log(locs[c(max(1L,idx-1L), min(n.bw,idx+1L))])
		}, envir=gradEnv.grid)
		range.start=fixup.range(gradEnv.grid$log.brack, gradEnv.grid$log.start)
		niter=n.bw
		sses=gradEnv.grid$sses
	}else range.start=fixup.range(log.bw.range, log(start.bw))

	opt.range=range.start$range; 
	log.start=range.start$mid

	if(FALSE){
		opt.rslt=suppressWarnings(optim(log.start, objfunc, method='Nelder-Mead'))
	}else if(FALSE){
		opt.rslt=nlminb(log.start, objfunc)
	}else {
		# first attempt
			opt.rslt=optimize(objfunc, opt.range)
		# second attempt
			if(any(abs(opt.rslt$minimum-opt.range)< 1e-3)){
				opt.rslt=optimize(objfunc, log.bw.range)
			}
		# third attempt (grid search, followed by optimize)
			if(any(abs(opt.rslt$minimum-log.bw.range[2L])< 1e-3)){
				tmp.bws=seq(log.bw.range[1], log.bw.range[2],length=25L)
				lapply(tmp.bws, objfunc)
				ord=order(locs)
				idx=which.min(sses[ord])
				if(idx>1 && idx<niter){
					opt.rslt=optimize(objfunc, exp(locs[ord[idx+c(-1L,1L)]]))
				}
			}
	}
	ord=order(locs)
	bw=locs[ord]
	sses=sses[ord]
	idx=which.min(sses)
	ans=bw[idx]
	if(verbose && ans >= exp(log.bw.range[2L]))warning(overshoot.msg)
	min.sse=sses[idx]

	if(ans<=strict.lb){
		if(verbose) warning(undershoot.msg)
		ans=lower.bound
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
	
	ans
	
}

bw.kdep=bw.kdep.optim.mrpp

bw.matchpdf=function(x, kernel=.kernels, pdf, bw = NULL, verbose=FALSE, title=NULL)
{## TODO: Minimize (approximate?) Hellinger Distance
	if(missing(bw) || is.null(bw)) bw = bw.range(x,kernel,length=7)
	if(is.function(pdf)) pdf=pdf(x)
	#x=sort.int(x, method='quick')
	#diffx=diff(x)
	bw=sort(bw)
	n.bw=length(bw)
	n=length(x)
	
	#sse=function(bw) sum((dkde(x, bw, kernel)(x)-pdf)^2)
	#ss = sapply(bw, sse)
	if(TRUE){ ## grid search: minimize sum of squared differences
		if(n.bw==1L) {
			bw=bw*10^seq(log10(.2), log10(5), length=7L)
			n.bw=7L
		}else if(n.bw==2L){
			bw=c(min(bw), median(bw), max(bw))
			n.bw=3L
		}
		niter=0L
		locs=numeric(0L)
		ss=numeric(0L)
		repeat{
			if(niter>500L) browser()
			niter=niter+n.bw
			pdfs=dkde(x, bw, kernel)
			locs=c(locs,bw)
			dens=sapply(pdfs, do.call, list(v=x))
			ss0=.colSums((dens-pdf)^2, n, n.bw)
			idx=which.min(ss0)
			ss=c(ss,ss0)
			if(idx==1L){# factor 5
				bw=c(exp(seq(log(bw[1L]/5),log(bw[1L]),length=n.bw-1L)), bw[2L])
				next
			}else if(idx==n.bw){ #factor 5: different than 3 to avoid infinite swabbing
				bw=c(bw[n.bw-1L], exp(seq(log(bw[n.bwL]),log(bw[n.bw]*5),length=n.bw-1L)))
				next
			}else break
		}
		ans=bw[idx]
		min.ss=ss0[idx]
	}#else{# golden sectioning
		#pdf0=dkde(x, median(bw), kernel)
		pdf0=dkde(x, ans, kernel)
		env=attr(pdf0, 'environment')
		this.env=environment()
		#niter=0L
		#locs=numeric(0L)
		#ss=numeric(0L)
		func=function(logbw){
			niter <<- niter+1L
			this.env$locs[[niter]] = thisbw = exp(logbw)
			zetalstar = env$ftkern(thisbw*env$sl) * env$Yl
			zetak=fft(zetalstar)*c(1,-1)
			kde.est=approx(env$tk0, pmax.int(0,Re(zetak)), x, ties='ordered')$y
			#rt.fg=sqrt(kde.est*pdf)
			#this.env$ss[[niter]]=1-sum(caTools::runmean(rt.fg,k=2,alg='fast',endrule='trim')*diffx)
			this.env$ss[[niter]] = sum((kde.est-pdf)^2)
		}
		opt.rg.const=c(.5*(1-sqrt(5)),1)*log(10)
#		if(length(bw)>1L){
#			log.range=log(range(bw))
#		}else{
#			log.range=log(bw)+opt.rg.const
#		}
		log.range=log(bw[idx+c(-1L,1L)])
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
#	}
	if(verbose){
		dev.new(width=10, height=5, noRStudioGD=TRUE)
		par(mfrow=1:2)
		plot(log10(bw), log10(ss), xlab='bandwidth', ylab='log10(SSE of p-value density)', type='o', axes=FALSE )
		if(!is.null(title))graphics::title(main=paste0('\nMethod: ', title)) else graphics::title(main='Matching Initial PDF')

		axis(1L, at = log10(ans), labels=sprintf('%.3g',ans), col='blue',col.ticks='blue', col.axis='blue', line=1)
		ats=axTicks(1)
		axis(1L, at=ats, labels=parse(text=paste0('10^',ats))  )

		abline(v=log10(ans), h=log10(min.ss), col='blue', lty=3L)

		axis(2L, at = log10(min.ss), labels=sprintf('%.3g',min.ss), col='blue',col.ticks='blue', col.axis='blue', line=1)
		ats=axTicks(2)
		axis(2L, at=ats, labels=parse(text=paste0('10^',ats))  )

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
	if(verbose && any( ans==range(bw) ) )warning('optimal bandwidth occurs at the boundary')
	ans
}

