bw.mse.pdf.asym=function(x,iter.max=1L,eps=1e-6,start.bw=bw.nrd, kernel='gaussian', verbose=FALSE)
{#require(ks)
    if(is.function(start.bw)) bw0=start.bw(x)
    if(is.numeric(start.bw)) bw0=start.bw
    if(is.character(start.bw)) bw0=call(start.bw, x)

    Rkern=density(x,bw0,from=x[1L],to=x[1L],n=1L,kernel=kernel, give.Rkern=TRUE)

	dkern=dkernel(kernel)
    n.iter=1L
    xdiff=x[1L]-x
    repeat{
#### these lines are using the ks:::drvkde and or ks::kdde functions and  are replaced by explicit calculations
#        f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
#        ddf=drvkde(x,2,bw0,se=FALSE)
#        ddf.1=approx(ddf$x.grid[[1]],ddf$est,xout=x[1])$y
            tmp=xdiff/bw0
            f.1=mean(dkern(tmp))/bw0
            ddf.1=mean(dkern(tmp)*(tmp*tmp-1))/bw0/bw0/bw0
        bw=(f.1/ddf.1/ddf.1*Rkern/length(x))^.2
        if(abs(bw-bw0)<eps || n.iter>=iter.max) return(bw)
        if(verbose && isTRUE(n.iter%%verbose==0))cat(bw,fill=TRUE)
        bw0=bw
        n.iter=n.iter+1
    }
}

bw.range=function(x, length=200, lower=.05, upper=.95, safety=100)
{
	uniqstat=sort(unique(x))
	lo = quantile(diff(uniqstat), lower)
	hi = diff(quantile(uniqstat, prob=c(lower, upper)))
	lf = log10(safety)
	10^seq(log10(lo)-lf, log10(hi)+lf, length=length)
}

bw.smoothp <-
function(y, permutedTrt, r=seq_len(NCOL(y)), bw = NULL, 
	distFunc=dist,  kernel='gaussian', 
	method=c('sym1','drop1','add1','keep1','dropadd1','dropaddsym1','ss.gradp','kde.mse1','match.pear3'), verbose=FALSE, ...)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
	method=match.arg(method)
    if(!is.matrix(y) && !is.data.frame(y)) y = as.matrix(y)
	R=NCOL(y)
	if(R==1L ) method='ss.gradp'

	distObj = distFunc(y[,r,drop=FALSE])
	lst=list(...)
	lst$y = distObj
	lst$permutedTrt=permutedTrt
	nms = setdiff(c(names(formals(mrpp.test.dist)), names(formals(mrpp))), '...')
	idx=names(lst)%in%nms
	mrpp = do.call('mrpp.test.dist', c(lst[idx],list(method='permutation')))
	if(method=='kde.mse1'){
		lst=list(...)
		lst$x=mrpp$all.statistics
		lst$kernel=kernel
		lst$verbose=verbose
		nms=names(lst)
		idx= nms=='' || nms%in%names(formals(bw.mse.pdf.asym))
		return(do.call('bw.mse.pdf.asym', lst[idx]))
	}

	weight.trt = attr(mrpp$parameter, 'weight.trt')
	
	if(missing(bw) || is.null(bw)) bw = bw.range(mrpp$all.statistics)	
	n.bw=length(bw)
	
	if(method=='match.pear3') {
		cums=cumulant(mrpp, order=1:3); cums[2L]=sqrt(cums[2L]); cums[3L]=cums[3L]/cums[2L]^3
		pear3pdf = dpearson3(mrpp$all.statistics, cums[1L], cums[2L], cums[3L])
		return(bw.matchpdf(mrpp$all.statistics, kernel=kernel, pdf=pear3pdf, bw=bw, verbose=verbose)
	}
	
	## pre-computing all.ddelta.dw in grad.smoothp
	gradEnv=new.env()
	gradEnv$mrpp.stats=mrpp$all.statistics
	gradEnv$scale=1
	eval(grad.smoothp.bw, envir=gradEnv)
	
	if(method=='ss.gradp'){
		ss = rowSums(gradEnv$ans^2)
		idx = which.max(ss)
		ans=bw[idx]
		if(verbose){
			plot(log10(bw), (ss), xlab='bandwidth', ylab='SS of p-value gradients', type='l', axes=FALSE )
			axis(2)
			axis(1, at = log10(ans), labels=sprintf('%.1g',ans), col='blue',col.ticks='blue', col.axis='blue')
			ats=axTicks(1)
			axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
			abline(v=log10(ans), h=(ss[idx]), col='blue', lty=3)
			box()
		}
		if(idx==1L || idx==length(bw))warning('optimal bandwidth occurs at the boundary')
		return(  exp(mean(log(ans)))  )
	}
	
	if(method%in%c('drop1','sym1','dropadd1', 'dropaddsym1')){
		drop1p = function(r.i){
			lst$y=distFunc(y[,r[-r.i],drop=FALSE])
			p.value(do.call('mrpp.test.dist', lst),type="midp")
		}
		drop1pval= sapply(seq_along(r), drop1p)
	}
	if(method%in%c('add1','sym1','dropadd1','dropaddsym1')){
		add1p = function(r.i){
			lst$y=distFunc(y[,c(r,r[r.i]),drop=FALSE])
			p.value(do.call('mrpp.test.dist', lst),type="midp")
		}
		add1pval=sapply(seq_along(r), add1p)
	}
	if(method%in%c('keep1')){
		keep1 = function(r.i){
			lst$y=distFunc(y[,r[r.i],drop=FALSE])
			p.value(do.call('mrpp.test.dist', lst),type="midp")
		}
		keep1pval=sapply(seq_along(r), keep1)
		#t(smps - gradEnv$ans%*%(1-diag(1, length(r), length(r))))
	}
	sses = switch(method, 
		drop1= .rowSums((gradEnv$ans - (p.value(mrpp, type='midp') - drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) ,
		add1 = .rowSums((gradEnv$ans - (add1pval - p.value(mrpp, type='midp'))[col(gradEnv$ans)])^2, n.bw, length(r)) , 
		sym1 = .rowSums((gradEnv$ans - .5*(add1pval- drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) , 
		keep1= .rowSums((gradEnv$ans%*%(1-diag(1, length(r), length(r)))-keep1pval[col(gradEnv$ans)])^2, n.bw, length(r)) , 
		dropadd1 = .rowSums((gradEnv$ans - (p.value(mrpp, type='midp') - drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) 
			+ .rowSums((gradEnv$ans - (add1pval - p.value(mrpp, type='midp'))[col(gradEnv$ans)])^2, n.bw, length(r)) ,
		dropaddsym1 = .rowSums((gradEnv$ans - (p.value(mrpp, type='midp') - drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) 
			+ .rowSums((gradEnv$ans - (add1pval - p.value(mrpp, type='midp'))[col(gradEnv$ans)])^2, n.bw, length(r))
			+ .rowSums((gradEnv$ans - .5*(add1pval- drop1pval)[col(gradEnv$ans)])^2, n.bw, length(r)) 
	)
	idx=which.min(sses)
	ans=bw[idx]
	min.sse=sses[idx]
	
	if(verbose){
		plot(log10(bw), sses, xlab='bandwidth', ylab='SS of p-value approx errors', type='l', main='', axes=FALSE)
		title(main=switch(method, drop1='Backward Difference', add1='Forward Difference', keep1='Keep 1 Variable', sym1='Central Difference', dropadd1='Total From Backward and Forward Difference', dropaddsym1='Total From Backward, Forward and Central Difference'))
		axis(3, at = log10(ans), labels=sprintf('%.1g',ans), col='red',col.ticks='red',col.axis='red')
		axis(2)
		ats=axTicks(1)
		axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
		abline(h=min.sse, col='red', lty=3); box()
		
		axis(1, at = log10(ans), labels=sprintf('%.1g',ans), col='blue',col.ticks='blue', col.axis='blue')
		abline(v=log10(ans), col='blue', lty=3)
	}
	ans
}

bw.matchpdf=function(x, kernel=.kernels, pdf, bw = NULL, verbose=FALSE)
{
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
		opt.rslt=optimize(func, log(range(bw)))
		ans=exp(opt.rslt$minimum)
		min.ss=opt.rslt$objective
		ord=order(locs)
		bw=locs[ord]
		ss=ss[ord]
	}
	if(verbose){
		plot(log10(bw), log10(ss), xlab='bandwidth', ylab='log10(SSE of p-value density)', type='l', axes=FALSE )
		axis(2)
		axis(1, at = log10(ans), labels=sprintf('%.1g',ans), col='blue',col.ticks='blue', col.axis='blue')
		ats=axTicks(1)
		axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
		abline(v=log10(ans), h=log10(min.ss), col='blue', lty=3)
		box()
	}
	if(any( ans==range(bw) ) )warning('optimal bandwidth occurs at the boundary')
	ans
}
formals(bw.matchpdf)$kernel=.kernels

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

	midp = midp.empirical(x, eps)
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
	
	xx=seq(min(x), max(x), length=1e5)
	midpKde=kde(x, ,bw.constr.midp(x), eval.points=xx)
	
}
