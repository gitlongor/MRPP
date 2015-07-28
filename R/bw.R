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
	method=c('drop1','keep1','dropkeep1','ss.gradp','kde.mse1'), verbose=FALSE, ...)
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
	nms = setdiff(names(formals(mrpp.test.dist)), '...')
	idx=names(lst)%in%nms
	mrpp = do.call('mrpp.test.dist', lst[idx])
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
	
	## pre-computing all.ddelta.dw in grad.smoothp
	gradEnv=new.env()
	gradEnv$mrpp.stats=mrpp$all.statistics
	gradEnv$scale=1
	eval(grad.smoothp.bw, envir=gradEnv)
	
	drop=switch(method, drop1 = TRUE, keep1 = FALSE, dropkeep1 = c(FALSE, TRUE), NA)
	if(!any(is.na(drop))) {
		smps = sapply(bw, function(bw) pkde(mrpp$all.statistics, bw=bw, kernel=kernel)(mrpp$all.statistics[1L]) )
		ans=numeric(length(drop))
		sses = matrix(NA_real_, length(drop), length(bw))
	
		for(d.i in seq_along(drop)) {
			if(drop[d.i]){
				permp = function(rr){
					lst$y=distFunc(y[,-rr,drop=FALSE])
					p.value(do.call('mrpp.test.dist', lst),type="midp")
				}
				approxp = t(smps - gradEnv$ans)
			}else{
				permp = function(rr){
					lst$y=distFunc(y[,rr,drop=FALSE])
					p.value(do.call('mrpp.test.dist', lst),type="midp")
				}
				approxp = t(smps - gradEnv$ans%*%(1-diag(1, length(r), length(r))))
			}
			unips=sapply(r, permp)
			sses[d.i,] = colSums((approxp - unips)^2)
			idx=which.min(sses[d.i,])
			ans[d.i]=bw[idx]

			if(idx==1L || idx==length(bw))warning('optimal bandwidth occurs at the boundary')
		}
		if(verbose){
			min.sse=apply(sses, 1L, min)
			plot(log10(bw), sses[1L,], xlab='bandwidth', ylab='SS of p-value approx errors', type='l', ylim=c(min(min.sse), min(max(sses),10*max(min.sse))), main='', axes=FALSE)
			title(main=switch(method, drop1='Drop 1 Variable', keep1='Keep 1 Variable', dropkeep1='Drop 1 and/or Keep 1 Variable'))
			axis(3, at = log10(ans[1L]), labels=sprintf('%.1g',ans[1L]), col='red',col.ticks='red',col.axis='red')
			for(d.i in safeseq(2L, length(drop), by=1L)){
				lines(log10(bw), sses[d.i,])
				axis(3, at = log10(ans[d.i]), labels=sprintf('%.1g',ans[d.i]), col='red',col.ticks='red',col.axis='red')
			}
			axis(2)
			ats=axTicks(1)
			axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
			abline(v=log10(ans), h=min.sse, col='red', lty=3); box()
			
			ans=exp(mean(log(ans)))
			axis(1, at = log10(ans), labels=sprintf('%.1g',ans), col='blue',col.ticks='blue', col.axis='blue')
			abline(v=log10(ans), col='blue', lty=3)
		}else ans=exp(mean(log(ans)))
	}else if(method=='ss.gradp'){
		ss = rowSums(gradEnv$ans^2)
		idx = which.max(ss)
		ans=bw[idx]
		if(verbose){
			plot(log10(bw), ss, xlab='bandwidth', ylab='SS of p-value gradients', type='l', axes=FALSE )
			axis(2)
			axis(1, at = log10(ans), labels=sprintf('%.1g',ans), col='blue',col.ticks='blue', col.axis='blue')
			ats=axTicks(1)
			axis(1, at=ats, labels=parse(text=paste0('10^',ats))  )
			abline(v=log10(ans), h=ss[idx], col='blue', lty=3)
			box()
		}
		if(idx==1L || idx==length(bw))warning('optimal bandwidth occurs at the boundary')
		ans = exp(mean(log(ans)))
	}else stop('unknown method argument')
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

	midp = midp(x, eps)
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
