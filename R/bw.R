midp = function(x, eps=1e-8)
{
	mean(x<=x[1L]-eps) + .5* mean(x>x[1L]-eps & x<x[1L]+eps)
}
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

pkernel=function(kernel= c("gaussian", "biweight", 'triweight', 'tricube',"logistic"), eps=1e-8)
{
	kernel=match.arg(kernel)
	switch(kernel,
	gaussian=pnorm,
	biweight = function(x){
		ans=numeric(0, length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=1/16 * x0 *(15 - 10 * x0 *x0  + 3 *x0^4)
		ans[x>=1]=1
		ans
	},
	triweight = function(x){
		ans=numeric(0, length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=1/32 * (35 * x0 - 35 * x0^3 + 21 * x0^5 - 5 * x0^7)
		ans[x>=1]=1
		ans
	},
	tricube = function(x){
		ans=numeric(0, length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=1/162 *(81 + 140 *x0 -sign(x0)* 105 *x0^4 + 60 *x0^7 -sign(x0)* 14 *x0^10)
		ans[x>=1]=1
		ans
	},
	logistic =plogis)
}

dkernel=function(kernel= c("gaussian", "biweight", 'triweight', 'tricube',"logistic"), eps=1e-8)
{
	kernel=match.arg(kernel)
	switch(kernel,
	gaussian=dnorm,
	biweight = function(x){
		ans=numeric(0, length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=15/16 * (1 - x0 *x0)^2
		ans
	},
	triweight = function(x){
		ans=numeric(0, length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=35/32 * (1  - x0 *x0)^3
		ans
	},
	tricube = function(x){
		ans=numeric(0, length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=70/81 * (1  - abs(x0^3))^3
		ans
	},
	logistic =dlogis)
}

bw.mse.pdf.asym=function(x,iter.max=1L,eps=1e-6,start.bw=bw.nrd, verbose=FALSE)
{#require(ks)
    if(is.function(start.bw)) bw0=start.bw(x)
    if(is.numeric(start.bw)) bw0=start.bw
    if(is.character(start.bw)) bw0=call(start.bw, x)

    Rkern=density(x,bw0,from=x[1L],to=x[1L],n=1L,give.Rkern=TRUE)

    n.iter=1L
    xdiff=x[1L]-x
    repeat{
#### these lines are using the ks:::drvkde and or ks::kdde functions and  are replaced by explicit calculations
#        f.1=density(x,bw0,from=x[1],to=x[1],n=1)$y
#        ddf=drvkde(x,2,bw0,se=FALSE)
#        ddf.1=approx(ddf$x.grid[[1]],ddf$est,xout=x[1])$y
            tmp=xdiff/bw0
            f.1=mean(dnorm(tmp))/bw0
            ddf.1=mean(dnorm(tmp)*(tmp*tmp-1))/bw0/bw0/bw0
        bw=(f.1/ddf.1/ddf.1*Rkern/length(x))^.2
        if(abs(bw-bw0)<eps || n.iter>=iter.max) return(bw)
        if(verbose && isTRUE(n.iter%%verbose==0))cat(bw,fill=TRUE)
        bw0=bw
        n.iter=n.iter+1
    }
}

bw.grad.smoothp <-
function(y, permutedTrt, r=seq_len(NCOL(y)), bw, 
	test=FALSE, distFunc=dist,  kernel='gaussian', wtmethod=0, 
	drop=FALSE, verbose=FALSE, ...)
## y=N-by-p data matrix; b=permutation index for the 1st trt; r=dimension index; 
{
	permp = if(drop){
		function(rr)mrpp.test.dist(distFunc(y[,-rr,drop=FALSE]), permutedTrt=permutedTrt, wtmethod=0, ...)$midp
	}else{
		function(rr)mrpp.test.dist(distFunc(y[,rr,drop=FALSE]), permutedTrt=permutedTrt, wtmethod=0, ...)$midp
	}
    if(!is.matrix(y) && !is.data.frame(y)) y = as.matrix(y)
	unips=sapply(r, permp)
	
	mrpp = mrpp.test(y[,r,drop=FALSE], permutedTrt=permutedTrt, wtmethod=0, ...)
	distObj = distFunc(y[,r,drop=FALSE])
	
	if(missing(bw)) {
		uniqstat=unique(mrpp$all.statistics)
		lo = quantile(diff(sort(uniqstat)), .05)
		hi = diff(quantile(uniqstat, prob=c(.05, .95)))
		fact=100; lf = log10(fact)
		bw = 10^seq(log10(lo)-lf, log10(hi)+lf, length=100)
	}
	
	sig = if(drop) 1 else 1-diag(1, length(r), length(r))
	unip.approx=function(bw)
	{
		smp = smoothp(mrpp$all.statistics, bw=bw, kernel=kernel)$p.value
		grad = grad.smoothp(y, permutedTrt=permutedTrt,r=r,kernel=kernel, wtmethod=wtmethod, bw=bw, distObj=distObj, mrpp.stats=mrpp$all.statistics)
		smp - base::drop(sig %*% grad)
	}
	approxp = sapply(bw, unip.approx)
	sse = colSums((approxp - unips)^2)
	idx=which.min(sse)
	ans=bw[idx]
	if(verbose){
		plot(log10(bw), sse, xlab='log10(bandwidth)', ylab='SSE', type='l', ylim=c(0, min(max(sse),10*sse[idx])))
		abline(v=log10(ans), h=sse[idx], col='gray', lty=3)
	}
	if(idx==1L || idx==length(bw))warning('optimal bandwidth occurs at the boundary')
	ans
}


if(FALSE){
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
