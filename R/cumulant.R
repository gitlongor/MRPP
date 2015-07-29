mean.mrpp=function(x, ...)
{
	cumulant(x, order=1L)
}

var=function(x,...)UseMethod('var')
var.default=local({
	ans=stats::var
	formals(ans)=c(formals(stats::var), alist(...=))
	ans
})
var.mrpp=function(x,...)
{
	cumulant(x, order=2L)
}

sd=function(x,...)UseMethod('sd')
sd.default=local({
	ans=stats::sd
	formals(ans)=c(formals(stats::sd), alist(...=))
	ans
})
sd.mrpp=function(x,...)
{
	call=match.call()
	call[[1L]]=as.name('var.mrpp')
	envir=parent.frame()
	sqrt(eval(call, envir=envir))
}

skewness=function(x,...)UseMethod('skewness')
skewness.default=local({
	ans=moments::skewness
	formals(ans)=c(formals(moments::skewness), alist(...=))
	ans
})
skewness.mrpp=function(x,...)
{
	cum=cumulant(x, order=2:3)
	cum[2L]/cum[1L]^1.5
}

kurtosis=function(x,excess=TRUE,...)UseMethod('kurtosis')
kurtosis.default=local({
	ans=moments::kurtosis
	fm =formals(moments::kurtosis)
	formals(ans)=c(fm, alist(excess=TRUE, ...=))[c('x','excess',setdiff(names(fm),c('x','...','excess')), '...')]
	ans
})
kurtosis.mrpp=function(x, excess=TRUE, ...)
{
	.NotYetImplemented()
}


moment=function(x, order, central,...)UseMethod('moment')
moment.default=function(x, order=1:3, central=FALSE, ...)
{
	order=as.integer(order)
	stopifnot(all(order>=0))
	mOrd=max(order)
	if(mOrd<2L){
		ans=rep(1,length(order))
		ans[order==1L]=if(central) 0 else mean(x)
		return(ans)
	}
	structure(moments::all.moments(x, order.max=mOrd, central=central,...)[order+1L], names=as.character(order))
}
moment.mrpp=function(x, order=1:3, central=FALSE,...)
{
	order=as.integer(order)
	mOrd=max(order)
	if(mOrd>=4) .NotYetImplemented()
	cum=cumulant(x, seq(mOrd))
	
	ans=numeric(length(order)); names(ans)=as.character(order)

	ans[order==0L]=1
	ans[order==1L]=if(central) 0 else cum[1L]
	ans[order==2L]=if(central) cum[2L] else cum[2L]+cum[1L]^2
	ans[order==3L]=if(central) cum[3L] else cum[3L] + 3*cum[1L]*cum[2L] + cum[1L]^3
	
	ans
}

cumulant = function(x, order, ...)UseMethod('cumulant')
cumulant.default=function(x, order=1:3, ...)
{
	order=as.integer(order)
	stopifnot(all(order>=0))
	mOrd=max(order)
	mu.raw = moments::all.moments(x, order.max=mOrd, central=FALSE, ...)
	ans = c(0, moment2cumulant(mu.raw[-1L]))
	structure(ans[order+1L], names=as.character(order))
}
cumulant.mrpp=function(x, order=1:3,...)
{
	order=as.integer(order)
	if(any(order)>=4) .NotYetImplemented()

	ans=seq_along(order)
	names(ans)=as.character(order)
	ans[order==0L]=0
	mOrd = max(order)
	
	distObj = x$distObj
	permutedTrt = x$permutedTrt
	weight.trt  = x$weight.trt
	n			= x$n
	N			= x$nobs
	
	Nc=function(C)factorial.rising(N-C+1L, C)
	nc=function(n, C)factorial.rising(n-C+1L, C)
	
	dmat = as.matrix(distObj)
	
	if(mOrd>=1L) {
		N2=Nc(2L)
		d1=2*sum(distObj)
		D1 =d1/N2
		ans[order==1L]=mu = as.numeric(D1)
	}
	
	if(mOrd>=2L){
		N3=Nc(3L)
		N4=Nc(4L)
		d2 = 2*sum(distObj^2)
		D2 = d2/N2
		
		d1I=rowSums(dmat)
		D2p=(sum(d1I^2)-d2)/N3
		D2pp=(d1^2-4*N3*D2p-2*d2)/N4
		
		w2=weight.trt^2
		n2=sapply(n, nc, C=2)
		ans[order==2L] = sig2 = 2 * ( sum( w2/n2 ) - as.numeric(1L/N2)) *
		as.numeric(D2 - 2 * D2p + D2pp) +
		4 * ( sum(w2 / n) - 1/N) * 	as.numeric(D2p - D2pp)
	}
	
	if(mOrd>=3L){
		d3=2*sum(distObj^3)
		D3=d3/N2
		N5=Nc(5L)
		N6=Nc(6L)
		
		d2I=rowSums(dmat^2)
		D3p=(sum(d1I*d2I)-d3)/N3
		D3pp=(d1*d2-4*N3*D3p-2*d3)/N4
		D3s=6*sum(combn(N,3L,function(idx)distObj[[ idx[1L], idx[2L] ]] * distObj[[ idx[1L], idx[3L] ]] * distObj[[ idx[2L], idx[3L] ]])) / N3
		D3ss=(sum(dmat*(d1I%o%d1I))-N3*(2*D3p+D3s)-d3)/N4
		D3s3=(sum(d1I^3)-3*N3*D3p-d3)/N4
		D3p3=(N3*(d1*D2p-4*D3p-2*D3s)-2*N4*(2*D3ss+D3s3))/N5
		D3p4=(N4*(d1*D2pp-4*D3pp-8*D3ss)-8*N5*D3p3)/N6
		
		n2.2=n2*n2
		n2.3=n2.2*n2
		n3=sapply(n, nc, C=3)
		n4=sapply(n, nc, C=4)
		n5=sapply(n, nc, C=5)
		n6=sapply(n, nc, C=6)
		w3=w2*weight.trt
		mean.delta3 = 
			4*D3*sum(w3/n2.2)+
			8*(3*D3p+D3s)*sum(w3*n3/n2.3)+
			8*(3*D3ss+D3s3)*sum(w3*n4/n2.3)+
			6*D3pp*sum(w2*(1-weight.trt+weight.trt*n4/n2.2)/n2)+
			12*D3p3*sum(w2*((1-weight.trt)*n3+weight.trt*n5/n2)/n2.2)+
			D3p4*sum(weight.trt*((1-weight.trt)*(1-2*weight.trt)+3*weight.trt*(1-weight.trt)*n4/n2.2+w2*n6/n2.3))
		ans[order==3L] = mean.delta3 - 3*mu*sig2-mu^3
	}
	ans
}
