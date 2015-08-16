pkernel=function(kernel= .kernels)
{## Optimized based on the fact that ^2 and * are faster than other ^ powers
	kernel=match.arg(kernel)
	switch(kernel,
	gaussian=pnorm,
	biweight = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		x02=x0^2
		ans[idx]=0.5 + x0 * ( 0.9375  - 0.625 * x02 + 0.1875 * x02 ^2)
		ans[x>=1]=1
		attributes(ans)=attributes(x)
		ans
	},
	triweight = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		x02 = x0 ^2; x03 = x02 * x0
		ans[idx]=0.5  -0.03125 * x0 * (35 * x02 - 21 * x02^2 + 5 * x03^2 -35) 
		ans[x>=1]=1
		attributes(ans)=attributes(x)
		ans
	},
	tricube = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		sig = sign(x0)
		x02 = x0^2; x03=x02*x0; x04=x02^2
		ans[idx]=0.00617284 * (81 + 140 * x0 - sig * 105 * x04 + 60 * x03*x04 - sig * 14 * x04 * x03^2)
		ans[x>=1]=1
		attributes(ans)=attributes(x)
		ans
	},
	logistic =plogis)
}
formals(pkernel)$kernel=.kernels

dkernel=function(kernel= .kernels)
{
	kernel=match.arg(kernel)
	switch(kernel,
	gaussian=dnorm,
	biweight = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=15/16 * (1 - x0 *x0)^2
		attributes(ans)=attributes(x)
		ans
	},
	triweight = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=35/32 * (1  - x0 *x0)^3
		attributes(ans)=attributes(x)
		ans
	},
	tricube = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=70/81 * (1  - abs(x0^3))^3
		attributes(ans)=attributes(x)
		ans
	},
	logistic =dlogis)
}
formals(dkernel)$kernel=.kernels

fourier.kernel=function(kernel= .kernels, root.2pi=TRUE)
{
	kernel=match.arg(kernel)
	root.2pi= if(root.2pi) sqrt(2*base::pi) else 1
	switch(kernel,
	gaussian=function(s) exp(-.5*s*s) /root.2pi,
	biweight = function(s){
		ss=sin(s); s2=s*s; s4=s2*s2
		ans=(-15*(3*s*cos(s) - 3*ss + s2*ss))/(s4*s) /root.2pi
		idx=which(abs(s)<5e-2) ## close to 0/0 region: taylor series
		ans[idx]=(1 - s2[idx]/14 + s4[idx]/504 - s4[idx]*s2[idx]/33264) / root.2pi
		ans
	},
	triweight = function(s){
		cs=cos(s); ss=sin(s);
		s2=s*s; s4=s2*s2; s6=s4*s2; 
		ans=(105*(-15*s*cs + s2*s*cs + 15*ss - 6*s2*ss))/(s6*s) /root.2pi
		idx=which(abs(s)<1e-1) ## close to 0/0 region
		ans[idx]=(1 - s2[idx]/18 + s4[idx]/792 - s6[idx]/61776) / root.2pi
		ans
	},
	tricube = function(s){
		cs=cos(s); ss=sin(s);
		s2=s*s; s4=s2*s2; s6=s4*s2
		ans=(280*(20160 - s6 + 9*(-2240 + 1120*s2 - 80*s4 + s6)*cs - 
			36*s*(560 - 90*s2 + 3*s4)*ss))/(9*s4*s6) / root.2pi
		idx=which(abs(s)<3.5e-1) ## close to 0/0 region
		ans[idx]=(1 - (35*s2[idx])/486 + s4[idx]/528 - s6[idx]/37440 + s4[idx]^2/4199040) / root.2pi
		ans
   },
	logistic =function(s){
		ps=base::pi*s
		ans=ps/sinh(ps) /root.2pi
		idx=which(abs(s)<1e-2) ## close to 0/0 region
		ps2=ps*ps; ps4=ps2*ps2
		ans[idx]=(1 - ps2[idx]/6 + (7*ps4[idx])/360 - (31*ps4[idx]*ps2[idx])/15120) / root.2pi
		ans
	}
	)
}
formals(fourier.kernel)$kernel=.kernels

pkde=function(x, bw=bw.nrd, kernel=.kernels)
{
    if(is.character(bw)) bw=get(bw, mode='function')
    if(is.function(bw))  bw=bw(x)
    stopifnot(is.numeric(bw))
	if(is.character(kernel)) {
		Kernel = pkernel(kernel)
	}else if(is.function(kernel)) Kernel = kernel
	stopifnot(Kernel(-Inf)==0 && Kernel(Inf)==1)

	knots=x; n=length(knots); rm(x)
	function(x)
	{	nx=length(x)
		xdiff=x-rep(knots, each=nx)
		dim(xdiff) = c(nx, n)
		ans=.rowMeans(Kernel(xdiff/bw), nx, n)
		attributes(ans)=attributes(x)
		ans
	}
}
formals(pkde)$kernel=.kernels

if(FALSE) {
dkde=function(x, bw=bw.nrd, kernel=.kernels)
{## naive but should always work
    if(is.character(bw)) bw=get(bw, mode='function')
    if(is.function(bw))  bw=bw(x)
    stopifnot(is.numeric(bw))
	if(is.character(kernel)) {
		Kernel = dkernel(kernel)
	}else if(is.function(kernel)) Kernel = kernel

	knots=x; n=length(knots); rm(x)
	function(x)
	{	nx=length(x)
		xdiff=x-rep(knots, each=nx)
		dim(xdiff) = c(nx, n)
		ans=.rowMeans(Kernel(xdiff/bw)/bw, nx, n)
		attributes(ans)=attributes(x)
		ans
	}
}
dkde=function(x, bw=bw.nrd, kernel=.kernels)
{## not sure what R does when kernel is not gaussian
	fit=density(x, bw=bw, kernel=kernel, from=min(x)-3*bw, to=max(x)+3*bw, n=nextn(length(x)*2L,2L))
	approxfun(fit$x, fit$y)
}
}

dkde=function(x, bw=bw.nrd, kernel=.kernels, from=min(x)-3*median(bw), to=max(x)+3*median(bw))
{# implementation based on FFT in Silverman 1984
    if(is.character(bw)) bw=get(bw, mode='function')
    if(is.function(bw))  bw=bw(x)
	n=length(x)
	from0=from; force(to)
	x=x-from; to=to-from; from=0
	M=nextn(n*2L, 2L)
	delta=(to-from)/M
	tk=seq(from=from, to=to, length=M+1L)
	#cuts=as.integer(cut(x, tk))
	cuts=findInterval(x, tk)
		buf=integer(M)
	ucut=.Call(radixSort_preallocMax, unique(cuts), buf, M)
	xik=numeric(M+1L)
		fcuts=factor(cuts,levels=ucut)
		#attr(fcuts, 'class')='ordered'
		#attr(fcuts,'levels')=ucut
		xik[ucut]= sapply(split((tk[cuts+1L]-x)/n/delta^2, fcuts), sum)
		xik[ucut+1L]=xik[ucut+1L] + sapply(split((x-tk[cuts])/n/delta^2, fcuts), sum)
	xik=xik[1:M]
	kk=seq_len(M)-1L
	l=-M/2L+kk
	
	#Yl.bak=colMeans(xik*exp(1i*2*pi*(kk %o% l)/M))
	#Yl=fft(xik/exp(1i*pi*kk), inverse = TRUE)/M
	Yl=fft(xik*c(1,-1), inverse = TRUE)/M
	
	sl=2*pi*l/(to-from)
	ftkern=fourier.kernel(kernel, root.2pi=FALSE)
	
	tk0=tk[1:M]+from0; 
	ans=vector('list', length(bw))
	for(i in seq_along(bw)){
		zetalstar=ftkern(bw[i]*sl) * Yl
		#zetak.bak = colSums(exp((l%o%kk)/M*(1i)*-2*pi)*zetalstar) ## correct but slow
		#zetak=fft(zetalstar)*exp(1i*pi*kk)
		zetak=fft(zetalstar)*c(1,-1)
		ans[[i]] = approxfun(tk0, pmax.int(0,Re(zetak)), ties='ordered')
	}
	structure(if(length(ans)==1L) ans[[1L]] else ans, 
		environment=environment()
	)
}
formals(dkde)$kernel=.kernels
