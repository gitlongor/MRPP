pkernel=function(kernel= .kernels)
{## Optimized based on the fact that ^2 and * are faster than other ^ powers
	kernel=match.arg(kernel)
	switch(kernel,
	gaussian=pnorm,
	cosine = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		ans[idx]=0.5 * ( 1 + sin(1.5707963267948966 *x[idx]))
		ans[x>=1]=1
		attributes(ans)=attributes(x)
		ans
	},
	uniform = , rectangular = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		ans[idx]=0.5 * ( 1 + x[idx])
		ans[x>=1]=1
		attributes(ans)=attributes(x)
		ans
	},
	triangular = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]; ans[idx] = tmp=.5*(1-abs(x0))^2
		idx0=which(x0>0)
		ans[idx[idx0]]=1-tmp[idx0]
		ans[x>=1]=1
		attributes(ans)=attributes(x)
		ans
	},
	epanechnikov = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=0.25 * ( 2 + 3*x0 - x0*x0*x0)
		ans[x>=1]=1
		attributes(ans)=attributes(x)
		ans
	},
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
#formals(pkernel)$kernel=.kernels

dkernel=function(kernel= .kernels)
{
	kernel=match.arg(kernel)
	switch(kernel,
	gaussian=dnorm,
	cosine = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=0.78539816339744828*cos(1.5707963267948966*x0)
		attributes(ans)=attributes(x)
		ans
	},
	uniform = , rectangular = function(x){
		ans=ifelse(abs(x)<1, .5, 0)
		attributes(ans)=attributes(x)
		ans
	},
	triangular = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=1 - abs(x0)
		attributes(ans)=attributes(x)
		ans
	},
	epanechnikov = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=3/4 * (1 - x0 *x0)^2
		attributes(ans)=attributes(x)
		ans		
	},
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
#formals(dkernel)$kernel=.kernels

fourier.kernel=local({
	iroot.2pi=one12=one14=one18=one280=one33264=
	one360=one37440=one4199040=one504=one528=one6=
	one61776=one792=pipi=seven360=three1d15120=
	three5d486=two80d9=NULL
function(kernel= .kernels, root.2pi=TRUE)
{
	kernel=match.arg(kernel)
	enclEnv=new.env(hash=TRUE, size=37L)
	enclEnv$iroot.2pi= if(root.2pi) 1/sqrt(2*base::pi) else 1
	enclEnv$pipi=base::pi^2
	enclEnv$one12=1/12; enclEnv$one360=1/360
	enclEnv$one280=1/280
	enclEnv$one14=1/14; enclEnv$one504=1/504; enclEnv$one33264=1/33264
	enclEnv$one18=1/18; enclEnv$one792=1/792; enclEnv$one61776=1/61776
	enclEnv$three5d486=35/486; enclEnv$one528=1/528; enclEnv$one37440=1/37440; enclEnv$one4199040=1/4199040
	enclEnv$one6=1/6; enclEnv$seven360=7/360; enclEnv$three1d15120=31/15120
	enclEnv$two80d9=280/9
	ans=switch(kernel,
	gaussian=function(s) exp(-.5*s*s) *iroot.2pi,
	cosine=function(s)	pipi*cos(s)/(pipi-4*s*s) * iroot.2pi, 
	uniform = , rectangular = function(s){
		sinc(s) * iroot.2pi
	},
	triangular = function(s){
		s2=s*s; 
		ans=2*(1-cos(s))/s2 * iroot.2pi
		idx=which(abs(s)<1e-2) ## close to 0/0 region: taylor series
		ans[idx]=(1 - s2[idx]*one12 + s2[idx]*s2[idx]*one360) * iroot.2pi
		ans	
	},
	epanechnikov = function(s){
		s2=s*s; 
		ans=3*(sin(s)-s*cos(s))/(s2*s) * iroot.2pi
		idx=which(abs(s)<1e-2) ## close to 0/0 region: taylor series
		ans[idx]=(1 - s2[idx]*0.1 + s2[idx]*s2[idx]*one280) * iroot.2pi
		ans
	},
	biweight = function(s){
		ss=sin(s); s2=s*s; s4=s2*s2
		ans=((45*(ss - s*cos(s)) - 15*s2*ss))/(s4*s) * iroot.2pi
		idx=which(abs(s)<5e-2) ## close to 0/0 region: taylor series
		ans[idx]=(1 - s2[idx]*one14 + s4[idx]*one504 - s4[idx]*s2[idx]*one33264) * iroot.2pi
		ans
	},
	triweight = function(s){
		cs=cos(s); ss=sin(s);
		s2=s*s; s4=s2*s2; s6=s4*s2; 
		ans=(105*(s2*s*cs + 15*(ss -s*cs) - 6*s2*ss))/(s6*s) *iroot.2pi
		idx=which(abs(s)<1e-1) ## close to 0/0 region
		ans[idx]=(1 - s2[idx]*one18 + s4[idx]*one792 - s6[idx]*one61776) * iroot.2pi
		ans
	},
	tricube = function(s){
		cs=cos(s); ss=sin(s);
		s2=s*s; s4=s2*s2; s6=s4*s2
		ans=(two80d9*(20160 - s6 - 9*(2240 - 1120*s2 + 80*s4 - s6)*cs - 
			s*(20160 - 3240*s2 + 108*s4)*ss))/(s4*s6) *iroot.2pi
		idx=which(abs(s)<3.5e-1) ## close to 0/0 region
		ans[idx]=(1 - s2[idx]*three5d486 + s4[idx]*one528 - s6[idx]*one37440 + s4[idx]^2*one4199040) *iroot.2pi
		ans
   },
	logistic =function(s){
		ps=base::pi*s
		ans=ps/sinh(ps) *iroot.2pi
		idx=which(abs(s)<1e-2) ## close to 0/0 region
		ps2=ps*ps; ps4=ps2*ps2
		ans[idx]=(1 - ps2[idx]*one6 + ps4[idx]*seven360 - ps4[idx]*ps2[idx]*three1d15120) *iroot.2pi
		ans
	}
	)
	environment(ans)=enclEnv
	ans
}
#formals(fourier.kernel)$kernel=.kernels
})

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
#formals(pkde)$kernel=.kernels

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
#formals(dkde)$kernel=.kernels
