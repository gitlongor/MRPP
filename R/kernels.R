skernel=eval(substitute(function(kernel= .kernels, eps=1e-3)
{# (effective) support of kernel dist'n
	kernel=match.arg(kernel)
	if(eps==0){
		switch(kernel,
		cosine = ,
		uniform = , rectangular = ,
		triangular = ,
		epanechnikov =, parabolic=,
		biweight = , quartic=,
		triweight = ,
		tricube = c(-1,1),
		gaussian=,normal=,
		logistic =,
		sech=c(-Inf,Inf))
	}else{
		k=1/eps
		fix.root=function(rt, im.eps=1e-5){
			im=Im(rt);re=Re(rt); 
			idx=abs(im)<=im.eps & re<=1 & re>0; 
			max(re[idx])
		}

		c(-1,1)*switch(kernel, 
			cosine=2*acos(fourdpi/k)*onedpi,
			uniform=, rectangular=1,
			triangular=(k-1)/k,
			epanechnikov=,parabolic=sqrt(1-fourThirds/k),
			biweight=, quartic=fix.root(polyroot(c(-16+ 15 *k, 0,- 30* k,0, 15* k))),
			triweight=fix.root(polyroot(c(32 - 35* k,0,105 *k ,0, - 105 *k,0 ,35* k))),
			tricube=fix.root(polyroot(c(81 - 70 *k,0,0, 210 *k,0,0, - 210 *k ,0,0, 70* k))),
			gaussian=,normal=sqrt(2*log(k*iroot.2pi)),
			logistic=log(.5* k + .5*sqrt((-4 + k)* k)-1),
			sech=log(k + sqrt(k*k - pipi))-logpi
		)
	}
} #of function
,constEnv ) # of substitute
) # of eval 
attr(skernel, 'srcref')=NULL


pkernel = eval(substitute(
function(kernel= .kernels)
{# cdf of kernel dist'n
## Optimized based on the fact that ^2 and * are faster than other ^ powers
	kernel=match.arg(kernel)
	ans=switch(kernel,
	gaussian=,normal=pnorm,
	cosine = function(q){
		ans=numeric( length(q))
		idx = which(abs(q)<1)
		ans[idx]=0.5 * ( 1 + sin(pid2 *q[idx]))
		ans[q>=1]=1
		attributes(ans)=attributes(q)
		ans
	},
	uniform = , rectangular = function(q){
		ans=numeric( length(q))
		idx = which(abs(q)<1)
		ans[idx]=0.5 * ( 1 + q[idx])
		ans[q>=1]=1
		attributes(ans)=attributes(q)
		ans
	},
	triangular = function(q){
		ans=numeric( length(q))
		idx = which(abs(q)<1)
		x0=q[idx]; ans[idx] = tmp=.5*(1-abs(x0))^2
		idx0=which(x0>0)
		ans[idx[idx0]]=1-tmp[idx0]
		ans[q>=1]=1
		attributes(ans)=attributes(q)
		ans
	},
	epanechnikov =, parabolic= function(q){
		ans=numeric( length(q))
		idx = which(abs(q)<1)
		x0=q[idx]
		ans[idx]=0.5 + x0* ( 0.75 - 0.25 *x0*x0)
		ans[q>=1]=1
		attributes(ans)=attributes(q)
		ans
	},
	biweight = , quartic=function(q){
		ans=numeric( length(q))
		idx = which(abs(q)<1)
		x0=q[idx]
		x02=x0^2
		ans[idx]=0.5 + x0*(0.9375 + x02*(0.1875*x02-0.625))
		ans[q>=1]=1
		attributes(ans)=attributes(q)
		ans
	},
	triweight = function(q){
		ans=numeric( length(q))
		idx = which(abs(q)<1)
		x0=q[idx]
		x02 = x0 ^2; #x03 = x02 * x0
		ans[idx]=0.5 + x0*(thr5d32 + x02*(x02*(two1d32 - five32*x02) -thr5d32))
		ans[q>=1]=1
		attributes(ans)=attributes(q)
		ans
	},
	tricube = function(q){
		ans=numeric( length(q))
		idx = which(abs(q)<1)
		x0=q[idx]
		sig = sign(x0)
		x03=x0*x0*x0; 
		ans[idx]=0.5+x0*(sev0d81 + x03 * (x03*(one0d27-sev81*sig*x03) -sig*thr5d54))
		ans[q>=1]=1
		attributes(ans)=attributes(q)
		ans
	},
	logistic =plogis,
	sech = function(q)twodpi*atan(exp(q))
	) # of switch
	attr(ans ,'srcref')=NULL
	if(identical(environment(ans),sys.frame(sys.nframe()))) environment(ans)=pkgEnv
	ans
},  # of function
constEnv) # of substitute
) # of eval
attr(pkernel, 'srcref')=NULL

dkernel=eval(substitute(function(kernel= .kernels)
{# pdf of kernel dist'n
	kernel=match.arg(kernel)
	ans=switch(kernel,
	gaussian=,normal=dnorm,
	cosine = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=pid4*cos(pid2*x0)
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
	epanechnikov =, parabolic= function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=0.75 * (1 - x0*x0)
		attributes(ans)=attributes(x)
		ans		
	},
	biweight = , quartic=function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x0=x[idx]
		ans[idx]=0.9375 * (1 - x0 *x0)^2
		attributes(ans)=attributes(x)
		ans
	},
	triweight = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		x02=x[idx]^2
		ans[idx]=1.09375 + x02*(x02*(3.28125 - 1.09375*x02)-3.28125)
		attributes(ans)=attributes(x)
		ans
	},
	tricube = function(x){
		ans=numeric( length(x))
		idx = which(abs(x)<1)
		ax0=abs(x[idx]); tmp=1-ax0*ax0*ax0
		ans[idx]=sev0d81*tmp*tmp*tmp
		attributes(ans)=attributes(x)
		ans
	},
	logistic =dlogis, 
	sech = function(x)#twodpi/(exp(x)+exp(-x))
		onedpi/cosh(x)
	) # of switch
	attr(ans, 'srcref')=NULL
	if(identical(environment(ans),sys.frame(sys.nframe()))) environment(ans)=pkgEnv
	ans	
}
,constEnv) # of substitute
) # of eval
attr(dkernel, 'srcref')=NULL

krRkernel=eval(substitute(function(kernel=.kernels)
{# Parzen (1962) constants used in asym. optimal bandwidth
	kernel=match.arg(kernel)
	
	switch(kernel,
	#gaussian=c(kr=.5, r=2,R=.5/sqrt(base::pi)),
	gaussian=,normal=c(kr=.5, r=2,R=halfIrootPi),
	#cosine =c(kr=.5-4/base::pi^2, r=2, R=base::pi^2/16),
	cosine =c(kr=halfm4dpipi, r=2, R=pipid16),
	uniform = , rectangular = c(kr=oneSixth, r=2, R=.5),
	triangular=c(kr=one12, r=2, R=twoThirds),
	epanechnikov=,parabolic=c(kr=.1,r=2,R=0.6),
	#biweight=c(kr=1/14,r=2,R=5/7),
	biweight=, quartic=c(kr=one14,r=2,R=five7),
	#triweight=c(kr=1/18,r=2,R=350/429),
	triweight=c(kr=one18,r=2,R=three50d429),
	#tricube=c(kr=35/486,r=2,R=175/247),
	tricube=c(kr=three5d486,r=2,R=one75d247),
	#logistic=c(kr=base::pi^2/6,r=2,R=1/6),
	logistic=c(kr=pipid6,r=2,R=oneSixth),
	#sech=c(kr=base::pi^2/8,r=2,R=2/base::pi^2)
	sech=c(kr=pipid8,r=2,R=twodpipi)
	)# of switch
}
, constEnv) # of substitute
) # of eval
attr(krRkernel,'srcref')=NULL

ddkernel=eval(substitute(function(kernel=.kernels)
{
	kernel=match.arg(kernel)
	ans=switch(kernel,
		gaussian=,normal=function(x){
			x2=x*x; #iroot.2pi=1/sqrt(2*base::pi)
			-x*exp(-.5*x2) * iroot.2pi
		},
		cosine=function(x){
			ans=numeric(length(x))
			ax=abs(x)
			idx=which(ax<1)
			ans[idx]=-pipid8*sin(pid2*x[idx])
			ans[ax==1]=NaN
			attributes(ans)=attributes(x)
			ans
		},
		uniform=, rectangular=function(x){
			ans=numeric(length(x))
			ans[abs(x)==1]=NaN
			attributes(ans)=attributes(x)
			ans
		}, 
		triangular=function(x){
			ans=numeric(length(x))
			ax=abs(x); idx=which(ax<1)
			ans[idx]=-sign(x[idx])
			ans[ax==1]=NaN
			ans[x==0]=NaN
			attributes(ans)=attributes(x)
			ans
		}, 
		epanechnikov=,parabolic=function(x){
			ans=numeric(length(x))
			ax=abs(x); idx=which(ax<1)
			ans[idx]=-1.5*x[idx]
			ans[ax==1]=NaN
			attributes(ans)=attributes(x)
			ans
		},
		biweight=, quartic=function(x){
			ans=numeric(length(x))
			attributes(ans)=attributes(x)
			ax=abs(x); idx=which(ax<1);x=x[idx]
			ans[idx]=one5d4 * (x* x-1)*x
			ans
		},
		triweight=function(x){
			ans=numeric(length(x))
			attributes(ans)=attributes(x)
			ax=abs(x); idx=which(ax<1); 
			x=x[idx]; x2=x^2
			ans[idx]=none05d16 * x*(x2 -1)^2
			ans
		},
		tricube=function(x){
			ax=abs(x); idx=which(ax<1); 
			ax=ax[idx]; ax2=ax^2; ax3=ax*ax2
			ans=numeric(length(x))
			ans[idx]=0.5*sign(x[idx])*none40d9* ax2*(ax3-1)^2
			attributes(ans)=attributes(x)
			ans
		},
		logistic=function(x){
			# -2*sinh(x/2)^4/sinh(x)^3 # this over-flows earlier
			sx2=sinh(x/2)
			ans=-2*sx2*(sx2/sinh(x))^3
			idx=abs(x)<.Machine$double.xmin*8
			ans[idx]=-x[idx]/8 # linear approximation for x extremely close to 0
			ans
		},
		sech=function(x){
			-onedpi*tanh(x)/cosh(x)
		}
	) # of switch
	attr(ans, 'srcref')=NULL
	if(identical(environment(ans),sys.frame(sys.nframe()))) environment(ans)=pkgEnv
	ans
}, # of function
constEnv) # of substitute
) # of eval 
attr(ddkernel,'srcref')=NULL

d2dkernel=eval(substitute(function(kernel=.kernels)
{
	kernel=match.arg(kernel)
	ans=switch(kernel,
		gaussian=,normal=function(x){
			x2=x*x; #iroot.2pi=1/sqrt(2*base::pi)
			exp(-.5*x2) * (x2-1)*iroot.2pi
		},
		cosine=function(x){
			ans=numeric(length(x))
			ax=abs(x)
			idx=which(ax<1)
			ans[idx]=npi3d16*cos(pid2*x[idx])
			ans[ax==1]=NaN
			attributes(ans)=attributes(x)
			ans
		},
		uniform=, rectangular=function(x){
			ans=numeric(length(x))
			ans[abs(x)==1]=NaN
			attributes(ans)=attributes(x)
			ans
		}, 
		triangular=function(x){
			ans=numeric(length(x))
			ans[abs(x)==1]=NaN
			ans[x==0]=NaN
			attributes(ans)=attributes(x)
			ans
		}, 
		epanechnikov=,parabolic=function(x){
			ans=numeric(length(x))
			ax=abs(x)
			ans[ax<1]=-1.5
			ans[ax==1]=NaN
			attributes(ans)=attributes(x)
			ans
		},
		biweight=, quartic=function(x){
			ans=numeric(length(x))
			attributes(ans)=attributes(x)
			ax=abs(x); x=x[ax<1]
			ans[ax==1]=NaN
			ans[ax<1]=one5d4 * (3*x* x-1)
			ans
		},
		triweight=function(x){
			ans=numeric(length(x))
			ax=abs(x); x2=x[ax<1]^2
			ans[ax<1]=none05d16 * (1 - 6 *x2 + 5 *x2*x2)
			attributes(ans)=attributes(x)
			ans
		},
		tricube=function(x){
			ax=abs(x); idx=which(ax<1); 
			ax=ax[idx]; ax3=ax^3
			ans=numeric(length(x))
			ans[idx]=none40d9* (ax - 5 *ax3*ax + 4 *ax3*ax3*ax)
			attributes(ans)=attributes(x)
			ans
		},
		logistic=function(x){
			# 1/8 *(cosh(x)-2)/cosh(x/2)^4 # this over-flows
			
			ilogit.x=ilogit(as.numeric(x))
			ilogit.x1_x=ilogit.x*(1-ilogit.x)
			
			ilogit.x1_x*((2*ilogit.x-1)^2-2*ilogit.x1_x)
		},
		sech=function(x){
			ax=abs(x) # even function; the code below works better for positive x's 
			ilogit.2x=ilogit(2*ax)
			ilogit.2x1_2x=ilogit.2x*(1-ilogit.2x)
			
			twodpi*ilogit.2x*((2*ilogit.2x-1)^2-4*ilogit.2x1_2x)/exp(ax)
		}
	) # of switch
	attr(ans, 'srcref')=NULL
	if(identical(environment(ans),sys.frame(sys.nframe()))) environment(ans)=pkgEnv
	ans
}, # of function
constEnv) # of substitute
) # of eval 
attr(d2dkernel,'srcref')=NULL

ckernel=fourier.kernel=eval(substitute(function(kernel= .kernels, root.2pi=FALSE, const=NULL)
{
	kernel=match.arg(kernel)
	iroot.2pi.this =if(root.2pi) iroot.2pi else 1
	
	ans=eval(substitute(switch(kernel,
	gaussian=,normal=function(s) exp(-.5*s*s) *iroot.2pi.this,
	dgaussian=,dnormal=function(s) exp(-.5*s*s)*s*(0-1i)*iroot.2pi.this,
	cosine=function(s)	pipi*cos(s)/(pipi-4*s*s) * iroot.2pi.this, 
	uniform = , rectangular = function(s){
		sinc(s) * iroot.2pi.this
	},
	triangular = function(s){
		s2=s*s; 
		ans=2*(1-cos(s))/s2 * iroot.2pi.this
		idx=which(abs(s)<1e-2) ## close to 0/0 region: taylor series
		ans[idx]=(1 + s2[idx]*(s2[idx]*one360 - one12)) * iroot.2pi.this
		ans	
	},
	epanechnikov =, parabolic= function(s){
		s2=s*s; 
		ans=3*(sin(s)-s*cos(s))/(s2*s) * iroot.2pi.this
		idx=which(abs(s)<1e-2) ## close to 0/0 region: taylor series
		ans[idx]=(1 + s2[idx]*(s2[idx]*one280-0.1)) * iroot.2pi.this
		ans
	},
	biweight = , quartic=function(s){
		ss=sin(s); s2=s*s; s4=s2*s2
		ans=((45*(ss - s*cos(s)) - 15*s2*ss))/(s4*s) * iroot.2pi.this
		idx=which(abs(s)<5e-2) ## close to 0/0 region: taylor series
		ans[idx]=(1 + s2[idx]*(s2[idx]*(one504 - s2[idx]*one33264)--one14)) * iroot.2pi.this
		ans
	},
	triweight = function(s){
		cs=cos(s); ss=sin(s);
		s2=s*s; s4=s2*s2; s6=s4*s2; 
		ans=(105*(s2*s*cs + 15*(ss -s*cs) - 6*s2*ss))/(s6*s) *iroot.2pi.this
		idx=which(abs(s)<1e-1) ## close to 0/0 region
		ans[idx]=(1 + s2[idx]*( s2[idx]*(one792 - s2[idx]*one61776) -one18))* iroot.2pi.this
		ans
	},
	tricube = function(s){
		cs=cos(s); ss=sin(s);
		s2=s*s; s4=s2*s2; s6=s4*s2
		ans=(two80d9*(20160 - s6 - 9*(2240 - 1120*s2 + 80*s4 - s6)*cs - 
			s*(20160 - 3240*s2 + 108*s4)*ss))/(s4*s6) *iroot.2pi.this
		idx=which(abs(s)<3.5e-1) ## close to 0/0 region
		ans[idx]=(1 + s2[idx]*(s2[idx]*(one528 - s2[idx]*(s2[idx]*one4199040-one37440)) -three5d486)) *iroot.2pi.this
		ans
   },
	logistic =function(s){
		ps=base::pi*s
		ans=ps/sinh(ps) *iroot.2pi.this
		idx=which(abs(s)<1e-2) ## close to 0/0 region
		ps2=ps*ps; ps4=ps2*ps2
		ans[idx]=(1 + ps2[idx]*(ps2[idx]*(seven360 - ps2[idx]*three1d15120) -one6)) *iroot.2pi.this
		ans
	},
	sech = function(s){
		1/cosh(pid2 *s) *iroot.2pi.this
	}
	)  # of switch 
	))
	ans.const=numeric(0L)
	if(!is.null(const)){
		tmp.const=krRkernel(kernel)
		
		if(isTRUE(const)) const=c('exponent','coefficient')
		if(is.character(const[1L])) const=pmatch(const, c('exponent','coefficient'), duplicates.ok=TRUE) else stop('"const" is misspecified.')
		if(any(1L==const)) ans.const=c(ans.const, exponent=tmp.const['r'])
		if(any(2L==const)) ans.const=c(ans.const, coefficient=tmp.const['kr'])
	}
	attr(ans, 'const')=ans.const 
	attr(ans ,'srcref')=NULL
	if(identical(environment(ans),sys.frame(sys.nframe()))) environment(ans)=pkgEnv
	ans
}
,constEnv) # of substitute
) # of eval 
attr(ckernel,'srcref')=attr(fourier.kernel, 'srcref')=NULL

cdkernel=fourier.dkernel=function(kernel= .kernels, root.2pi=FALSE, const=NULL)
{
	ans=fourier.kernel(kernel, root.2pi, const)
	function(s) (0-1i)*s*ans(s)
}

mkernel=eval(substitute(function(kernel=.kernels, order=1:4, R=TRUE)
{
	kernel=match.arg(kernel)
	if(any(order>4L)) order=order[order<=4]
	ans=switch(kernel,
		gaussian=,normal=c(0,1,0,3),
		cosine=c(0, onem8dpi2, 0, cos.4m),		
		uniform=,rectangular=c(0,oneThird,0,.2),
		triangular=c(0,oneSixth,0,oneFifteenth),
		epanechnikov=,parabolic=c(0,.2,0,thrd35),
		biweight=, quartic=c(0,oneSeventh,0,oned21),
		triweight=c(0,oneNinth,0,oned33),
		tricube=c(0,three5d243,0,oned22),
		logistic=c(0, pipid3, 0, pi4.7d15),
		sech=c(0, pipid4, 0, pi4.5d16)
	)
	ans=ans[order]; names(ans)=as.character(order)
	if(isTRUE(R) || length(order)==0L){
		ans=c(ans, krRkernel(kernel)['R'])
	}
	ans
} # of function
,constEnv) # of substitute
) # of eval
attr(mkernel,'srcref')=NULL


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
	function(q)
	{	nx=length(q)
		xdiff=q-rep(knots, each=nx)
		dim(xdiff) = c(nx, n)
		ans=.rowMeans(Kernel(xdiff/bw), nx, n)
		attributes(ans)=attributes(q)
		ans
	}
}


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
	tk=seq.int(from=from, to=to, length.out=M+1L)
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
		ans[[i]] = approxfun(tk0, pmax.int(0,Re(zetak)), ties='ordered',yleft=0, yright=0)
	}
	structure(if(length(ans)==1L) ans[[1L]] else ans, 
		environment=environment()
	)
}



ddkde=function(x, bw=bw.nrd, kernel=.kernels, from=min(x)-3*median(bw), to=max(x)+3*median(bw))
{# implementation based on FFT in Silverman 1984
    if(is.character(bw)) bw=get(bw, mode='function')
    if(is.function(bw))  bw=bw(x)
	n=length(x)
	from0=from; force(to)
	x=x-from; to=to-from; from=0
	M=nextn(n*2L, 2L)
	delta=(to-from)/M
	tk=seq.int(from=from, to=to, length.out=M+1L)
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
	ftkern=fourier.dkernel(kernel, root.2pi=FALSE)
	
	tk0=tk[1:M]+from0; 
	ans=vector('list', length(bw))
	for(i in seq_along(bw)){
		zetalstar=ftkern(bw[i]*sl) * Yl
		#zetak.bak = colSums(exp((l%o%kk)/M*(1i)*-2*pi)*zetalstar) ## correct but slow
		#zetak=fft(zetalstar)*exp(1i*pi*kk)
		zetak=fft(zetalstar)*c(1,-1)
		ans[[i]] = approxfun(tk0, Re(zetak)/bw[i], ties='ordered',yleft=0, yright=0)
	}
	structure(if(length(ans)==1L) ans[[1L]] else ans, 
		environment=environment()
	)
}

