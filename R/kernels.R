
pkernel=function(kernel= c("gaussian", "biweight", 'triweight', 'tricube',"logistic"))
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

dkernel=function(kernel= c("gaussian", "biweight", 'triweight', 'tricube',"logistic"))
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

pkde=function(x, bw=bw.nrd, kernel=c("gaussian", "biweight", 'triweight', 'tricube',"logistic"))
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

