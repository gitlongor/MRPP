transform <- function(y,
	method=c('nbinom.vst',	'nbinom.quadLL',	'nbinom.symm',
			 'pois.vst',  	'pois.quadLL',  	'pois.symm', 
			 'binom.vst', 	'binom.quadLL',  	'binom.symm', 
			 'chisq.vst',	'chisq.quadLL',		'chisq.symm',
			 'fisher.z',), 
	...)
{
	stopifnot(is.numeric(y))
	att=attributes(y)
	if(!exists(method[1L], mode='function')) method=match.arg(method)
	ans=do.call(method[1L], list(y, ...))
	attributes(ans)=att
	ans
}

constEnv=new.env()
constEnv$oneThird = 1/3
constEnv$oneSixth = 1/6
constEnv$twoThirds = 2/3
constEnv$beta13.13=beta(1/3,1/3)
constEnv$beta23.23=beta(2/3,2/3)

nbinom.vst = function(y, size, tau=1/size)
{
	d=dim(y); 
	if(length(d)==2L && length(tau)==d[2L]) tau=rep(tau,each=d[1L])
	structure(2/sqrt(tau)*asinh(sqrt(y*tau)), dim=d)
}


nbinom.quadLL = local({  # just to shut up CRAN checker
	oneThird = 1/3
	oneSixth = 1/6
	beta13.13=beta(oneThird,oneThird)
	function(y, size, tau=1/size, standardize = FALSE)
	{
		d=dim(y); 
		if(length(d)==2L && length(tau)==d[2L]) tau=rep(tau,each=d[1L])
		ans = pbeta(tau*y/(1+tau*y), oneThird, oneThird) * beta13.13/tau^oneThird
		structure(
			if(standardize) ans *(y*(1+tau*y))^oneSixth else ans
			,dim=d
		)
	}
})
environment(nbinom.quadLL) = constEnv

nbinom.symm = local({  # just to shut up CRAN checker
	oneThird = 1/3
	twoThirds = 2/3
	beta23.23=beta(2/3,2/3)
	function(y, size, tau=1/size)
	{
		d=dim(y); 
		if(length(d)==2L && length(tau)==d[2L]) tau=rep(tau,each=d[1L])
		ans = 3*y^twoThirds/(1+tau*y)^oneThird - beta23.23/tau^twoThirds * pbeta(tau*y/(1+tau*y), twoThirds, twoThirds)
		structure(ans, dim=d)
	}
})
environment(nbinom.symm) = constEnv

pois.vst = sqrt
pois.quadLL = local({oneThird = 1/3; function(y) y^oneThird})
pois.symm = local({twoThirds = 2/3; function(y) y^twoThirds})
environment(pois.quadLL)=environment(pois.symm) = constEnv

binom.vst = function(y) asin(sqrt(y))
chisq.symm = local{oneThird = 1/3; function(y) y^oneThird})

binom.quadLL = 
binom.symm = 
function(y, ...).NotYetImplemented()


fisher.z = atanh


if(FALSE){  ## penglh's original version
	Transform <- function(y,a,method="Log"){
	  if (is.vector(y)){
		if(method=="Log"){
		  result <- log(y+1)
		}else if(method=="VS"){
		  result <- 2/sqrt(a)*asinh(sqrt(y*a))
		}else if(method=="Qua1"){
		  result <- pbeta(a*y/(1+a*y),1/3,1/3)
		}else{
		  result <- pbeta( a*y/(1+a*y),1/3,1/3)*(y/a*(1+a*y)/a)^(1/6)*beta(1/3,1/3)
		}
		return(result)
	  }else{
		n <- dim(y)[2]
		p <- dim(y)[1]
		if(method=="Log"){
		  result <- apply(y,2,function(x) log(x+1))
		}else if(method=="VS"){
		  result <- apply(y,2,function(x) 2/sqrt(a)*asinh(sqrt(x*a)))
		}else if(method=="Qua1"){
		  result <- apply(y,2,function(x) pbeta(a*x/(1+a*x),1/3,1/3))
		}else{
		  result <- apply(y,2,function(x) pbeta( a*x/(1+a*x),1/3,1/3)*(x/a*(1+a*x)/a)^(1/6)*beta(1/3,1/3))
		}
		return(result)
	  }
	}
}