known.transform.methods=c(
	'nbinom.vst',	'nbinom.quadLL',	'nbinom.symm',
	 'pois.vst',  	'pois.quadLL',  	'pois.symm', 
	 'binom.vst', 	'binom.quadLL',  	'binom.symm', 
	 'norm.vst', 	'norm.quadLL',  	'norm.symm', 
	 'chisq.vst',	'chisq.quadLL',		'chisq.symm',
	 'fisher.z'
)
			 
transform.matrix <- function(`_data`,	method='I', ...)
{
	stopifnot(is.numeric(`_data`))
	att=attributes(`_data`)
	if(!exists(method[1L], mode='function')) method=match.arg(method, known.transform.methods)
	if(length(formals(args(get(method[1L],mode='function'))))==1L){
			ans=do.call(method[1L], list(`_data`))
	}else 	ans=do.call(method[1L], list(`_data`,...))
	attributes(ans)=att
	ans
}

transform.numeric <- function(`_data`, ...)
{
	dim(`_data`)=c(NROW(`_data`), 1L)
	.Class=c(.Class, 'matrix')
	drop(NextMethod(.Generic))
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

pois.vst = function(y) 2*base::sqrt(y)
pois.quadLL = local({oneThird = 1/3; function(y) y^oneThird})
pois.symm = local({twoThirds = 2/3; function(y) y^twoThirds})
environment(pois.quadLL)=environment(pois.symm) = constEnv

norm.vst = 
norm.quadLL = 
norm.symm = base::I


binom.vst = function(y) asin(sqrt(y))

chisq.symm = local({oneThird = 1/3; function(y) y^oneThird})

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
