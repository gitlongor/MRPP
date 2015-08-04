known.transform.methods=c(
	'nbinom.vst',	'nbinom.quadLL',	'nbinom.symm',
	 'pois.vst',  	'pois.quadLL',  	'pois.symm', 
	 'binom.vst', 	'binom.quadLL',  	'binom.symm', 
	 'norm.vst', 	'norm.quadLL',  	'norm.symm', 
	 'chisq.vst',				 		'chisq.symm',
	 'fisher.z'
)
			 
transform.matrix <- function(`_data`,	method='I', ...)
{
	stopifnot(is.numeric(`_data`))
	att=attributes(`_data`)
	if(!exists(method[1L], mode='function')) method=match.arg(method, known.transform.methods)
	toCallArgs=names(formals(args(get(method[1L],mode='function'))))
	ddd=list(...)
	if(! ('...'%in%toCallArgs) ) {
		nmddd=names(ddd)
		ddd = ddd[nmddd=='' | nmddd %in% toCallArgs]
	}
	if(length(ddd) == 0L){
			ans=do.call(method[1L], list(`_data`))
	}else 	ans=do.call(method[1L], c(list(`_data`), ddd))
	attributes(ans)=att
	ans
}

transform.numeric <- function(`_data`, ...)
{
	dim(`_data`)=c(NROW(`_data`), 1L)
	.Class=c(.Class, 'matrix')
	drop(NextMethod(.Generic))
}

nbinom.vst = local({
	twoThree.192=23/192
	twoThree.96 =23/96
	
	function(y, size = 1/invsize, invsize=1/size, adjust=c('none','anscombe1','anscombe2','anscombe3'))
	{
		stopifnot(all(size==1/invsize || invsize==1/size))
		d=dim(y); 
		if(length(d)==2L && length(invsize)==d[2L]) invsize=rep(invsize,each=d[1L])
		adjust = match.arg(adjust)
		if(adjust == 'anscombe2' && size < 2) adjust='anscombe3'
		if(adjust == 'anscombe3' && size < 1) adjust='none'
		if(adjust == 'anscombe1' && size < 1) adjust='none'
		
		ans = switch(adjust, 
			none      = 2/sqrt(invsize)*asinh(sqrt(invsize*y)), 
			anscombe1 = 2/trigamma(size) *asinh(sqrt((y+0.375)/(size-0.75))), 
			anscombe2 = asinh(sqrt((y+0.375+twoThree.192*invsize)/(size-0.75-twoThree.96*invsize))) / trigamma(size), 
			anscombe3 = log(y+0.5*size) / trigamma(size)
		)
		structure(ans, dim=d)
	}
})
environment(nbinom.vst) = constEnv


nbinom.quadLL = local({  # just to shut up CRAN checker
	oneThird = 1/3
	oneSixth = 1/6
	beta13.13=beta(oneThird,oneThird)
	function(y, size = 1/invsize, invsize=1/size,  standardize = FALSE)
	{
		stopifnot(all(size==1/invsize || invsize==1/size))
		d=dim(y); 
		if(length(d)==2L && length(invsize)==d[2L]) invsize=rep(invsize,each=d[1L])
		ans = pbeta(invsize*y/(1+invsize*y), oneThird, oneThird) 
		structure(
			if(standardize) ans *(y*(1+invsize*y))^oneSixth * beta13.13/invsize^oneThird else ans
			,dim=d
		)
	}
})
environment(nbinom.quadLL) = constEnv

nbinom.symm = local({  # just to shut up CRAN checker
	oneThird = 1/3
	twoThirds = 2/3
	beta23.23=beta(2/3,2/3)
	function(y,size = 1/invsize, invsize=1/size)
	{
		stopifnot(all(size==1/invsize || invsize==1/size))
		d=dim(y); 
		if(length(d)==2L && length(invsize)==d[2L]) invsize=rep(invsize,each=d[1L])
		ans = 3*y^twoThirds/(1+invsize*y)^oneThird - beta23.23/invsize^twoThirds * pbeta(invsize*y/(1+invsize*y), twoThirds, twoThirds)
		structure(ans, dim=d)
	}
})
environment(nbinom.symm) = constEnv

pois.vst = function(y, adjust=c('anscombe','none'))
{
	adjust=match.arg(adjust)
	const = if(adjust=='anscombe') 0.375 else 0
	2*base::sqrt(y+const)
}
pois.quadLL = local({
	oneThird = 1/3; 
	oneSixth = 1/6
	function(y, standardize=FALSE)
	{
		ans = y^oneThird
		if(standardize) ans * y^oneSixth * 3 else ans
	}
})
pois.symm = local({
	twoThirds = 2/3; 
	function(y, adjust=c('none','anscombe'))
	{
		adjust=match.arg(adjust)
		const = if(adjust=='anscombe') twoThirds else 0
		(y + const)^twoThirds
	}
})
environment(pois.quadLL)=environment(pois.symm) = constEnv

norm.vst = 
norm.quadLL = 
norm.symm = base::I


binom.vst = function(y, size, adjust=c('anscombe','none'))
{
	adjust=match.arg(adjust)
	if(adjust=='anscombe') {
			2*sqrt(size+.5)*asin(sqrt((y+0.375)/(size+0.75)))
	}else 	2*sqrt(size)*asin(sqrt(y/size))
}
binom.quadLL = local({
	oneThird = 1/3
	oneSixth = 1/6
	oneNinth = 1/9
	beta13.13=beta(oneThird,oneThird)
	function(y, size, adjust=c('none','borges'), standardize=FALSE)
	{
		adjust=match.arg(adjust)
		if(adjust=='borges'){
			y=y+oneSixth; size=size+oneThird
			const = oneSixth
		}else const = 0
		d=dim(y)
		if(length(d)==2L && length(size)==d[2L]) size=rep(size, each=d[1L])
		ans = pbeta(y/size, oneThird, oneThird) 
		structure(
			if(standardize) ans * beta13.13 * (y *(size-y)*size)^oneSixth * sqrt(1+const/size) else ans,
			dim=d
		)
	}
})
environment(binom.quadLL) = constEnv
binom.symm   = local({
	oneSixth = 1/6
	oneThird = 1/3; 
	twoThirds = 2/3; 
	function(y, size, adjust=c('none','borges'), standardize=FALSE)
	{
		adjust=match.arg(adjust)
		d=dim(y)
		if(length(d)==2L && length(size)==d[2L]) size=rep(size, each=d[1L])
		if(adjust=='borges') {
			const.numerator = oneSixth
			const.denominator = oneThird
		}else{
			const.numerator =
			const.denominator = 0
		}
		ans = pbeta((y+const.numerator)/(size+const.denominator), twoThirds, twoThirds) 
		structure(
			if(standardize) ans * beta23.23 / (y/size *(1-y/size))^oneSixth * sqrt(size+const.denominator) else ans,
			dim=d
		)
	}
})
environment(binom.symm) = constEnv


chisq.vst = function(y, df, ncp = 0)
{
	if(ncp!=0) .NotYetImplemented()
	sqrt(2*y) - sqrt(2*df-1)
}
chisq.symm = local({
	oneThird = 1/3; 
	oneSixth= 1/6; 
	sqrt2.3 = sqrt(2)/3
	function(y, df, ncp=0, adjust=c('abdel-aty','none'), standardize=FALSE)
	{
		adjust = match.arg(adjust)
		if(adjust=='abdel-aty') {
			dfncp=df+ncp
			df2ncp=dfncp+ncp
			ans=(y/dfncp)^oneThird
			if(standardize) {
				ans =(ans - 1 +2*df2ncp/9/dfncp^2 ) *
					3*dfncp/sqrt(2*df2ncp)
			}else ans
		}else{
			if(ncp!=0).NotYetImplemented()
			ans = y^oneThird
			if(standardize) {
				sqrt(df)/sqrt2.3 ->tmp
				ans*df^oneSixth/sqrt2.3 - tmp+1/tmp 
			}else ans
		}
	}
})
environment(chisq.symm) = constEnv

#chisq.quadLL = function(y, ...).NotYetImplemented()


fisher.z = local({
	fourThirds =4/3
	function(y, df, adjust=c('none','hotelling'))
	{
		ans = atanh(y)
		adjust=match.arg(adjust)
		if(adjust=='hotelling') (ans-y/2/(df-1.5))*sqrt(df-fourThirds) else ans 
	}
})

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
