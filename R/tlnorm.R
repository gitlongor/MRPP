## Affine transformed Log-Normal Distribtion that Matches 
## the First 3 Specified Moments: 
## Identically distributed as 
## 		sd * sign(skew) * (X - mu_ln) + mean
## where X ~ Log-Normal random variable, with 
## skewness |skew|, variance 1, and mean mu_ln
.tlnorm.common=eval(bquote(quote({
	.x1xm1=function(x)x+1/x-1
	N=max(length(arg),length(mean),length(sd),length(skew))
	arg=rep_len(arg,N); mean=rep_len(mean,N);sd=rep_len(sd,N);skew=rep_len(skew,N)
	skew2=skew*skew
	expsig2 = .x1xm1(( 1+.5*skew2+sqrt(skew2+.25*skew2*skew2) )^.(1/3))
	exp2mu = 1/expsig2/(expsig2-1)
	ln.mu = sqrt(exp2mu * expsig2)
	meanlog = .5*log(exp2mu)
	sdlog = sqrt(log(expsig2))
})))
.tlnorm.dp=quote({
	std=(arg-mean)/sd
	trans=sign(skew)*std + ln.mu
})
.tlnorm.rq=quote({
	ifelse(skew!=0,
		sign(skew)*sd*(arg-ln.mu)+mean, 
		sd*arg+mean
	)
})

dtlnorm=eval(bquote(function(x, mean, sd, skew, log=FALSE)
{
	arg=x
	.(.tlnorm.common)
	.(.tlnorm.dp)
	ans=ifelse(skew!=0,
		dlnorm(trans, meanlog, sdlog, log),
		dnorm(std, log=log)
	)
	if(log) ans-log(sd) else ans/sd
}
))

ptlnorm=eval(bquote(function(q, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
{
	arg=q
	.(.tlnorm.common)
	.(.tlnorm.dp)
	ifelse(skew!=0,
		plnorm(trans, meanlog, sdlog, lower.tail=xor(skew<0, lower.tail), log.p=log.p), 
		pnorm(std, lower.tail=lower.tail, log.p=log.p)
	)
}
))

qtlnorm=eval(bquote(function(p, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
{
	arg=p
	.(.tlnorm.common)
	arg=qlnorm(arg, meanlog, sdlog, lower.tail=xor(skew<0,lower.tail), log.p=log.p)
	.(.tlnorm.rq)
}
))

rtlnorm=eval(bquote(function(n, mean, sd, skew)
{
	arg=rep(n,n)
	.(.tlnorm.common)
	arg=rlnorm(n, meanlog, sdlog)
	.(.tlnorm.rq)
}
))
