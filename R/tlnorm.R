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

ExKurtTlnorm=eval(bquote(function(skew)
{
	arg=mean=sd=skew
	.(.tlnorm.common)
	as.function(polynomial(c(-6,0,3,2,1)))(expsig2)
}
))

## Note: Suppose conditioning on Z~Bernoulli(p), 
## and X|Z=1 ~F_1 and X|Z=0 ~ F_0.  Let k_n denote the nth cumulant. 
## If E(F_1) = E(F_0), then 
##     E(X)=p E(F_1) + (1-p) E(F_0), 
##     V(X)=p V(F_1) + (1-p) V(F_0), 
##     k_3(X) = p k_3(F_1) + (1-p) k_3(F_0), 
##     k_4(X) = p k_4(F_1) + (1-p) k_4(F_0) + 3p(1-p) ( V(F_1) - V(F_2) )^2 .
## If further V(F_1) = V(F_0), then further
##     k_4(X) = p k_4(F_1) + (1-p) k_4(F_0), 
##     skew(X) = p skew(F_1) + (1-p) skew(F_0),
##     kurt(X) = p kurt(F_1) + (1-p) kurt(F_0),
##     exKurt(X) = p exKurt(F_1) + (1-p) exKurt(F_0).
dMixP3Tln=eval(bquote(function(x, mean, sd, skew, exkurt, log=FALSE, proper=FALSE)
{
	ek.tln=ExKurtTlnorm(skew)
	ek.p3=ExKurtPearson3(skew)
	prop.p3=(exkurt-ek.tln)/(ek.p3-ek.tln)
	if(!proper || (0<prop.p3 && prop.p3<1)) {
		ans=dpearson3(x,mean,sd,skew)*prop.p3+dtlnorm(x,mean,sd,skew)*(1-prop.p3)
		if(isTRUE(log))log(ans) else ans 
	}else if(prop.p3<=0){
		dtlnorm(x,mean,sd,skew,log)
	}else dpearson3(x,mean,sd,skew,log)
}
))

pMixP3Tln=eval(bquote(function(q, mean, sd, skew, exkurt, lower.tail=TRUE, log.p=FALSE, proper=TRUE)
{
	ek.tln=ExKurtTlnorm(skew)
	ek.p3=ExKurtPearson3(skew)
	prop.p3=(exkurt-ek.tln)/(ek.p3-ek.tln)
	if(!proper || (0<prop.p3 && prop.p3<1)) {
		ans=ppearson3(q,mean,sd,skew,lower.tail)*prop.p3+ptlnorm(q,mean,sd,skew,lower.tail)*(1-prop.p3)
		if(isTRUE(log.p))log(ans) else ans 
	}else if(prop.p3<=0){
		ptlnorm(q,mean,sd,skew,lower.tail,log.p)
	}else ppearson3(q,mean,sd,skew,lower.tail,log.p)
}
))
