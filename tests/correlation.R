ar1mat=function(n=5, rho=.5)
{
	ans=diag(n)
	for(i in seq_len(n-1)+1L){
		seq.row=seq.int(from=1L, to=n-i+1L, by=1L)
		seq.col=seq.int(from=i,to=n,by=1L)
		if(length(seq.row)>1L){
			diag(ans[seq.row, seq.col])=diag(ans[seq.col,seq.row])=rho^(i-1)
		}else ans[seq.row, seq.col]=ans[seq.col,seq.row]=rho^(i-1)
	}
	ans
}

xdat=function(n1=10, n2=n1, mu0=3,rho0=mu0)
{# simulating "X" shaped bi-variate data that can be detected more easily by MRPP
	require(MASS)
	n=n1+n2
	dat=matrix(NA_real_, n, 2L)
	ind=sample(c(TRUE,FALSE), n1, replace=TRUE)
	f=rep(FALSE, n2)
	if(sum(ind)>0)dat[c(ind,f),]=mvrnorm(sum(ind), mu=c(mu0,mu0), Sigma=ar1mat(2,rho0))
	if(sum(ind)<n1)dat[c(!ind,f),]=mvrnorm(n1-sum(ind), mu=-c(mu0,mu0), Sigma=ar1mat(2,rho0))
	ind=sample(c(TRUE,FALSE), n2, replace=TRUE)
	f=rep(FALSE, n1)
	if(sum(ind)>0)dat[c(f,ind),]=mvrnorm(sum(ind), mu=c(-mu0,mu0), Sigma=ar1mat(2,-rho0))
	if(sum(ind)<n2)dat[c(f,!ind),]=mvrnorm(n2-sum(ind), mu=-c(-mu0,mu0), Sigma=ar1mat(2,-rho0))
	attr(dat, 'trt')=as.factor(rep(1:2, c(n1,n2)))
	attr(dat, 'mu0')=mu0
	attr(dat, 'rho0')=rho0
	attr(dat, 'n')=c(n1,n2)
	class(dat)='xdat'
	dat
}
anova.xdat=function(object, ...)
{
	trt=attr(object, 'trt')
	anova(lm(object~trt), ...)
}
mrpp.test.xdat=function(y, ...)
{	trt=attr(y, 'trt')
	y=unclass(y)
	mrpp.test(y~trt, ...)
}
plot.xdat=function(x,y,...)
{
	x=unclass(x)
	plot(x, col=rep(1:2, attr(x,'n')), pch=19,...)
}

cor.xdat=function(object, mu0=attr(object, 'mu0'), rho0=attr(object, 'rho0'))
{
	c(1,-1)*(rho0+mu0^2)/(1+mu0^2)
}

mu0=function(rho0, mu0.rho0.ratio=1)
{
	
	ans=c(mu0=(1 - sqrt(1 + 4 * mu0.rho0.ratio^2*rho0 - 4 *mu0.rho0.ratio^2*rho0^2))/(2 *mu0.rho0.ratio*(-1 + rho0)))
	attr(ans, 'rho0')=ans/mu0.rho0.ratio
	stopifnot(abs(attr(ans, 'rho0'))<=1)
	ans
}


pvals0=replicate(1e3, {tmp=mu0(.9, 5); mrpp.test(xdat(mu0=tmp, rho0=attr(tmp, 'rho0')))$p.value})


pvals1=replicate(1e3,{dat=rbind(mvrnorm(10,c(0,0),ar1mat(2,.9)),mvrnorm(10,c(0,0),ar1mat(2,-.9))); mrpp.test(dat~trt)$p.value})


ks.test(pvals0, pvals1, alternative='greater')


dat=xdat(n1=10, mu0=3, rho0=.5); trt=attr(dat, 'trt')
plot(dat, col=rep(1:2, each=nrow(dat)/2))
anova(lm(dat~trt))
mrpp.test(dat~trt)


dat=xdat(n1=50, mu0=3, rho0=.8); trt=attr(dat, 'trt'); pmtrt=permuteTrt(trt)
plot(dat, col=rep(1:2, each=nrow(dat)/2))
anova(lm(dat~trt))
mrpp.test(dat~trt)
(imps=grad.smoothp(dat, pmtrt, adjust='log scale'))


dat.noise=cbind(dat, matrix(rnorm(100*50),100))
mrpp.test(dat.noise~trt)
(imps.noise=grad.smoothp(dat.noise, pmtrt, adjust='log scale'))

dat.noise2=cbind(dat.noise, matrix(rnorm(100*50),100))
mrpp.test(dat.noise2~trt)
(imps.noise2=grad.smoothp(dat.noise2, pmtrt, adjust='log scale'))

