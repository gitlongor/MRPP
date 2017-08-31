pearson3.env=new.env(hash=FALSE)
pearson3.env$pearson3.common=expression({
	n=max(length(y),length(mean),length(sd),length(skew))
	twoskew=2/skew
	twoskew2=twoskew*twoskew
	skew2=.5*skew
	skew=rep_len(skew, n)
})
pearson3.env$rq.common=expression({
	ifelse(skew>0, 
		(y-twoskew)*sd+mean,
	ifelse(skew<0,
		mean-(y+twoskew)*sd,
		y*sd+mean
	))
})
pearson3.env$twoskew=pearson3.env$twoskew2=pearson3.env$skew2=NULL

dpearson3=function(x, mean, sd, skew, log=FALSE)
{
	y=x
	eval(pearson3.common)
	std=(y-mean)/sd
	ans=ifelse(skew!=0,
		dgamma(sign(skew)*(std+twoskew), twoskew2, scale=abs(skew*.5), log=log), 
		dnorm(std, log=log)
	)
	if(log) ans-log(sd) else ans/sd 
}

d2dpearson3=function(x, mean, sd, skew)
{
	y=x
	eval(pearson3.common)
	std=(x-mean)/sd
	shape=k=twoskew2; scale=theta=abs(skew*.5)
	y=sign(skew)*(std+twoskew)
	ifelse(skew!=0,
		dgamma(y, shape=k, scale=theta)
		*(theta^-2-2*(k-1)/theta/y+(k-1)*(k-2)/y^2)
		/sd^3, 
		dnorm(std)*(std^2-1)/sd^3
	)
}

ppearson3=function(q, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
{
	y=q
	eval(pearson3.common)
	std=(y-mean)/sd
	ifelse(skew>0,
		pgamma(std+twoskew, twoskew2, scale=skew2, lower.tail=lower.tail, log.p=log.p), 
	ifelse(skew<0,
		pgamma(-std-twoskew, twoskew2, scale=-skew2, lower.tail=!lower.tail, log.p=log.p), 
		pnorm(std, lower.tail=lower.tail, log.p=log.p)
	))
}

qpearson3=function(p, mean, sd, skew, lower.tail=TRUE, log.p=FALSE)
{
	y=p
	eval(pearson3.common)
	y=qgamma(p, twoskew2, scale=abs(skew2), lower.tail=lower.tail, log.p=log.p)
	eval(rq.common)
}

rpearson3=function(n, mean, sd, skew)
{
	y=rep_len(n,n); n0=n
	eval(pearson3.common)
	y=rgamma(n0, twoskew2, scale=abs(skew2))
	eval(rq.common)
}

environment(dpearson3)=
environment(ppearson3)=
environment(qpearson3)=
environment(d2dpearson3)=
environment(rpearson3)=pearson3.env


dpearson3gca=function(x, mean, sd, skew, exkurt, log=FALSE)
{
	y=x
	eval(pearson3.common)
	std=(y-mean)/sd
	stdmoms=cumulant2moment(c(abs(twoskew),1,abs(skew[1L]),exkurt[1L]))
	ans=ifelse(skew>0,
		dapx_gca(std+twoskew, raw.moments=stdmoms,support=c(0,Inf),basis='gamma',basepar=list(shape= twoskew2, scale=skew2), log=log), 
	ifelse(skew<0,
		dapx_gca(-std-twoskew, raw.moments=stdmoms, support=c(0,Inf),basis='gamma', basepar=list(shape=twoskew2, scale=-skew2), log=log), 
		dapx_gca(std, raw.moments=stdmoms, support=c(-Inf,Inf),basis='normal', log=log)
	))
	if(log) ans-log(sd) else ans/sd
}


ppearson3gca=function(q, mean, sd, skew, exkurt, lower.tail=TRUE, log.p=FALSE)
{
	y=q
	eval(pearson3.common)
	std=(y-mean)/sd
	stdmoms=cumulant2moment(c(abs(twoskew),1,abs(skew[1L]),exkurt[1L]))
	ifelse(skew>0,
		papx_gca(std+twoskew, raw.moments=stdmoms,support=c(0,Inf),basis='gamma',basepar=list(shape= twoskew2, scale=skew2), lower.tail=lower.tail, log.p=log.p), 
	ifelse(skew<0,
		papx_gca(-std-twoskew, raw.moments=stdmoms, support=c(0,Inf),basis='gamma', basepar=list(shape=twoskew2, scale=-skew2), lower.tail=!lower.tail, log.p=log.p), 
		papx_gca(std, raw.moments=stdmoms, support=c(-Inf,Inf),basis='normal', lower.tail=lower.tail, log.p=log.p)
	))
}

environment(dpearson3gca)=
environment(ppearson3gca)=pearson3.env

dgammagca=function(x, mean, sd, skew, exkurt, log=FALSE)
{
	moms=cumulant2moment(c(mean, sd*sd, skew*sd*sd*sd, exkurt*sd^4))
	if(moms[2L]==0)moms[]=0
	dapx_gca(x, moms, support=c(0,Inf), basis='gamma', log=log)
}
pgammagca=function(q, mean, sd, skew, exkurt, lower.tail = TRUE, 
    log.p = FALSE)
{
	moms=cumulant2moment(c(mean, sd*sd, skew*sd*sd*sd, exkurt*sd^4))
	if(moms[2L]==0)moms[]=0
	papx_gca(q, moms, support=c(0,Inf), basis='gamma',lower.tail = TRUE, 
    log.p = FALSE)
}

if(FALSE){
dgammagca0=function(x, mean, sd, skew, exkurt, log=FALSE)
{# This follows from Berberan-Santos, 2006, Journ Math. Chemistry, 42, 585--593 
 # as well as eq(12.52) on page 25 of Johnson,Kotz,Balakrishnan (1994) Continuous Univariate Distributions, 2ed.
 # But this does not seem to work well. 
	a=(mean/sd)^2; b=sd*sd/mean
	third=
	((-6 + 11*a - 6*a^2 + a^3)*b^3 - 
  3*(2 - 3*a + a^2)*b^2*x + 3*(-1 + a)*b*x^2 - 
  x^3)/(b^3*x^3)
    fourth=
	((24 - 50*a + 35*a^2 - 10*a^3 + a^4)*b^4 - 
  4*(-6 + 11*a - 6*a^2 + a^3)*b^3*x + 
  6*(2 - 3*a + a^2)*b^2*x^2 - 
  4*(-1 + a)*b*x^3 + x^4)/(b^4*x^4)
    ans1=dgamma(x,a,scale=b)
	ans2=pmax(-Inf, 1-(skew*sd^3-2*a*b^3)/6*third + (exkurt*sd^4-6*a*b^4)/24*fourth)
	ans1*ans2
}
dgammagca0=function(x, mean, sd, skew, exkurt, log=FALSE)
{# Bowers (1966) Transactions of Society of Actuaries, 18, 125-147
 # this agrees with dgammagca
	a=(mean/sd)^2; b=sd*sd/mean
	y=x/b
	cp = cumprod(a+2:0)
	a3=cp[3L]
	third=as.function(polynomial(c(c(-1,3,-3)*rev(cp), 1)))(y)
	cp=c(a+3, (a+3)*cp )
	a4=cp[4L]
	fourth=as.function(polynomial(c(c(1,-4,6,-4)*rev(cp), 1)))(y)
	
	sd2=sd/b*sd/b
	c3=skew*sd2*sd/b
	c4=(exkurt+3)*sd2*sd2
	ans1=dgamma(y,a,scale=1,log=log)
	ans2=pmax(-Inf, 1+
		1/a3/6*(c3-2*a)*third +
		1/a4/24*(c4-12*c3-3*a*a+18*a)*fourth
		)
	if(log) ans1+log(ans2)-log(b) else ans1*ans2/b
}


mom0=function(x,ek=1.5,pdf=dgammagca)pdf(x,3,2,4/3,ek)
mom1=function(x,ek=1.5,pdf=dgammagca)x*pdf(x,3,2,4/3,ek)
mom2=function(x,ek=1.5,pdf=dgammagca)(x-3)^2*pdf(x,3,2,4/3,ek)
mom3=function(x,ek=1.5,pdf=dgammagca)(x-3)^3/2^3*pdf(x,3,2,4/3,ek)
mom4=function(x,ek=1.5,pdf=dgammagca)((x-3)^4/2^4-3)*pdf(x,3,2,4/3,ek)
mom=function(ek=1.5, pdf=dgammagca)c(
	integrate(mom0, -Inf, Inf,ek=ek, pdf=pdf)$value,
	integrate(mom1, -Inf, Inf,ek=ek, pdf=pdf)$value,
	integrate(mom2, -Inf, Inf,ek=ek, pdf=pdf)$value,
	integrate(mom3, -Inf, Inf,ek=ek, pdf=pdf)$value,
	integrate(mom4, -Inf, Inf,ek=ek, pdf=pdf)$value
)
mom(ek=17/3-3,pdf=dgammagca)
mom(ek=17/3-3,pdf=dgammagca0)
mom(ek=2,pdf=dgammagca)
mom(ek=2,pdf=dgammagca0)

integrate(tmp1, -Inf, Inf,ek=2)
integrate(tmp2, -Inf, Inf,ek=2)
integrate(tmp3, -Inf, Inf,ek=2)
integrate(tmp4, -Inf, Inf,ek=2)

dpearson3gca0=function(x, mean, sd, skew, exkurt, log=FALSE)
{## direct implementation according to Berberan-Santos, 2006, Journ Math. Chemistry, 42, 585--593
	y=x
	eval(pearson3.common)
	std=(y-mean)/sd
	stdmoms=cumulant2moment(c(abs(twoskew),1,abs(skew[1L]),exkurt[1L]))
	ask=abs(skew)
	y=sign(skew)*(std+twoskew)
	expr=ifelse(y>0,
	(1/(ask^8* y^4))*8* (32 + 
	   ask* (-64* y + 
		  ask* (-80 + 3* ask^6 + 6* ask^5* y + 48 *y^2 + 4 *ask^3* y *(-11 + y^2) - 
			 16 *ask *y *(-6 + y^2) + ask^4 *(-25 + 6* y^2) + 
			 2 *ask^2 *(35 - 18 *y^2 + y^4))))
	,0)
	ans1=dgamma(y, twoskew2, scale=abs(skew2), log=log)
	ans2=pmax(0, 1 + (exkurt-1.5*ask*ask)/24 * expr)
	if(log) ans1+log(ans2)-log(sd)  else ans1*ans2/sd
}
environment(dpearson3gca0)=MRPP:::pearson3.env

kurtosis(rpearson3(5e6, 3,2,1)) ## should be 1.5 * skew^2

integrate(tmp1, -Inf, Inf,pdf=dpearson3gca0,ek=2)
integrate(tmp2, -Inf, Inf,pdf=dpearson3gca0,ek=2)
integrate(tmp3, -Inf, Inf,pdf=dpearson3gca0,ek=2)
integrate(tmp4, -Inf, Inf,pdf=dpearson3gca0,ek=2)

integrate(tmp01, -Inf, Inf, ek=1.2)
integrate(tmp02, -Inf, Inf, ek=1.2)
integrate(tmp03, -Inf, Inf, ek=1.2)
integrate(tmp04, -Inf, Inf, ek=1.2)

tmp1=function(x)x*dpearson3(x,3,2,1)
tmp2=function(x)(x-3)^2*dpearson3(x,3,2,1)
tmp3=function(x)(x-3)^3/2^3*dpearson3(x,3,2,1)
tmp4=function(x)((x-3)^4/2^4-3)*dpearson3(x,3,2,1)
integrate(tmp1, -Inf, Inf)
integrate(tmp2, -Inf, Inf)
integrate(tmp3, -Inf, Inf)
integrate(tmp4, -Inf, Inf)

tmp1=function(x,ek=1.5)x*dpearson3gca(x,3,2,1,ek)
tmp2=function(x,ek=1.5)(x-3)^2*dpearson3gca(x,3,2,1,ek)
tmp3=function(x,ek=1.5)(x-3)^3/2^3*dpearson3gca(x,3,2,1,ek)
tmp4=function(x,ek=1.5)((x-3)^4/2^4-3)*dpearson3gca(x,3,2,1,ek)
integrate(tmp1, -Inf, Inf)
integrate(tmp2, -Inf, Inf)
integrate(tmp3, -Inf, Inf)
integrate(tmp4, -Inf, Inf)

integrate(tmp1, -Inf, Inf,ek=1.2)
integrate(tmp2, -Inf, Inf,ek=1.2)
integrate(tmp3, -Inf, Inf,ek=1.2)
integrate(tmp4, -Inf, Inf,ek=1.2)
}