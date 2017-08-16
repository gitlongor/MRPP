mean.mrpp=function(x, ...)
{
	cumulant(x, order=1L)
}

var=function(x,...)UseMethod('var')
var.default=local({
	ans=stats::var
	formals(ans)=c(formals(stats::var), alist(...=))
	ans
})
var.mrpp=function(x,...)
{
	cumulant(x, order=2L)
}

sd=function(x,...)UseMethod('sd')
sd.default=local({
	ans=stats::sd
	formals(ans)=c(formals(stats::sd), alist(...=))
	ans
})
if(FALSE)sd.mrpp=function(x,...)
{
	call=match.call()
	call[[1L]]=as.name('var')
	envir=parent.frame()
	sqrt(eval(call, envir=envir))
}
sd.mrpp=function(x,...)
{	.Generic='var'
	sqrt(NextMethod())
}

skewness=function(x,...)UseMethod('skewness')
skewness.default=local({
	ans=moments::skewness
	formals(ans)=c(formals(moments::skewness), alist(...=))
	ans
})
skewness.mrpp=function(x,...)
{
	cum=cumulant(x, order=2:3)
	cum[2L]/cum[1L]^1.5
}

kurtosis=function(x,excess=TRUE,...)UseMethod('kurtosis')
kurtosis.default=function(x, excess=TRUE, na.rm=FALSE, ...)
	moments::kurtosis(x,na.rm,...) - if(excess) 3 else 0
kurtosis.mrpp=function(x, excess=TRUE, ...)
{
	cum=cumulant(x, order=c(2L,4L))
	cum[2L]/cum[1L]/cum[1L] + if(!excess) 3 else 0
}


moment=function(x, order, central,...)UseMethod('moment')
moment.default=function(x, order=1:4, central=FALSE, ...)
{
	order=as.integer(order)
	stopifnot(all(order>=0))
	mOrd=max(order)
	if(mOrd<2L){
		ans=rep(1,length(order))
		ans[order==1L]=if(central) 0 else mean(x)
		return(ans)
	}
	structure(moments::all.moments(x, order.max=mOrd, central=central,...)[order+1L], names=as.character(order))
}
moment.mrpp=function(x, order=1:4, central=FALSE,...)
{
	order=as.integer(order)
	mOrd=max(order)
	if(mOrd>=5L) .NotYetImplemented()
	cum=cumulant(x, seq_len(mOrd))
	
	ans=numeric(length(order)); names(ans)=as.character(order)

	ans[order==0L]=1
	ans[order==1L]=if(central) 0 else cum[1L]
	ans[order==2L]=if(central) cum[2L] else cum[2L]+cum[1L]^2
	ans[order==3L]=if(central) cum[3L] else cum[3L] + 3*cum[1L]*cum[2L] + cum[1L]^3
	ans[order==4L]=if(central) cum[4L]+3*cum[2L]*cum[2L] else cum[4L] + 4*cum[1L]*cum[3L] + 3*cum[2L]*cum[2L] +6*cum[2L]*cum[1L]*cum[1L]  + cum[1L]^4
	ans
}

cumulant = function(x, order, ...)UseMethod('cumulant')
cumulant.default=function(x, order=1:4, ...)
{
	order=as.integer(order)
	stopifnot(all(order>=0))
	mOrd=max(order)
	mu.raw = switch(as.character(mOrd), 
		'0'=1, '1'=c(1, mean(x)),
		moments::all.moments(x, order.max=mOrd, central=FALSE, ...))
	ans = c(0, moment2cumulant(mu.raw[-1L]))
	structure(ans[order+1L], names=as.character(order))
}
cumulant.mrpp=eval(substitute(	
	function(x, order=1:4,...)
	{
		order=as.integer(order)
		if(any(order)>=5) .NotYetImplemented()

		ans=seq_along(order)
		names(ans)=as.character(order)
		ans[order==0L]=0
		mOrd = max(order)
		if(mOrd==0L) return(ans)
		
		distObj = x$distObj
		#permutedTrt = x$permutedTrt
		weight.trt  = x$weight.trt
		n			= x$n
		N			= x$nobs
		ntrt		= length(n)
		
		if(max(n)<=1L){
			return(rep(0,length(order)))
		}
		
		Nc=function(C)factorial.rising(max(0,N-C+1L), C)
		nc=function(n, C)factorial.rising(max(0,n-C+1L), C)
		`%//%`=function(numer, denom)ifelse(denom==0L, 0, as.numeric(numer/denom))
		Ncs=integer(8L)
		if( N > rising.fact.int.bound[switch(mOrd, 2L, 4L, 6L, 8L)] ) 
			Ncs=as.bigq(Ncs)
		if(FALSE && max(n) > rising.fact.int.bound[switch(mOrd, 2L, 4L, 6L, 8L)] ){
			`The "cumulant" function for large sample sizes`=function().NotYetImplemented()
			`The "cumulant" function for large sample sizes`()
		}
		
		dmat = as.matrix(distObj)
		
		if(mOrd>=1L) {
			Ncs[[2L]]=Nc(2L)
			d1=2*sum(distObj)
			D1 =d1/Ncs[[2L]]
			ans[order==1L]=mu = as.numeric(D1)
		}
		
		if(mOrd>=2L){
			Ncs[[3L]]=Nc(3L)
			Ncs[[4L]]=Nc(4L)
			w2=weight.trt^2
			n2=sapply(n, nc, C=2)

			d2 = 2*sum(distObj^2)
			D2 = d2%//%Ncs[[2L]]
			
			d1I=rowSums(dmat)
			sum.d1I2=sum(d1I^2)
			D2p=(sum.d1I2-d2)%//%Ncs[[3L]]
			D2pp=(d1^2-4*Ncs[[3L]]*D2p-2*d2)%//%Ncs[[4L]]
			
			ans[order==2L] = sig2 = as.numeric(2 * ( sum(w2%//%n2)  - 1L/Ncs[[2L]]) *
			(D2 - 2 * D2p + D2pp) +
			4 * ( sum(w2 %//% n) - 1/N) * (D2p - D2pp)  )
		}
		
		if(mOrd>=3L){
			Ncs[[5L]]=Nc(5L)
			Ncs[[6L]]=Nc(6L)
			w3=w2*weight.trt
			n2.2=n2*n2
			n2.3=n2.2*n2
			n3=sapply(n, nc, C=3)
			n4=sapply(n, nc, C=4)
			n5=sapply(n, nc, C=5)
			n6=sapply(n, nc, C=6)
			
			d3=2*sum(distObj^3)
			D3=d3%//%Ncs[[2L]]
			dmat2=dmat^2
			ddmat=dmat%*%dmat
			d2I=rowSums(dmat2)
			sum.d1I.d2I=sum(d1I*d2I)
			D3p=(sum.d1I.d2I-d3)%//%Ncs[[3L]]
			D3pp=(d1*d2-4*Ncs[[3L]]*D3p-2*d3)%//%Ncs[[4L]]
			#D3s=6*sum(combn(N,3L,function(idx)distObj[[ idx[1L], idx[2L] ]] * distObj[[ idx[1L], idx[3L] ]] * distObj[[ idx[2L], idx[3L] ]])) %//% Ncs[[3L]]
			D3s = sum(dmat*ddmat)%//%Ncs[[3L]] # numerator = term S_{4}^{(3)} in Siemiatycki 1978
			D3ss=(sum(dmat*(d1I%o%d1I))-Ncs[[3L]]*(2*D3p+D3s)-d3)%//%Ncs[[4L]]
			D3s3=(sum(d1I^3)-3*Ncs[[3L]]*D3p-d3)%//%Ncs[[4L]]
			D3p3=(Ncs[[3L]]*(d1*D2p-4*D3p-2*D3s)-2*Ncs[[4L]]*(2*D3ss+D3s3))%//%Ncs[[5L]]
			D3p4=(Ncs[[4L]]*(d1*D2pp-4*D3pp-8*D3ss)-8*Ncs[[5L]]*D3p3)%//%Ncs[[6L]]

			mean.delta3 = 
				4*D3*sum(w3%//%n2.2)+
				8*(3*D3p+D3s)*sum(w3*n3%//%n2.3)+
				8*(3*D3ss+D3s3)*sum(w3*n4%//%n2.3)+
				6*D3pp*sum(w2*(1-weight.trt+weight.trt*n4%//%n2.2)%//%n2)+
				12*D3p3*sum(w2*((1-weight.trt)*n3+weight.trt*n5%//%n2)%//%n2.2)+
				D3p4*sum(weight.trt*((1-weight.trt)*(1-2*weight.trt)+3*weight.trt*(1-weight.trt)*n4%//%n2.2+w2*n6%//%n2.3))
			ans[order==3L] = as.numeric(mean.delta3 - 3*mu*sig2-mu^3)
		}

		if(mOrd>=4L){
			Ncs[[7L]]=Nc(7L); 
			Ncs[[8L]]=Nc(8L)
			w4=w2*w2
			wow=tcrossprod(weight.trt)
			n2.4=n2.2*n2.2
			n2.o.n2=tcrossprod(n2, n2)
			n3.o.n2=tcrossprod(n3, n2)
			n4.o.n2=tcrossprod(n4, n2)
			n5.o.n2=tcrossprod(n5, n2)
			n6.o.n2=tcrossprod(n6, n2)
			n3.o.n3=tcrossprod(n3)
			n4.o.n3=tcrossprod(n4, n3)
			n4.o.n4=tcrossprod(n4, n4)
			n7=sapply(n, nc, C=7)
			n8=sapply(n, nc, C=8)
		
			d1It=tcrossprod(d1I, d1I)
			zeros = numeric(12L) 
			if(!is.integer(Ncs)) zeros =as.bigq(zeros)
			
			Sa4=c(
				sum(distObj^4)*2, 
				sum(rowSums(dmat2*dmat)*d1I), 
				sum(d2I^2),
				sum(dmat2 * ddmat),
				sum(d2I*d1I^2), 
				sum(dmat2*d1It),
				sum(dmat *tcrossprod(d2I, d1I)), 
				sum(ddmat^2),
				sum(d1I * colSums(dmat*ddmat)),
				sum((dmat%*%d1I)^2), 
				sum(d1I^4),
				sum(dmat*tcrossprod(d1I^2, d1I)), 
				d2^2,
				sum.d1I2^2,
				d2*sum.d1I2, 
				d2*d1*d1, 
				sum.d1I2*d1*d1, 
				d1*d3,
				d1*sum.d1I.d2I,
				d1* sum(dmat * d1It), 
				d1*sum(dmat * ddmat), 
				d1*sum(d1I^3)
			)
			Pa4 = drop(.order4.S2P.mat %*% Sa4)
			Pa4 = c(Pa4,  d1^4 - .order4.f.alpha %*% Pa4 )
			tmp= # avoiding division by 0
			Ds=c(zeros, Pa4[.order4.P2D.ord]%//%Ncs[rep(2:8,c(1L,3L,7L,6L,4L,1L,1L))])
			
			Ez3100=tcrossprod(w3%//%n2.3, weight.trt%//%n2) * (
				4*n2.o.n2*Ds[17L] +
				8*n3.o.n2*(3*Ds[25L] + Ds[29L]) +
				2*n4.o.n2*(3*Ds[30L] + 4*Ds[31L] + 12*Ds[33L]) +
			   12*n5.o.n2*Ds[34L]+
				  n6.o.n2*Ds[35L]
			)
			Ez3100 =structure(as.double(Ez3100), dim=dim(Ez3100))
			Ez3100 = sum(Ez3100) - sum(diag(Ez3100))

			Ez4000 = sum(
				w4%//%n2.4*(
					 8*n2*Ds[13L] +
					16*n3*(4*Ds[14L] + 3*Ds[15L] + 6*Ds[16L]) +
					 4*n4*(4*Ds[17L] + 3*Ds[18L]) +
					48*n4*(2*Ds[19L] + 2*Ds[20L] + 4*Ds[21L] + Ds[22L] + 4*Ds[23L]) +
					16*n5*(3*Ds[24L] + 6*Ds[25L] + Ds[26L]) +
					32*n5*(6*Ds[27L] + 6*Ds[28L] + Ds[29L]) +
					 4*n6*(3*Ds[30L] + 8*Ds[31L]) +
					48*n6*(Ds[32L] + 2*Ds[33L]) + 
					24*n7*Ds[34L] +
					   n8*Ds[35L]
				)
			)

			Ez2200=tcrossprod(w2%//%n2.2, w2%//%n2.2) * (
					4*n2.o.n2* Ds[18L] +
					8*(n3.o.n2+t(n3.o.n2))*Ds[24L] +
					2*(n4.o.n2+t(n4.o.n2))*Ds[30L] +
					16*n3.o.n3*Ds[32L] +
					4*(n4.o.n3+t(n4.o.n3))*Ds[34L] +
					   n4.o.n4*Ds[35L]
			)
			Ez2200 =structure(as.double(Ez2200), dim=dim(Ez2200))
			Ez2200=sum(Ez2200)-sum(diag(Ez2200))
			
			Ez2110= 0 
			if(ntrt>=3L) for(ii in seq_len(ntrt)) {
				wow0=wow; wow0[ii,]=wow0[,ii]=0
				tmp = w2[ii]%//%n2.2[ii] * wow0 * (
					2*n2[ii] * Ds[30L] +
					4*n3[ii] * Ds[34L] +
					  n4[ii] * Ds[35L] 
				)
				#tmp[ii,]=tmp[,ii]=0  	# tmp[,ii]=0 part does not work for bigq matrices; so set corresponding wow elements to zero as above.
				tmp=structure(as.double(tmp), dim=dim(tmp))
				Ez2110=Ez2110+sum(tmp)-sum(diag(tmp))
			}
			
			Ez1111=if(ntrt<4L) 0 else 
				sum(combn(ntrt, 4L, function(ii)prod(weight.trt[ii]))) *24*Ds[35L]
			
			mean.delta4=Ez4000+4*Ez3100+3*Ez2200+6*Ez2110+Ez1111
			mu2=mu*mu
			mean.delta2=sig2+mu2
			ans[order==4L] = as.numeric(mean.delta4 - 4*mean.delta3*mu-3*mean.delta2^2 + 12*mean.delta2*mu2-6*mu2*mu2)
		}
		ans
	} ## of function
  ,constEnv) ## of substitute
) ## of eval

if(FALSE){
naiveEnumSum=function(ord,n, idx.List, delta, verbose=FALSE) ## enumeration of the summation in D functions
{
	seqn=seq_len(n)

	idx.List=lapply(idx.List, sprintf,fmt='J%d')
	idx.List=sapply(idx.List, paste0, collapse=',')
	
	txt=c(
		'seqn=seq_len(n);\n',
		'ans=0;\n'
	)
	i=length(txt)+1L
	for(r in seq_len(ord)){
		txt[[i]]=sprintf('for(J%d in seqn){\n', r)
		if(r==ord){
			i=i+1L
			txt[[i]]=sprintf('if(length(unique(c(%s)))!=%d) next\n', 
				paste0(paste0('J',seq_len(ord)), collapse=','), 
				ord
			)
			tmpc='ans=ans+1'
			for(tmp in seq_along(idx.List))
				tmpc=paste(tmpc,'*delta[',idx.List[tmp], ']', collapse='')
			i=i+1L
			txt[[i]]=paste0(tmpc, '\n', collapse='')
		}else i=i+1L
	}
	for(r in seq_len(ord)){
		i=i+1L
		txt[[i]]='}\n'
	}
	txt[length(txt)+1L]='ans'
	if(verbose)cat(txt)
	eval(parse(text=txt))
}

nc=function(n, C)factorial.rising(max(0,n-C+1L), C)
d13=function(distMat) naiveEnumSum(ord=2, n=NROW(distMat), idx.List=c(rep(list(1:2),4)    ), distMat,FALSE)/nc(NROW(distMat), 2)
d14=function(distMat) naiveEnumSum(ord=3, n=NROW(distMat), idx.List=c(rep(list(1:2),3),list(c(1,3))    ), distMat,FALSE)/nc(NROW(distMat), 3)
d15=function(distMat) naiveEnumSum(ord=3, n=NROW(distMat), idx.List=c(rep(list(1:2),2),list(c(1,3), c(1,3))    ), distMat,FALSE)/nc(NROW(distMat), 3)
d16=function(distMat) naiveEnumSum(ord=3, n=NROW(distMat), idx.List=c(rep(list(1:2),2),list(c(1,3), c(2,3))    ), distMat,FALSE)/nc(NROW(distMat), 3)
d17=function(distMat) naiveEnumSum(ord=4, n=NROW(distMat), idx.List=c(rep(list(1:2),3),list(c(3,4))    ), distMat,FALSE)/nc(NROW(distMat), 4)
d18=function(distMat) naiveEnumSum(ord=4, n=NROW(distMat), idx.List=c(rep(list(1:2),2),rep(list(c(3,4)),2)    ), distMat,FALSE)/nc(NROW(distMat), 4)
d19=function(distMat) naiveEnumSum(ord=4, n=NROW(distMat), idx.List=c(rep(list(1:2),2),list(c(1,3),c(1,4))    ), distMat,FALSE)/nc(NROW(distMat), 4)
d20=function(distMat) naiveEnumSum(ord=4, n=NROW(distMat), idx.List=c(rep(list(1:2),2),list(c(1,3),c(2,4))    ), distMat,FALSE)/nc(NROW(distMat), 4)
d21=function(distMat) naiveEnumSum(ord=4, n=NROW(distMat), idx.List=c(rep(list(1:2),2),list(c(1,3),c(3,4))    ), distMat,FALSE)/nc(NROW(distMat), 4)
d22=function(distMat) naiveEnumSum(ord=4, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(2,4),c(3,4)), distMat,FALSE)/nc(NROW(distMat), 4)
d23=function(distMat) naiveEnumSum(ord=4, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(2,3),c(1,4)), distMat,FALSE)/nc(NROW(distMat), 4)
d24=function(distMat) naiveEnumSum(ord=5, n=NROW(distMat), idx.List=list(c(1,2),c(1,2),c(3,4),c(3,5)), distMat,FALSE)/nc(NROW(distMat), 5)
d25=function(distMat) naiveEnumSum(ord=5, n=NROW(distMat), idx.List=list(c(1,2),c(1,2),c(1,3),c(4,5)), distMat,FALSE)/nc(NROW(distMat), 5)
d26=function(distMat) naiveEnumSum(ord=5, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(1,4),c(1,5)), distMat,FALSE)/nc(NROW(distMat), 5)
d27=function(distMat) naiveEnumSum(ord=5, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(3,4),c(4,5)), distMat,FALSE)/nc(NROW(distMat), 5)
d28=function(distMat) naiveEnumSum(ord=5, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(1,4),c(2,5)), distMat,FALSE)/nc(NROW(distMat), 5)
d29=function(distMat) naiveEnumSum(ord=5, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(2,3),c(4,5)), distMat,FALSE)/nc(NROW(distMat), 5)
d30=function(distMat) naiveEnumSum(ord=6, n=NROW(distMat), idx.List=list(c(1,2),c(1,2),c(3,4),c(5,6)), distMat,FALSE)/nc(NROW(distMat), 6)
d31=function(distMat) naiveEnumSum(ord=6, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(1,4),c(5,6)), distMat,FALSE)/nc(NROW(distMat), 6)
d32=function(distMat) naiveEnumSum(ord=6, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(4,5),c(4,6)), distMat,FALSE)/nc(NROW(distMat), 6)
d33=function(distMat) naiveEnumSum(ord=6, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(2,4),c(5,6)), distMat,FALSE)/nc(NROW(distMat), 6)
d34=function(distMat) naiveEnumSum(ord=7, n=NROW(distMat), idx.List=list(c(1,2),c(1,3),c(4,5),c(6,7)), distMat,FALSE)/nc(NROW(distMat), 7)
d35=function(distMat) naiveEnumSum(ord=8, n=NROW(distMat), idx.List=list(c(1,2),c(3,4),c(5,6),c(7,8)), distMat,FALSE)/nc(NROW(distMat), 8)

}