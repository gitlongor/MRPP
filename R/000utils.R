constEnv=new.env(hash=TRUE, size=73L)
pkgEnv=parent.env(constEnv)
constEnv$oneThird = 1/3
constEnv$oneSixth = 1/6
constEnv$twoThirds = 2/3
constEnv$oneNinth = 1/9
constEnv$beta13.13=beta(1/3,1/3)
constEnv$beta23.23=beta(2/3,2/3)
constEnv$twoThree.192=23/192
constEnv$twoThree.96 =23/96
constEnv$sqrt2.3 = sqrt(2)/3
constEnv$fourThirds =4/3
constEnv$one120 = 1/120
constEnv$thr5d32 = 35/32
constEnv$two1d32 = 21/32
constEnv$five32 = 5/32
constEnv$sev0d81 = 70/81
constEnv$one0d27 = 10/27
constEnv$sev81 = 7/81
constEnv$thr5d54 = 35/54
# for fourier.kernels
constEnv$iroot.2pi=1/sqrt(2*base::pi)
constEnv$pipi=base::pi^2
constEnv$one12=1/12; constEnv$one360=1/360
constEnv$one280=1/280
constEnv$one14=1/14; constEnv$one504=1/504; constEnv$one33264=1/33264
constEnv$one18=1/18; constEnv$one792=1/792; constEnv$one61776=1/61776
constEnv$three5d486=35/486; constEnv$one528=1/528; constEnv$one37440=1/37440; constEnv$one4199040=1/4199040
constEnv$one6=1/6; constEnv$seven360=7/360; constEnv$three1d15120=31/15120
constEnv$two80d9=280/9
# for krRkernel
constEnv$halfIrootPi=.5/sqrt(base::pi)
constEnv$halfm4dpipi=.5-4/base::pi^2
constEnv$five7=5/7
constEnv$three50d429=350/429
constEnv$one75d247=175/249
constEnv$pipid6=base::pi^2/6
constEnv$pipid8=base::pi^2/8
constEnv$pipid16=base::pi^2/16
constEnv$twodpipi=2/base::pi^2
# for d2dkernel
constEnv$one5d4=15/4
constEnv$none05d16=-105/16
constEnv$none40d9=-140/9
constEnv$twodpi=2/base::pi
constEnv$npi3d16=-base::pi^3/16
# for dkernel
constEnv$pid4=base::pi/4
constEnv$pid2=base::pi/2
constEnv$onedpi=1/base::pi
# for skernel
constEnv$logpi=log(base::pi)
constEnv$fourdpi=4/base::pi
constEnv$onedpi=1/base::pi
# for mkernel
constEnv$onem8dpi=1-8/base::pi
constEnv$cos.4m=with(constEnv, (384 - 48 *pipi + pipi*pipi)/pipi/pipi)
constEnv$oneFifteenth=1/15
constEnv$thrd35=3/35
constEnv$oneSeventh=1/7
constEnv$oned21=1/21
constEnv$oned33=1/33
constEnv$three5d243=35/243
constEnv$oned22=1/22
constEnv$pipid3=base::pi^2/3
constEnv$pipid4=base::pi^2/4
constEnv$pi4.7d15=base::pi^4*7/15
constEnv$pi4.5d16=base::pi^4*5/16

## constants used by 4th order cumulants
constEnv$.order4.S2P.mat=matrix(c(
	# (transposed) matrix on page 28 of Siemiatycki (1978)
	1,		-1, -1, 0,		02, 01, 01, 01, 00, 02, 02,		-2, -6, -2, -4, 00, -4,		016, 012, 008, 010,		-48,
	0,		01, 00, 0,		-2, -2, -1, 00, 00, -4, 00,		04, 08, 02, 08, 00, 04,		-32, -28, -16, -16,		112,
	0,		00, 01, 0,		-1, 00, -1, -2, 00, 00, -4,		01, 03, 03, 02, 00, 06,		-16, -06, -08, -13,		048,
	0,		00, 00, 1,		00, -1, -1, 00, -2, 00, 00,		04, 00, 03, 02, 06, 02,		-08, -12, -14, -12,		064,
	0,		00, 00, 0,		01, 00, 00, 00, 00, 00, 00,		-1, -6, -1, -2, 00, -2,		016, 012, 008, 010,		-72,
	0,		00, 00, 0,		00, 01, 00, 00, 00, 00, 00,		-2, 00, 00, -2, 00, 00,		008, 012, 006, 004,		-48,
	0,		00, 00, 0,		00, 00, 01, 00, 00, 00, 00,		-1, 00, -2, -2, 00, -4,		016, 006, 008, 012,		-56,
	0,		00, 00, 0,		00, 00, 00, 01, 00, 00, 00,		00, 00, -1, 00, 00, 00,		000, 000, 002, 002,		-08,
	0,		00, 00, 0,		00, 00, 00, 00, 01, 00, 00,		-2, 00, -2, 00, -6, 00,		000, 006, 012, 008,		-56,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 01, 00, 00, 00,		000, 000, -04, -04,		024,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 01, 00, 00, 00, 00,		000, -02, 000, -01,		008,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		01, 00, 00, 00, 00, 00,		000, -06, -04, -04,		040,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 01,		00, 00, 00, 00, 00, -1,		002, 000, 000, 001,		-02,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 00, 00, 00,		000, 000, 000, 001,		-04,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 00, 00, 01,		-04, 000, 000, -02,		006,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 00, 00, 00,		001, 000, 000, 000,		-01,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 00, 00, 00,		000, 000, 000, 000,		001,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 01, 00,		00, 00, 00, -1, 00, 00,		004, 002, 001, 000,		-08,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 01, 00, 00,		-08, -03, -02, 000,		020,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 00, 00, 00,		000, 000, 001, 000,		-08,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 00, 01, 00,		000, 000, -01, 000,		004,
	0,		00, 00, 0,		00, 00, 00, 00, 00, 00, 00,		00, 00, 00, 00, 00, 00,		000, 001, 000, 000,		-04
	), nrow=22L, ncol=22L
)
constEnv$.order4.f.alpha=c(
	# column f (r=4 part) of table on page 25 of Siemiatycki (1978)
	8,64,48,96,96,96,192,48,192,16,12,192,16,192,96,32,48,12,32,96,48,24
)
constEnv$.order4.P2D.ord=order(
	# handwritten note on page 28 of Siemiatycki (1978)
	c(13:16,19:23,17:18,28,26:27,25,29,24,30:31,33,32,34:35)
)
constEnv$rising.fact.int.bound=c(2147483647L,46341L,1291L,216L,75L,38L,24L,18L)

safeseq =
function (from = 1L, to = 1L, by = 1L, ...) 
{
    disc = by * (from - to)
    if (disc > 0) {
        vector(class(disc), 0L)
    }
    else seq(from = from, to = to, by = by, ...)
}


midp.empirical = function(x, eps=1e-8)
{
	mean(x<x[1L]-eps) + .5* mean(x>=x[1L]-eps & x<=x[1L]+eps)
}


.factorial.rising=	function(start, nterms)
{
	stopifnot(nterms==round(nterms))
	if(nterms==0)return(do.call(paste0(c('as',class(start)), collapse='.'),list(1L)))
	if(nterms<0){
		tmp=Recall(start+nterms, -nterms)
		return(1/tmp)
	}
	if(is.numeric(start) && start==round(start) ){
		ans = prod.bigz(as.bigz(seq.int(from=start, by=1L, length.out=nterms)))
		if(is.integer(start)) {
			tryCatch(as.integer(ans), warning=function(w)ans)
		}else as.numeric(ans)
	}else if(inherits(start, c('bigz','bigq','numeric'))){
		ans = prod.bigq(as.bigq(start+seq.int(from=0, by=1, length.out=nterms)))
		do.call(paste0(c('as',class(start)), collapse='.'),list(ans))
	}else stop('unsupported class of "start"')
}
.factorial.rising=Vectorize(.factorial.rising,SIMPLIFY=FALSE)
factorial.rising=function(start, nterms)
{
	rslt = .factorial.rising(start, nterms)
	classes=sapply(rslt, class)
	if('bigq'%in%classes){
		rslt=lapply(rslt, 'as.bigq')
	}else if('bigz'%in%classes){
		rslt=lapply(rslt, 'as.bigz')
	}
	do.call('c', rslt)
}

`[[.dist`=function(x, i,j)
{
	if(i==j) return(0)
	if(i>j){tmp=i; i=j; j=tmp}
	unclass(x)[attr(x,'Size')*(i-1L) - i*(i-1L)/2L + j-i]
}

vech=function(x){
	d = NROW(x)
	stopifnot(d != NCOL(x))
	dim(x)=c(d,d)
	x[lower.tri(x, diag=TRUE)]
}		


.smooth2.kernels=c("gaussian", 'normal','triweight', 'tricube',"logistic",'sech')
.smooth1.kernels=c(.smooth2.kernels, 'biweight','quartic')
.smooth0.kernels=c(.smooth1.kernels, 'epanechnikov','parabolic','triangular','cosine')
.kernels=c(.smooth0.kernels, 'uniform', 'rectangular')

sinc =eval(substitute(function(x){
	s2=x*x; 
	ans=sin(x)/x 
	idx=which(abs(x)<1e-2) ## close to 0/0 region: taylor series
	ans[idx]=(1 - s2[idx]*oneSixth + s2[idx]*s2[idx]*one120)
	ans		
},
constEnv))
attr(sinc, 'srcref')=NULL

ilogit=make.link('logit')$linkinv
