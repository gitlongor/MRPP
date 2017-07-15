constEnv=new.env()
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
		ans = prod.bigz(as.bigz(seq(from=start, length=nterms)))
		if(is.integer(start)) {
			tryCatch(as.integer(ans), warning=function(w)ans)
		}else as.numeric(ans)
	}else if(inherits(start, c('bigz','bigq','numeric'))){
		ans = prod.bigq(as.bigq(start+seq(from=0, length=nterms)))
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

.kernels=c("gaussian", "uniform", "rectangular", "epanechnikov", "triangular", "biweight", 'triweight', 'tricube',"logistic","cosine")
.smooth.kernels=c("gaussian", "biweight", 'triweight', 'tricube',"logistic","cosine")

sinc =function(x){
	s2=x*x; 
	ans=sin(x)/x 
	idx=which(abs(x)<1e-2) ## close to 0/0 region: taylor series
	ans[idx]=(1 - s2[idx]*oneSixth + s2[idx]*s2[idx]*one120)
	ans		
}
environment(sinc)=constEnv