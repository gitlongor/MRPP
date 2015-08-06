R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != ''


library(MRPP)
set.seed(2340)
x=matrix(rnorm(10*5),10)
apply(x, 2, moment)
apply(x, 2, cumulant)
apply(x, 2, mean)
apply(x, 2, sd)
apply(x, 2, skewness)
apply(x, 2, kurtosis)



for(i in 1:5){
	if(R_CHECK_TIMINGS_ && i>3L) next 
	
	if(i==1){
		trt=gl(2,8)
	}else if (i==2){
		trt=gl(3,5)
	}else if (i==3){
		trt=gl(4,3)
	}else if (i==4){
		trt=gl(4,8)
	}else if (i==5){
		trt=gl(4,9)
	}
    x=matrix(rnorm(length(trt)*5), length(trt))
	
	mrpp.obj = mrpp(x, trt, B=if(i>3) 1e5L else Inf)
	mrpp.tst = mrpp.test(mrpp.obj)
	
	eps= if(i>3) 2*sd(mrpp.tst$all.statistics)/sqrt(1e5L) else 1e-10
	
	stopifnot(abs(
	mean(mrpp.obj) - mean(mrpp.tst$all.statistics)
	)<eps)

	stopifnot(abs(
	var(mrpp.obj) - mean((mrpp.tst$all.statistics-mean(mrpp.tst$all.statistics))^2)
	)<eps)

	stopifnot(abs(
	cumulant(mrpp.obj, order=3) - mean((mrpp.tst$all.statistics-mean(mrpp.tst$all.statistics))^3)
	)<eps)

	stopifnot(abs(
	moment(mrpp.obj, order=4, central=TRUE) - moment(mrpp.tst$all.statistics, order=4, central=TRUE)
	)<eps)

}

moment(mrpp.obj, order=1:4, central=FALSE)
moment(mrpp.obj, order=1:4, central=TRUE)
cumulant(mrpp.obj, order=1:4)

skewness(mrpp.obj)
kurtosis(mrpp.obj)
