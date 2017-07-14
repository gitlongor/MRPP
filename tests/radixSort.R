(R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != '')


library(MRPP)
for(i in seq(if(R_CHECK_TIMINGS_) 10L else 1e4L)){
	dat=sample(rpois(1, 30L)+1);
	sdat=sort(dat);
	d2=dat; d2[1];
	stopifnot(all(sdat==.Call(MRPP:::radixSort_prealloc, dat, d2)));
	m=tail(sdat,1L)
	stopifnot(all(sdat==.Call(MRPP:::radixSort_preallocMax, dat, d2, m)));	
}

print(.libPaths())
print(Sys.getenv('R_LIBS'))

if(!R_CHECK_TIMINGS_){
library(microbenchmark)
microbenchmark(times=10, sort(dat, method='quick'))
microbenchmark(times=10, sort(dat, method='shell'))
microbenchmark(times=10, dat[sort.list(dat, method='shell',na.last=NA)])
microbenchmark(times=10, dat[sort.list(dat, method='quick',na.last=NA)])
microbenchmark(times=10, dat[sort.list(dat, method='radix',na.last=NA)])
microbenchmark(times=10, sort.list(dat, method='radix'))
microbenchmark(.Call(MRPP:::radixSort_prealloc, dat, d2))
microbenchmark(.Call(MRPP:::radixSort_preallocMax, dat, d2, m))
microbenchmark(.Call(MRPP:::radixSort_preallocMax, dat, d2, max(dat)))
}
