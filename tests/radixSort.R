R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != ''


library(MRPP)
for(i in seq(1e4L)){
	dat=sample(rpois(1, 30L)+1);
	sdat=sort(dat);
	d2=dat; d2[1];
	stopifnot(all(sdat==.Call('radixSort_prealloc', dat, d2)));
	m=tail(sdat,1L)
	stopifnot(all(sdat==.Call('radixSort_preallocMax', dat, d2, m)));	
}

print(.libPaths())
print(Sys.getenv('R_LIBS'))
library(microbenchmark)
microbenchmark(sort(dat, method='quick'))
microbenchmark(sort(dat, method='shell'))
microbenchmark(dat[sort.list(dat, method='shell',na.last=NA)])
microbenchmark(dat[sort.list(dat, method='quick',na.last=NA)])
microbenchmark(dat[sort.list(dat, method='radix',na.last=NA)])
microbenchmark(sort.list(dat, method='radix'))
microbenchmark(.Call('radixSort_prealloc', dat, d2))
microbenchmark(.Call('radixSort_preallocMax', dat, d2, m))
microbenchmark(.Call('radixSort_preallocMax', dat, d2, max(dat)))
