R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != ''

require(MRPP)
set.seed(2340)
trt=gl(2,10)
nparts(table(trt)) == 92378L		## 92378 partitions =  choose(20,10)/2

urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 1e3L)		## use 1000 random permutations
if(!R_CHECK_TIMINGS_){
	urand.bigz(0,seed=1032940L) # init seed
	pmat1=permuteTrt(trt, 1e3L)		## use 1000 random permutations
	stopifnot(identical(pmat, pmat1))


	urand.bigz(0,seed=1032940L) # init seed
	pmat.ind=permuteTrt(trt, 1e3L, idxOnly=TRUE)		## use 1000 random permutations

	stopifnot(ncol(pmat[[1L]]) == length(attr(pmat.ind, 'idx')) )

	allmats = sapply(attr(pmat.ind, 'idx'), function(x)dec2permvec(as.bigz(x), N=length(trt)))
	allmats = apply(allmats, 2L, function(x)unlist(by(x, trt, sort)))

	stopifnot(all(allmats == do.call('rbind', pmat)))


	pmat=permuteTrt(trt, 1e6L)		## use all partitions, as 1e6L >= 92378 
	stopifnot(ncol(pmat[[1L]]) == 92378L)
}