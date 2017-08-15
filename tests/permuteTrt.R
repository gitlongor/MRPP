(R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != '')

library(MRPP)
set.seed(2340)
trt=ordered(sample(gl(2,10)), levels=c('1','2'))
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
	
	
	for(i in seq_len(50)){
		trt=sample(gl(3,3))
		pmat=permuteTrt(trt, 1e6L)
		stopifnot(ncol(pmat[[1L]]) == nparts(rep(3,3)))
		
		trt=sample(gl(4,2))
		pmat=permuteTrt(trt, 1e6L)
		stopifnot(ncol(pmat[[1L]]) == nparts(rep(2,4)))
		
		n=sample(2:5,1)
		trt=sample(rep(head(letters,n), rpois(n, 2)+1))
		n=table(trt);
		np=nparts(n);
		if(np <1e7){print(n); print(np)
			ordn = order(-n, names(n))
			trt=ordered(trt, levels=names(n)[ordn])
			part0=split(1:sum(n),trt)
			pmat=permuteTrt(trt, 1e7L)
			stopifnot(ncol(pmat[[1L]])==np)
			stopifnot(identical(part0, lapply(pmat, '[',,1L)))
			
			trt0=trt.permutedTrt(pmat)
			pmat2=permuteTrt(trt0, 1e7L)
			stopifnot(identical(pmat, pmat2))
		}
		
	}
}