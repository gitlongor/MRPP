(R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != '')

library(MRPP)
trt=ordered(sample(gl(2,10)), levels=c('1','2'))
nparts(table(trt)) == 92378L		## 92378 partitions =  choose(20,10)/2


set.seed(1032940L)
urand.bigz(0,seed=1032940L) # init seed
pmat=permuteTrt(trt, 1e3L)		## use 1000 random permutations
if(!R_CHECK_TIMINGS_){

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