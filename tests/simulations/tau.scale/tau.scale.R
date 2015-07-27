library(MRPP)
reps=1e4L
R=2L; n=5L; ntrt=2L; N=ntrt*n 
mu1=rep(-1:0, c(n, N-n))
trt=gl(ntrt, n)
pmat=permuteTrt(trt,B=Inf)
(B = nperms.permutedTrt(pmat))
mh.trt=gl(ntrt, reps)
mh.blk=factor(rep(seq(reps), ntrt))

simTau=function(sds=c(1,1), standardize=FALSE){
	y=cbind(rnorm(N, mu1, sds[1L]), replicate(R-1L, rnorm(N,,sds[-1L])))
	dists = lapply(seq(NCOL(y)), function(x) dist(y[,x])^2/(dist(y)+dist(y[,-x])) )
	sapply(dists, function(dd)scale(mrpp.test(dd, permutedTrt=pmat)$all.statistics, scale=standardize)[1L])
}

RNGkind("Mersenne-Twister", "Inversion")
set.seed(2016565L)

seed.bak=.Random.seed
taus.1.1=t(replicate(reps, simTau(c(1,1))))
.Random.seed=seed.bak
taus.1.1std=t(replicate(reps, simTau(c(1,1), standardize=TRUE)))
mantelhaen.test(factor(c(taus.1.1[,1]<taus.1.1[,2], taus.1.1std[,1]<taus.1.1std[,2])), mh.trt, mh.blk)

seed.bak=.Random.seed
taus.1.10=t(replicate(reps, simTau(c(1,10))))
.Random.seed=seed.bak
taus.1.10std=t(replicate(reps, simTau(c(1,10), standardize=TRUE)))
mantelhaen.test(factor(c(taus.1.10[,1]<taus.1.10[,2], taus.1.10std[,1]<taus.1.10std[,2])), mh.trt, mh.blk)

seed.bak=.Random.seed
taus.1.100=t(replicate(reps, simTau(c(1,100))))
.Random.seed=seed.bak
taus.1.100std=t(replicate(reps, simTau(c(1,100), standardize=TRUE)))
mantelhaen.test(factor(c(taus.1.100[,1]<taus.1.100[,2], taus.1.100std[,1]<taus.1.100std[,2])), mh.trt, mh.blk)

seed.bak=.Random.seed
taus.1.01=t(replicate(reps, simTau(c(1,.1))))
.Random.seed=seed.bak
taus.1.01std=t(replicate(reps, simTau(c(1,.1), standardize=TRUE)))
mantelhaen.test(factor(c(taus.1.01[,1]<taus.1.01[,2], taus.1.01std[,1]<taus.1.01std[,2])), mh.trt, mh.blk)

par(mfcol=c(2,4))
bx=function(z){boxplot(z, main=sprintf('%%good=%.3f', mean(z[,1]<z[,2])));invisible(NULL)}
sapply(list(taus.1.01, taus.1.01std, taus.1.1, taus.1.1std, taus.1.10, taus.1.10std, taus.1.100, taus.1.100std), bx)

save.image('tau.scale.RData')
q('no')
