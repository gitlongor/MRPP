R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != ''
if(R_CHECK_TIMINGS_) q('no')


library(edgeR)
library(SimSeq)
library(MRPP)

### Load Data
RNGkind("Mersenne-Twister", "Inversion")
set.seed(156165L)
data(kidney); 
kidney$counts=kidney$counts[sample(nrow(kidney$counts)),]
counts <- kidney$counts
treatment <- kidney$treatment

### Set simulation variables
filter.mean <- 0 # lower bound of average read count for simulated genes
filter.nonzero <- 1 # lower bound for nonzero read counts for simulated genes

### Remove low count genes
keep.counts <- which(( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero ))
counts <- counts[keep.counts, ]

### Subset down to just tumor counts for dispersion estimation in tumor group
counts.tumor <- counts[ , treatment == "Tumor"]

### Use edgeR package's tagwise dispersion estimate
### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts.tumor)
y <- calcNormFactors(y)
y <- estimateDisp(y)
nbdisp.tumor <- y$tagwise.dispersion
norm.tumor = y$samples$norm.factors
mean.tumor = rowMeans(sweep(counts.tumor, 2L, norm.tumor, '/'))

### Subset down to just non-tumor counts for dispersion estimation in non-tumor group
counts.nontumor <- counts[ , treatment == "Non-Tumor"]

### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts.nontumor)
y <- calcNormFactors(y)
y <- estimateDisp(y)
nbdisp.nontumor <- y$tagwise.dispersion
mean.nontumor = rowMeans(sweep(counts.nontumor, 2L, y$samples$norm.factors, '/'))



################### start to simulate 
taus = quantile(c(nbdisp.tumor, nbdisp.nontumor), 0:10/10)
mus = quantile(c(mean.tumor, mean.nontumor), c(.01,1:9/10, .99))
n.mu=length(mus)
n.tau = length(taus)
reps=1e4L
N=1L
trans = c('I', 'log1p', 'nbinom.vst','nbinom.quadLL','nbinom.symm')
#trans.nargs=c(1L,1L, 2L ,2L, 2L)
n.trans = length(trans)

plot.factor=500
for(tau.i in seq_along(taus)) {
	png(sprintf('kidney.qq.tau%d.png',tau.i), width = plot.factor*n.trans, height = plot.factor * n.mu, pointsize=30)
		par(mfcol=c(n.mu, n.trans), mai=rep(2,4))
	y=rnbinom(reps*n.mu, size=1/taus[tau.i], mu=mus)
	dim(y)=c(n.mu, reps)
	
	for(tr.i in seq_along(trans)){
		#if(trans.nargs[tr.i]==1L){
		#	transy=MRPP::transform(y, method=trans[tr.i])
		#}else {
			transy=MRPP::transform(y, method=trans[tr.i], tau=taus[tau.i])
		#}
		for(mu.i in seq_along(mus)){
			qqnorm(transy[mu.i,], main=sprintf('mu=%.2f; tau=%.2f', mus[mu.i], taus[tau.i]),
				xlab='N(0,1) quantile', ylab=sprintf('%s quantile', trans[tr.i]), pch='.',cex=3)
			qqline(transy[mu.i,])
		}
	}
	dev.off()
}


rm(y,transy, counts, kidney)

save.image('NBtransformation.RData', compress='bzip2')
