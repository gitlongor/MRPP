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
norm.tumor = y$samples$norm.factors
scaled.tumor = sweep(counts.tumor, 2, y$samples$norm.factor, '/')
mean.tumor = rowMeans(scaled.tumor)
nbdisp.tumor <- y$tagwise.dispersion


################### start to simulate 
invsizes = quantile(c(nbdisp.tumor), 0:10/10)
ins.blk = cut(nbdisp.tumor, invsizes)
n.invsize = length(levels(ins.blk ))

mus = quantile(c(mean.tumor), c(0,1:9/10, 1))
mu.blk = cut(mean.tumor, mus)
n.mu=length(levels(mu.blk))

trans = c('I', 'log1p', 
	'nbinom.vst','nbinom.vst','nbinom.vst','nbinom.vst',
	'nbinom.quadLL','nbinom.quadLL',
	'nbinom.symm'
)
trans.args=list(I=NULL, log1p=NULL, 
	nbinom.vst=c(adjust='none'),  nbinom.vst=c(adjust='anscombe1'), 	nbinom.vst=c(adjust='anscombe2'), 	nbinom.vst=c(adjust='anscombe3'), 
	nbinom.quadLL = c(standardize=FALSE), nbinom.quadLL = c(standardize=TRUE), 
	nbinom.symm=NULL
)	
n.trans = length(trans)
trans.name = function(i)
{
	ans = gsub('^[^.]*\\.','',trans[i])
	if(is.null(trans.args[[i]])) return(ans)
	sprintf("%s\n%s=%s", ans, names(trans.args[[i]]), as.character(trans.args[[i]]))
}
trans.names=sapply(seq_len(n.trans), trans.name)

plot.factor=500
for(invsize.i in seq_len(n.invsize)) {
	png(sprintf('kidney.tumor.qq.invsize%d.png',invsize.i), width = plot.factor*n.trans, height = plot.factor * n.mu, pointsize=30)
		par(mfrow=c(n.mu, n.trans), mai=rep(2,4))

	for(mu.i in seq_len(n.mu)){
		this.idx = which(as.integer(ins.blk) == invsize.i & as.integer(mu.blk) == mu.i)
		if(length(this.idx)==0L){
			for(tr.i in seq_along(trans)) plot(0,0,xlab='',ylab='',main='',axes=FALSE, type='n')
			next
		}
		this.dist=as.matrix(dist(log1p(scaled.tumor[this.idx,,drop=FALSE])))
		ok.idx=this.idx[which.min(rowSums(this.dist))]
		y = scaled.tumor[ok.idx,]
		this.mean = mean.tumor[ok.idx]
		this.nbdisp=nbdisp.tumor[ok.idx]
	
		for(tr.i in seq_along(trans)){
			lst = c(list(`_data`=y, method=trans[tr.i], invsize=invsizes[invsize.i]), as.list( trans.args[[tr.i]]))
			transy=do.call('transform', lst)
			qqnorm(transy, sub=sprintf('mu=%.2f; invsize=%.2f', this.mean, this.nbdisp),
				xlab='N(0,1) quantile', ylab=sprintf('%s quantile', trans[tr.i]),  
				main=trans.names[tr.i])
			qqline(transy)
		}
	}
	dev.off()
}


### Subset down to just non-tumor counts for dispersion estimation in non-nontumor group
counts.nontumor <- counts[ , treatment == "Non-Tumor"]

### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts.nontumor)
y <- calcNormFactors(y)
y <- estimateDisp(y)
norm.nontumor = y$samples$norm.factors
scaled.nontumor = sweep(counts.nontumor, 2, y$samples$norm.factor, '/')
mean.nontumor = rowMeans(scaled.nontumor)
nbdisp.nontumor <- y$tagwise.dispersion



################### start to simulate 
invsizes = quantile(c(nbdisp.nontumor), 0:10/10)
ins.blk = cut(nbdisp.nontumor, invsizes)
n.invsize = length(levels(ins.blk ))

mus = quantile(c(mean.nontumor), c(0,1:9/10, 1))
mu.blk = cut(mean.nontumor, mus)
n.mu=length(levels(mu.blk))

trans = c('I', 'log1p', 
	'nbinom.vst','nbinom.vst','nbinom.vst','nbinom.vst',
	'nbinom.quadLL','nbinom.quadLL',
	'nbinom.symm'
)
trans.args=list(I=NULL, log1p=NULL, 
	nbinom.vst=c(adjust='none'),  nbinom.vst=c(adjust='anscombe1'), 	nbinom.vst=c(adjust='anscombe2'), 	nbinom.vst=c(adjust='anscombe3'), 
	nbinom.quadLL = c(standardize=FALSE), nbinom.quadLL = c(standardize=TRUE), 
	nbinom.symm=NULL
)	
n.trans = length(trans)
trans.name = function(i)
{
	ans = gsub('^[^.]*\\.','',trans[i])
	if(is.null(trans.args[[i]])) return(ans)
	sprintf("%s\n%s=%s", ans, names(trans.args[[i]]), as.character(trans.args[[i]]))
}
trans.names=sapply(seq_len(n.trans), trans.name)

plot.factor=500
for(invsize.i in seq_len(n.invsize)) {
	png(sprintf('kidney.nontumor.qq.invsize%d.png',invsize.i), width = plot.factor*n.trans, height = plot.factor * n.mu, pointsize=30)
		par(mfrow=c(n.mu, n.trans), mai=rep(2,4))

	for(mu.i in seq_len(n.mu)){
		this.idx = which(as.integer(ins.blk) == invsize.i & as.integer(mu.blk) == mu.i)
		if(length(this.idx)==0L){
			for(tr.i in seq_along(trans)) plot(0,0,xlab='',ylab='',main='',axes=FALSE, type='n')
			next
		}
		this.dist=as.matrix(dist(log1p(scaled.nontumor[this.idx,,drop=FALSE])))
		ok.idx=this.idx[which.min(rowSums(this.dist))]
		y = scaled.nontumor[ok.idx,]
		this.mean = mean.nontumor[ok.idx]
		this.nbdisp=nbdisp.nontumor[ok.idx]
	
		for(tr.i in seq_along(trans)){
			lst = c(list(`_data`=y, method=trans[tr.i], invsize=invsizes[invsize.i]), as.list( trans.args[[tr.i]]))
			transy=do.call('transform', lst)
			qqnorm(transy, sub=sprintf('mu=%.2f; invsize=%.2f', this.mean, this.nbdisp),
				xlab='N(0,1) quantile', ylab=sprintf('%s quantile', trans[tr.i]),  
				main=trans.names[tr.i])
			qqline(transy)
		}
	}
	dev.off()
}


rm(counts, kidney)

save.image('SimSeqTransformation.RData', compress='bzip2')
