base.dir=file.path('D:/work/dev/Rstudio/MRPP/tests/simulations/transvoom')
check.dir=function(d)
	if(file.exists(d)){
		if(any(!file.info(d)$isdir)) stop("file already exists")
	}else dir.create(d)
out.dir=file.path(base.dir, 'out'); check.dir(out.dir)
prg.dir=file.path(base.dir, 'prg'); check.dir(prg.dir)
dat.dir=file.path(base.dir, 'dat'); check.dir(dat.dir)
img.dir=file.path(base.dir, 'img'); check.dir(img.dir)

if(FALSE){
	library(parallel)
	RNGkind("L'Ecuyer-CMRG")
	cl=makeCluster(8L)
	clusterEvalQ(cl, {
		RNGkind("L'Ecuyer-CMRG"); 
		options(contrasts=c('contr.sum', 'contr.sum'))
	})
	clusterSetRNGStream(cl , 9992722L)
}
RNGkind("Mersenne-Twister", "Inversion")
set.seed(9992722L)

{### attaching custom packages here

### Load Bioconductor Packages
require(edgeR)
### Load CRAN Packages
require(SimSeq)

}

prg.name=paste('NB_disp_estimates', Sys.Date(), sep='_')
img.name=file.path(img.dir, paste(prg.name, 'RData', sep='.'))
out.file=function(suffix, subdir) if(missing(subdir)) file.path(out.dir, paste(prg.name, suffix, sep='_')) else file.path(out.dir, subdir, paste(prg.name, suffix, sep='_')) 

setwd(out.dir)

options(contrasts=c('contr.sum', 'contr.sum'))

.sessionInfo=list(sessionInfo=sessionInfo(), Sys.info=Sys.info(),
  Sys.getenv=Sys.getenv(names=TRUE), capabilities=capabilities(),
  options=options(), RNGkind=RNGkind(), .libPaths=.libPaths(), 
  date=date(), wd=getwd(),.Random.seed=.Random.seed)


####################################### END OF HEADER ################################################
	

### Load Data
data(kidney); 
kidney$counts=kidney$counts[sample(nrow(kidney$counts)),]
counts <- kidney$counts
treatment <- kidney$treatment

### Set simulation variables
filter.mean <- 0 # lower bound of average read count for simulated genes
filter.nonzero <- 1 # lower bound for nonzero read counts for simulated genes

### Remove low count genes
keep.counts <- ( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero )
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

### Subset down to just non-tumor counts for dispersion estimation in non-tumor group
counts.nontumor <- counts[ , treatment == "Non-Tumor"]

### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts.nontumor)
y <- calcNormFactors(y)
y <- estimateDisp(y)
nbdisp.nontumor <- y$tagwise.dispersion

### Save Results
saveRDS(nbdisp.tumor, file = file.path(out.dir, "nbdisp_tumor_edgeR.RDS"))
saveRDS(nbdisp.nontumor, file = file.path(out.dir, "nbdisp_nontumor_edgeR.RDS"))


save.image(img.name)
