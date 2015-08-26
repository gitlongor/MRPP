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

{
### Load Bioconductor Packages
library(DESeq2)
library(edgeR)
library(limma)
### Load CRAN Packages
library(QuasiSeq)
#library(samr) ## this package change options('error')
library(fdrtool)
library(SimSeq)
library(MRPP)
}


prg.name=paste('transvoom.n5', Sys.Date(), sep='_')
img.name=file.path(img.dir, paste(prg.name, 'RData', sep='.'))
out.file=function(suffix, subdir) if(missing(subdir)) file.path(out.dir, paste(prg.name, suffix, sep='_')) else file.path(out.dir, subdir, paste(prg.name, suffix, sep='_')) 

setwd(out.dir)

options(contrasts=c('contr.sum', 'contr.sum'))

.sessionInfo=list(sessionInfo=sessionInfo(), Sys.info=Sys.info(),
  Sys.getenv=Sys.getenv(names=TRUE), capabilities=capabilities(),
  options=options(), RNGkind=RNGkind(), .libPaths=.libPaths(), 
  date=date(), wd=getwd(),.Random.seed=.Random.seed)

####################################### END OF HEADER ################################################
### rnbinom parameterization Var(Y) = mu + mu^2/size
### edgeR parameterization   Var(Y) = mu + mu^2 * invsize
nb.disp.tumor <- 1 / readRDS(file.path(out.dir, "nbdisp_tumor_edgeR.RDS"))
nb.disp.nontumor <- 1 / readRDS(file.path(out.dir, "nbdisp_nontumor_edgeR.RDS"))


## list of transformations and associated arguments
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
trans.names=sapply(seq(n.trans), trans.name)


### Set simulation variables
n.iter <- 200     # Number of iterations
k.ind <- 5    # Sample size in each simulated treatment group 
n.genes <- 8500  # No of genes in each simulated matrix
n.diff <- 2000   # No of DE genes in each simulated matrix
n.genes.trim <- 5000 # No of genes to trim down to
pi0=.8
n.diff.trim <- round(n.genes.trim*(1-pi0)) # No of DE genes to trim down to
pi0=n.diff.trim/n.genes.trim
filter.mean <- 0 # lower bound of average read count for simulated genes
filter.nonzero <- 1 # lower bound for nonzero read counts for simulated genes

trt=gl(2, k.ind)
design <- model.matrix(~trt)

### Load Data
data(kidney); 
kidney$counts=kidney$counts[sample(nrow(kidney$counts)),] ## not truly necessary

counts <- kidney$counts
tumor <- kidney$treatment
replic <- kidney$replic

### Remove low count genes
keep.counts <- ( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero )
counts <- counts[keep.counts, ]

### Preprocessing Steps to speed up SimData function from SimSeq package

### Compute normalization factors to use in SimData function
### Effective library size is product of library size and size factors
### from calcNormFactors
lib.sizes <- apply(counts, 2, sum)
nf <- calcNormFactors(counts) * lib.sizes

### Compute weights to sample DE genes in SimData function
probs <- CalcPvalWilcox(counts, treatment = tumor, replic = replic,
                        sort.method = "paired", sorted = TRUE, nf,
                        exact = FALSE)
wghts <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr

### Initialize matrix of p-value output for each statistical method
pvals=array(NA_real_, dim=c(n.genes.trim, n.iter, 2L, n.trans))
dimnames(pvals)=list(seq(n.genes.trim), seq(n.iter), c('simseq','nb'), trans.names)


### Preprocessing steps for simulating NB data
lambdas <- matrix(NA, nrow = nrow(counts), ncol = 2)

sum.nf.nontumor <- sum(nf[tumor == "Non-Tumor"])
sum.nf.tumor <- sum(nf[tumor == "Tumor"])
lambdas[, 1] <- rowSums(counts[, tumor == "Non-Tumor"])/sum.nf.nontumor
lambdas[, 2] <- rowSums(counts[, tumor == "Tumor"])/sum.nf.tumor




for(i in 1:n.iter){
  
    {### Simulate matrix of read counts from SimSeq
   
    counts.simseq.list <- SimData(counts = counts, treatment = tumor, replic = replic, 
                         sort.method = "paired", k.ind = k.ind, n.genes = n.genes,
                         n.diff = n.diff, norm.factors = nf, weights = wghts, switch.trt = TRUE)
    
    counts.simseq <- counts.simseq.list$counts # Simulated Count matrix from SimSeq
    genes.samp <- counts.simseq.list$genes.subset # Genes sampled from source matrix
    de.genes <- counts.simseq.list$DE.genes # DE genes sampled from source matrix
    ee.genes <- genes.samp[ ! (genes.samp %in% de.genes) ] # EE genes sampled from source matrix
    samp.col <- counts.simseq.list$col # Columns sampled in SimSeq algorithm
    de.genes.sim <- counts.simseq.list$genes.subset %in% de.genes # logical vector giving which genes are DE in simulted matrix
    
    ### Compute matrix of means for simulating from NB model
    mu.samp <- matrix(NA, nrow = n.genes, ncol = 2 * k.ind)
    nf.samp <- nf[samp.col]
    
    ### Use normalization factors from Tumor group to match SimSeq algorithm
    mu.samp[de.genes.sim, 1:k.ind] <- lambdas[de.genes, 1, drop = FALSE] %*% nf.samp[ (k.ind + 1):(2 * k.ind) ]
    mu.samp[de.genes.sim, (k.ind + 1):(2 * k.ind)] <- lambdas[de.genes, 2, drop = FALSE] %*% nf.samp[ (2 * k.ind + 1):(3 * k.ind) ]
    mu.samp[ !de.genes.sim, ] <-  lambdas[ee.genes, 2, drop = FALSE] %*% nf.samp[ (k.ind + 1):(3 * k.ind) ]
  
    ### Set dispersion estimates
    disp.samp <- matrix(NA, nrow = n.genes, ncol = 2)
    disp.samp[!de.genes.sim, 1] <- disp.samp[!de.genes.sim, 2] <- nb.disp.tumor[ee.genes]
    disp.samp[de.genes.sim, 1] <- nb.disp.nontumor[de.genes]
    disp.samp[de.genes.sim, 2] <- nb.disp.tumor[de.genes]
    
    ### Simulate matrix of read counts from NB model
    counts.nb <- matrix(NA, nrow = n.genes, ncol = 2 * k.ind)
    for(jj in 1:n.genes){
      counts.nb[jj, 1:k.ind] <- rnbinom(k.ind, size = disp.samp[jj, 1], mu = mu.samp[jj, 1:k.ind])
      counts.nb[jj, -(1:k.ind)] <- rnbinom(k.ind, size = disp.samp[jj, 2], mu = mu.samp[jj, -(1:k.ind)])
    }
    
    ### Apply filtering rules to both simulated datasets and only keep genes who pass both filters
    keep.genes.simseq <- ( rowMeans(counts.simseq) >= filter.mean ) & ( rowSums(counts.simseq > 0) >= filter.nonzero )
    keep.genes.nb <- ( rowMeans(counts.nb) >= filter.mean ) & ( rowSums(counts.nb > 0) >= filter.nonzero )
    keep <- keep.genes.simseq & keep.genes.nb
    
    ee.genes <- sample(which(!de.genes.sim & keep), n.genes.trim - n.diff.trim)
    de.genes <- sample(which(de.genes.sim & keep), n.diff.trim)
    counts.simseq <- counts.simseq[c(ee.genes, de.genes), ]
    counts.nb <- counts.nb[c(ee.genes, de.genes), ]
	
  } # end of sim dat
  
    ### Compute Normalization Factors
    nf.simseq <- calcNormFactors(counts.simseq)
    nf.nb <- calcNormFactors(counts.nb)
	
	### scaling count data
	scaled.simseq=ceiling(sweep(counts.simseq, 2, nf.simseq, '/'))
	scaled.nb=ceiling(sweep(counts.nb, 2, nf.nb, '/'))
	
	for(sim.type in c('simseq','nb')){
		this.dat = get(paste0(c('scaled',sim.type),collapse='.'))
		tmp <- DGEList(counts = this.dat)
		tmp <- calcNormFactors(tmp)
		tmp <- estimateDisp(tmp)
		nbdisp <- tmp$tagwise.dispersion
		for(tr.i in seq_along(trans)){
			lst = c(list(`_data`=this.dat, method=trans[tr.i], invsize=nbdisp), as.list( trans.args[[tr.i]]))
			transy=do.call('transform', lst)
			
			lib.size=rep(max(colSums(transy)), length=ncol(transy))
			 ### fit voom
			 v <- voom(transy, design=design, lib.size=lib.size, normalize.method='none', plot=FALSE, span=.5)
			 fit <- lmFit(v,design)
			 fit <- eBayes(fit, proportion = 1-pi0, stdev.coef.lim=c(0,Inf))
			 pvals[, i, sim.type, tr.i] <- fit$p.value[, 2]
		}
	}
}

### Save Results

mainDir <- getwd()
subDir <- "pval_output"
if( !file.exists(file.path(mainDir, subDir)) )
{
  dir.create(file.path(mainDir, subDir)) 
}
if( !file.exists(file.path(mainDir, subDir, paste0("ss", k.ind))) )
{
  dir.create(file.path(mainDir, subDir, paste0("ss", k.ind))) 
}
#saveRDS(pvals.deseq2.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_simseq.RDS"))
#saveRDS(pvals.quasiseq.simseq.ql, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_simseq_ql.RDS"))
#saveRDS(pvals.quasiseq.simseq.spline, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_simseq_spline.RDS"))
#saveRDS(pvals.quasiseq.simseq.bart, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_simseq_bart.RDS"))
#saveRDS(pvals.edger.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_simseq.RDS"))
#saveRDS(pvals.samseq.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_simseq.RDS"))
#saveRDS(pvals.voom.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_simseq.RDS"))

saveRDS(fold.quasiseq.simseq.fit1, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit1.RDS"))
saveRDS(fold.quasiseq.simseq.fit2, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit2.RDS"))
saveRDS(fold.quasiseq.simseq.fit3, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit3.RDS"))
saveRDS(fold.voom.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_simseq.RDS"))

saveRDS(fold.quasiseq.nb.fit1, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit1.RDS"))
saveRDS(fold.quasiseq.nb.fit2, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit2.RDS"))
saveRDS(fold.quasiseq.nb.fit3, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit3.RDS"))
saveRDS(fold.voom.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_nb.RDS"))

saveRDS(fold.trulyDEgeneList, file.path(mainDir, subDir, paste0("ss", k.ind), "fold_trulyDEgeneList.RDS"))

#saveRDS(pvals.deseq2.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_nb.RDS"))
#saveRDS(pvals.quasiseq.nb.ql, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_nb_ql.RDS"))
#saveRDS(pvals.quasiseq.nb.spline, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_nb_spline.RDS"))
#saveRDS(pvals.quasiseq.nb.bart, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_nb_bart.RDS"))
#saveRDS(pvals.edger.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_nb.RDS"))
#saveRDS(pvals.samseq.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_nb.RDS"))
#saveRDS(pvals.voom.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_nb.RDS"))


save.image('simulation_output_n5.RData', compress='xz')
q('no')