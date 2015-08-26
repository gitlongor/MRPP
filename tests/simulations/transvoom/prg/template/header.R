base.dir=file.path('C:/Users/Long/XXXXXXXXXXXXXXXXXXXXXXXXX')
check.dir=function(d)
	if(file.exists(d)){
		if(any(!file.info(d)$isdir)) stop("file already exists")
	}else dir.create(d)
out.dir=file.path(base.dir, 'out'); check.dir(out.dir)
prg.dir=file.path(base.dir, 'prg'); check.dir(prg.dir)
dat.dir=file.path(base.dir, 'dat'); check.dir(dat.dir)
img.dir=file.path(base.dir, 'img'); check.dir(img.dir)

if(TRUE){ ## Set to FALSE if parallel running is not needed
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


}

prg.name=paste('XXXXXXXXXXXXXXXXXXXXXXXXX', Sys.Date(), sep='_')
img.name=file.path(img.dir, paste(prg.name, 'RData', sep='.'))
out.file=function(suffix, subdir) if(missing(subdir)) file.path(out.dir, paste(prg.name, suffix, sep='_')) else file.path(out.dir, subdir, paste(prg.name, suffix, sep='_')) 

setwd(out.dir)

options(contrasts=c('contr.sum', 'contr.sum'))

.sessionInfo=list(sessionInfo=sessionInfo(), Sys.info=Sys.info(),
  Sys.getenv=Sys.getenv(names=TRUE), capabilities=capabilities(),
  options=options(), RNGkind=RNGkind(), .libPaths=.libPaths(), 
  date=date(), wd=getwd(),.Random.seed=.Random.seed)


####################################### END OF HEADER ################################################
	