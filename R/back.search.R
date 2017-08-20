mrppBVS <-
function(y,permutedTrt, 
         importance=c('grad.smoothp','grad.energy','p.grad.dist','approx.keep1'),
         inc.thresh=NULL, exc.thresh=0, size.inc=0L, stepwise=FALSE, 
		 niter=Inf, verbose=FALSE, ...)
## y is a data matrix, with col's being variables and rows being observations
{
    if(!is.matrix(y))y=as.matrix(y)
    importance=match.arg(importance)
    if(is.null(inc.thresh)) inc.thresh=switch(importance,
		grad.smoothp = , 
		grad.energy = 0, 
		approx.keep1 = , 
		p.grad.dist = 0.05)
    N=nrow(y)
#    if(missing(cperm.mat)) cperm.mat=apply(permutedTrt,2,function(kk)(1:N)[-kk])
    B=nperms.permutedTrt(permutedTrt)
	ddd.mrpp=list(...)
	ddd.mrpp$y=y
	ddd.mrpp$permutedTrt=permutedTrt
	ddd.mrpp$trt=NULL
	ddd.mrpp$B=NULL
	ddd.mrpp.idx=names(ddd.mrpp)%in%names(formals(mrpp.matrix))
	mrpp.obj=do.call(mrpp, ddd.mrpp[ddd.mrpp.idx])

    ans=vector('list'); attr(ans, 'parameter')=list(importance=importance, inc.thresh=inc.thresh, exc.thresh=exc.thresh, size.inc=size.inc, stepwise=stepwise, niter=niter)
    R=NCOL(y)
    idx=1:R  # inclusion set
    xcl=integer(0L)
    i=1L
    imptnc.threshold=Inf
    imptnc=rep(-Inf, R)
    next.deleted.p=numeric(0L)
	ret.msg=c('All remaining variables meet importance requirements',
			  'Deleted variables fail to meet unimportance requirements if deletion continues', 
			  'Minimum number of selected variables reached',
			  'Maximum number of iterations reached')
	returnBVS=function(ret)structure(ret, class='mrppBVS')
	
	ddd.mrppt=list(...); ddd.mrppt$y=mrpp.obj; ddd.mrppt$method='permutation'
	ddd.mrppt.idx=names(ddd.mrppt)%in%names(formals(mrpp.test.mrpp))
	ddd.mrppt=ddd.mrppt[ddd.mrppt.idx]
	var.sign = var.rank = numeric(R)

	ddd.grad=list(...); 
	ddd.grad$y=matrix(NA_real_,0L,0L)
	ddd.grad$permutedTrt=permutedTrt
	ddd.grad$distObj=mrpp.obj$distObj
	ddd.grad$mrpp.stats=numeric(0L)
	if(importance=='grad.energy') {
		ddd.grad$bw=Inf
		ddd.grad$adjust='weighted.mean'
	}else ddd.grad$adjust='log scale'
	ddd.grad.idx=names(ddd.grad)%in%names(formals(grad.smoothp))
	ddd.grad=ddd.grad[ddd.grad.idx]
	
    repeat{
        if(verbose && (i%%verbose==0L)) {cat('iteration',i-1L,'...')
                    time0=proc.time()[3L]}
		var.sign[idx]=-1*(imptnc<=inc.thresh); var.sign[idx[imptnc==inc.thresh]]=0
		var.rank[idx]=rank(imptnc, ties.method='average')
        idx=idx[imptnc<=imptnc.threshold]
		
        if(length(idx)==0){
            ans[[i]]=list(iter=i-1L,var.idx=integer(0L), influence=numeric(0L), 
                        p.value=numeric(0L),deleted.p.value=ans[[1L]]$p.value, var.sign=var.sign, var.rank=var.rank)
            if(verbose){
				cat('\b\b\b:\t','no variables left; mrpp.p = ; deleted.mrpp.p =',ans[[i]]$deleted.p.value,
                        ';', proc.time()[3L]-time0,'seconds passed;\n')
				message('no variables left in the inclusion set')
			}
			attr(ans, 'status')=c(attr(ans, 'status'), ret.msg[3])
            return(returnBVS(ans))
        }
        dist0=dist(y[,idx,drop=FALSE])
		#mrpp.obj$distObj=dist0; mrpp.obj$R=length(idx)
		ddd.mrppt$y$distObj=dist0; ddd.mrppt$y$R=length(idx)
        mrpp.rslt=do.call(mrpp.test.mrpp, ddd.mrppt)
        mrpp.stats0=mrpp.rslt$all.statistics
        imptnc=switch(importance,
			grad.smoothp=, 
			grad.energy=,
			approx.keep1 ={
				ddd.grad$y=y[,idx,drop=FALSE]
				ddd.grad$distObj=dist0
				ddd.grad$mrpp.stats=mrpp.stats0
				do.call(grad.smoothp, ddd.grad)
			},
            p.grad.dist =get.p.dd.dw(y[,idx,drop=FALSE],permutedTrt,...) # CHECKME
		)
		if(importance=='approx.keep1') imptnc=p.value(imptnc, type='keep1')
        var.ord=order(imptnc)
        ans[[i]]=list(iter=i-1L, var.idx=idx[var.ord], influence=imptnc[var.ord],
                      p.value=mrpp.rslt$p.value,
                      deleted.p.value=next.deleted.p, var.sign=var.sign, var.rank=var.rank)
        imptnc.threshold=if(stepwise) {tmp=max(imptnc); .5*(max(-Inf,imptnc[imptnc<tmp])+tmp)} else inc.thresh
#        idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<inc.thresh]
        if(exc.thresh>=0) {
            xcl=c(xcl, idx[imptnc>imptnc.threshold])
            dist.del=dist(y[,xcl,drop=FALSE])
            if(all(!is.na(dist.del)) && length(xcl)>0L)
#                ans[[i]]$deleted.p.value=mrpp.test.dist(dist.del, permutedTrt=permutedTrt)$p.value
                next.deleted.p=mrpp.test.dist(dist.del, permutedTrt=permutedTrt)$p.value
        }
        if(verbose && (i%%verbose==0L)) {
          cat('\b\b\b:\t',length(idx),'variable(s) left; mrpp.p =',ans[[i]]$p.value,';', 
                        'deleted.mrpp.p =',ans[[i]]$deleted.p.value,
                        ';', proc.time()[3L]-time0,'seconds passed;\n')
        }
        if( all(imptnc<=inc.thresh) ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[1L] )
		 if( i-1L>=niter ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[4L] )
        if( (next.deleted.p<=exc.thresh) ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[2L] )
        if( (R-length(xcl)<size.inc) ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[3L] )
		 if (!is.null(attr(ans, 'status')))  return(returnBVS(ans))
        i=i+1L
    }
	stop('unexpected exit')
}

print.mrppBVS=function(x, ...)
{
	cat("\t\tMRPP Backward Variable Selection\n")
	pars=attr(x, 'parameter')
	cat("importance measure: ", pars$importance, '\t')
	cat("minimum #{inclusion}: ", pars$size.inc, '\n')
	cat("inclusion threshold: ", pars$inc.thresh, '\t')
	cat("exclusion threshold: ", pars$exc.thresh, '\n')
	cat("Iteration History:\n")
	for(i in seq_along(x)){
		cat('\titeration',x[[i]]$iter,'\t', length(x[[i]]$var.idx),'variable(s) left; mrpp.p =',x[[i]]$p.value,';', 
                        'deleted.mrpp.p =',x[[i]]$deleted.p.value,
                        '\n')
	}
	conv=attr(x, 'status')
	cat("Convergence Criteri", if(length(conv)==1L) 'on' else 'a', ":\n",sep='')
	for(cc in conv)cat("\t",cc,"\n",sep='')
	invisible(x)
}
