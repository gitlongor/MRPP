mrppBVS <-
function(y,permutedTrt, 
         importance=c('approx.keep1', 'grad.smoothp','p.grad.dist'),
         alpha.inc=NULL, alpha.exc=0, size.inc=0L, stepwise=FALSE, 
		 niter=Inf, verbose=FALSE, ...)
## y is a data matrix, with col's being variables and rows being observations
{
    if(!is.matrix(y))y=as.matrix(y)
    importance=match.arg(importance)
    if(is.null(alpha.inc)) alpha.inc=switch(importance,
		grad.smoothp = 0, 
		approx.keep1 = , 
		p.grad.dist = 0.05)
    N=nrow(y)
#    if(missing(cperm.mat)) cperm.mat=apply(permutedTrt,2,function(kk)(1:N)[-kk])
    B=nperms.permutedTrt(permutedTrt)

    ans=vector('list'); attr(ans, 'parameter')=list(importance=importance, alpha.inc=alpha.inc, alpha.exc=alpha.exc, size.inc=size.inc, stepwise=stepwise, niter=niter)
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
    repeat{
        if(verbose && (i%%verbose==0L)) {cat('iteration',i-1L,'...')
                    time0=proc.time()[3L]}
        idx=idx[imptnc<=imptnc.threshold]
        if(length(idx)==0){
            ans[[i]]=list(iter=i-1L,var.idx=integer(0L), influence=numeric(0L), 
                        p.value=numeric(0L),deleted.p.value=ans[[1L]]$p.value)
            if(verbose){
				cat('\b\b\b:\t','no variables left; mrpp.p = ; deleted.mrpp.p =',ans[[i]]$deleted.p.value,
                        ';', proc.time()[3L]-time0,'seconds passed;\n')
				message('no variables left in the inclusion set')
			}
			 attr(ans, 'status')=c(attr(ans, 'status'), ret.msg[3])
            return(returnBVS(ans))
        }
        dist0=dist(y[,idx,drop=FALSE])
		ddd=list(...); ddd$y=dist0; ddd$permutedTrt=permutedTrt
        mrpp.rslt=do.call("mrpp.test.dist", ddd)
        mrpp.stats0=mrpp.rslt$all.statistics
        imptnc=switch(importance,
			grad.smoothp=grad.smoothp(y[,idx,drop=FALSE],permutedTrt,distObj=dist0,mrpp.stats=mrpp.stats0, ...),
			approx.keep1 = p.value(grad.smoothp(y[,idx,drop=FALSE],permutedTrt,distObj=dist0,mrpp.stats=mrpp.stats0, ...), type='keep1'),
            p.grad.dist =get.p.dd.dw(y[,idx,drop=FALSE],permutedTrt,...)
		)
        var.rank=order(imptnc)
        ans[[i]]=list(iter=i-1L, var.idx=idx[var.rank], influence=imptnc[var.rank],
                      p.value=mrpp.rslt$p.value,
                      deleted.p.value=next.deleted.p)
        imptnc.threshold=if(stepwise) {tmp=max(imptnc); .5*(max(-Inf,imptnc[imptnc<tmp])+tmp)} else alpha.inc
#        idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<alpha.inc]
        if(alpha.exc>=0) {
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
        if( all(imptnc<=alpha.inc) ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[1L] )
		 if( i-1L>=niter ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[4L] )
        if( (next.deleted.p<=alpha.exc) ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[2L] )
        if( (R-length(xcl)<size.inc) ) attr(ans, 'status') = c( attr(ans, 'status') , ret.msg[3L] )
		 if (!is.null(attr(ans, 'status')))  return(returnBVS(ans))
        i=i+1L
#        if(stepwise) idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<alpha.inc]
#        if(length(idx)==0) {
#            warning('not converged'); 
#            ans[[i]]=list(iter=i-1, var.idx=numeric(0), influence=numeric(0),
#                      p.value=numeric(0),
#                      deleted.p.value=ans[[1]]$p.value)
#            return(returnBVS(ans)}
    }
}

print.mrppBVS=function(x, ...)
{
	cat("\t\tMRPP Backward Variable Selection\n")
	pars=attr(x, 'parameter')
	cat("importance measure: ", pars$importance, '\t')
	cat("minimum #{inclusion}: ", pars$size.inc, '\n')
	cat("inclusion threshold: ", pars$alpha.inc, '\t')
	cat("exclusion threshold: ", pars$alpha.exc, '\n')
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
