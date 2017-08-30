mrppBVS <-
function(y, #permutedTrt, 
         importance=c('grad.smoothp','grad.energy','p.grad.dist','approx.keep1'),
         inc.thresh=NULL, exc.thresh=0, size.inc=0L, stepwise=FALSE, 
		 niter=Inf, verbose=FALSE, ...)
## y is an mrpp obj
{
	y$data.env$y.bakBVS=y$y
	on.exit({ # recover original data.env
			if(!is.null(y$data.env$y.bakBVS)){
				y$data.env$y=y$data.env$y.bakBVS
				rm(y.bakBVS, envir=y$data.env)
			}
	})
#print(proc.time())
    importance=match.arg(importance)
    if(is.null(inc.thresh)) inc.thresh=switch(importance,
		grad.smoothp = 1, 
		grad.energy = 0, 
		approx.keep1 = , 
		p.grad.dist = 0.05)
#    N=y$nobs

#    B=nperms.permutedTrt(permutedTrt)
#	ddd.mrpp=list(...)
#	ddd.mrpp$y=y
#	ddd.mrpp$permutedTrt=permutedTrt
#	ddd.mrpp$trt=NULL
#	ddd.mrpp$B=NULL
#	ddd.mrpp.idx=names(ddd.mrpp)%in%names(formals(mrpp.matrix))
#	mrpp.obj=do.call(mrpp, ddd.mrpp[ddd.mrpp.idx])
	distFunc=y$distFunc

    R=y$R
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
	
	ddd=list(...)
	
	frm.mrppt=formals(mrpp.test.mrpp)
	test.method=if(is.null(ddd$test.method)) frm.mrppt$test.method else ddd$test.method
	eps=if(is.null(ddd$eps)) frm.mrppt$eps else ddd$eps
	
	frm.gradp=formals(grad.smoothp)
	bw=if(is.null(ddd$bw)) frm.gradp$bw else ddd$bw 
	kernel=if(is.null(ddd$kernel)) frm.gradp$kernel else ddd$kernel
	adjust=if(is.null(ddd$adjust)) frm.gradp$adjust else ddd$adjust
	
	ddd[['test.method']]=
	ddd[['eps']]=
	ddd[['bw']]=
	ddd[['kernel']]=
	ddd[['adjust']]=NULL
	
	if(length(ddd)>0L) warning('Unknown argument(s): ', names(ddd))
	
#	ddd.mrppt=list(...); ddd.mrppt$y=y; # ddd.mrppt$method='permutation'
#	ddd.mrppt.idx=names(ddd.mrppt)%in%names(formals(mrpp.test.mrpp))
#	ddd.mrppt=ddd.mrppt[ddd.mrppt.idx]
#	ddd.mrppt.method = if('test.method'%in%names(ddd.mrppt)) ddd.mrppt$test.method else NULL
	var.sign = integer(R)
	var.rank = numeric(R)

#	ddd.grad=list(...); 
#	ddd.grad$y=y
#	ddd.grad$permutedTrt=permutedTrt
#	ddd.grad$distObj=mrpp.obj$distObj
#	ddd.grad$mrpp.stats=numeric(0L)
	if(importance=='grad.energy') {
#		ddd.grad$bw=Inf
#		ddd.grad$adjust='weighted.mean'
		bw=Inf
		adjust='weighted.mean'
	}else if(is.null(adjust)) adjust='log scale'

    ans=vector('list'); 
	attr(ans, 'parameters')=list(
		importance=importance, inc.thresh=inc.thresh, exc.thresh=exc.thresh, size.inc=size.inc, stepwise=stepwise, niter=niter, 
		mrpp.test.options=list(test.method=test.method, eps=eps), 
		grad.smoothp.options=list(bw=bw, adjust=adjust, kernel=kernel))
	
#	ddd.grad.idx=names(ddd.grad)%in%names(formals(grad.smoothp))
#	ddd.grad=ddd.grad[ddd.grad.idx]
	
#print(proc.time())
    repeat{
        if(verbose && (i%%verbose==0L)) {cat('iteration',i-1L,'...')
                    time0=proc.time()[3L]}
        idx=idx[imptnc<=imptnc.threshold]
		
        if(length(idx)==0){
            ans[[i]]=list(iter=i-1L,`var.idx(sorted)`=integer(0L), `importance(sorted)`=numeric(0L), 
                        p.value=numeric(0L),deleted.p.value=ans[[1L]]$p.value, var.sign=var.sign, var.rank=var.rank)
            if(verbose){
				cat('\b\b\b:\t','no variables left; mrpp.p = ; deleted.mrpp.p =',ans[[i]]$deleted.p.value,
                        ';', proc.time()[3L]-time0,'seconds passed;\n')
				message('no variables left in the inclusion set')
			}
			attr(ans, 'status')=c(attr(ans, 'status'), ret.msg[3])
            return(returnBVS(ans))
        }
		
		y$data.env$y=y$data.env$y.bakBVS[,idx,drop=FALSE]
		y$R=length(idx); 
		y$distObj=dist0=distFunc(y$y)
		mrpp.stats0=.Call(mrppstats, y$distObj, y$permutedTrt, y$weight.trt,PACKAGE='MRPP')
        mrpp.rslt=mrpp.test.mrpp(y,test.method=test.method, eps=eps)
        imptnc=switch(importance,
			grad.smoothp=, 
			grad.energy=,
			approx.keep1 ={
				grad.smoothp.mrpp(y, bw=bw, kernel=kernel, adjust=adjust, mrpp.stats=mrpp.stats0, r=seq_len(y$R), test=FALSE)
#				ddd.grad$y=y[,idx,drop=FALSE]
#				ddd.grad$distObj=dist0
#				ddd.grad$mrpp.stats=mrpp.stats0
#				do.call('grad.smoothp', ddd.grad)
			},
            p.grad.dist =get.p.dd.dw(y$y,y$permutedTrt,...) # CHECKME
		)
		if(importance=='approx.keep1') imptnc=p.value(imptnc, type='keep1')
        var.ord=order(imptnc)
		
		
        ans[[i]]=list(iter=i-1L, `var.idx(sorted)`=idx[var.ord], `importance(sorted)`=imptnc[var.ord],
                      p.value=mrpp.rslt$p.value,
                      deleted.p.value=next.deleted.p#, var.sign=var.sign, var.rank=var.rank
					  )
        imptnc.threshold=if(stepwise) {tmp=max(imptnc); .5*(max(-Inf,imptnc[imptnc<tmp])+tmp)} else inc.thresh
#        idx=idx[imptnc<max(imptnc)] else idx=idx[imptnc<inc.thresh]
        if(exc.thresh>=0) {
            xcl=c(xcl, idx[imptnc>imptnc.threshold])
            dist.del=distFunc(y$data.env$y.bakBVS[,xcl,drop=FALSE])
            if(all(!is.na(dist.del)) && length(xcl)>0L)
#                ans[[i]]$deleted.p.value=mrpp.test.dist(dist.del, permutedTrt=permutedTrt)$p.value
                #next.deleted.p=mrpp.test.dist(dist.del, permutedTrt=permutedTrt)$p.value
				
#				ddd.mrppt$y$distObj=dist.del; ddd.mrppt$y$R=length(xcl)
#				next.deleted.p=do.call(mrpp.test.mrpp, ddd.mrppt)$p.value
				y$distObj=dist.del; y$R=length(xcl)
				next.deleted.p=mrpp.test.mrpp(y, test.method=test.method, eps=eps)$p.value
        }
		
		var.sign[idx]=-2L*(imptnc<=inc.thresh)+1L; var.sign[idx[imptnc==inc.thresh]]=0L
		var.rank[idx]=rank(imptnc, ties.method='average')
		ans[[i]]$var.sign=var.sign
		ans[[i]]$var.rank=var.rank
		
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
