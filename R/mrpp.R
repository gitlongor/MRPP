
.mrpp.expr=expression({
	if(missing(trt)&&missing(permutedTrt))
		stop('"trt" and "permutedTrt" cannot be both missing')
    if(missing(trt)) {  ## recovering trt from the first permutation
      trt=trt.permutedTrt(permutedTrt)
    }
	permutedTrt.env=new.env(hash=FALSE, parent=constEnv)
    if(missing(permutedTrt)) {
		delayedAssign('permutedTrt', permuteTrt(trt,B, idxOnly), eval.env=environment(), assign.env=permutedTrt.env)
	}else{
		permutedTrt.env$permutedTrt = permutedTrt
		idxOnly = !is.na(attr(permutedTrt, 'idx')[1L])
	}
	delayedAssign('B', nperms.permutedTrt(permutedTrt), eval.env=permutedTrt.env, assign.env=permutedTrt.env)
		n=table(trt)
		cn=cumsum(n)
		ntrt=length(n)
		N=cn[ntrt]
		ordn=order(-n, names(n), decreasing=FALSE)
		trt=ordered(trt, levels=names(n)[ordn])
		trtc=levels(trt)
		tabtrt = n[ordn]
	tmp=mrpp.weight.trt(weight.trt, as.factor(trt))
	wtmethod=tmp$wtmethod; weight.trt=tmp$weight.trt[trtc]
})

mrpp <- function(y, trt, B=as.integer(min(nparts(table(trt)), 1e4L)), permutedTrt, weight.trt='df', distFunc=dist, idxOnly=FALSE)
{
	eval(.mrpp.expr)
    if(missing(permutedTrt)) {
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else{
		dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permuted treatment', deparse(substitute(permutedTrt)))
	}
if(FALSE){	
    if(missing(trt)) {  ## recovering trt from the first permutation
      trt=trt.permutedTrt(permutedTrt)
    }
	permutedTrt.env=new.env(hash=FALSE, parent=constEnv)
    if(missing(permutedTrt)) {
		delayedAssign('permutedTrt', permuteTrt(trt,B, idxOnly), eval.env=environment(), assign.env=permutedTrt.env)
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else{
		permutedTrt.env$permutedTrt = permutedTrt
		dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permuted treatment', deparse(substitute(permutedTrt)))
		idxOnly = !is.na(attr(permutedTrt, 'idx')[1L])
	}
	delayedAssign('B', nperms.permutedTrt(permutedTrt), eval.env=permutedTrt.env, assign.env=permutedTrt.env)

	{
		n=table(trt)
		cn=cumsum(n)
		ntrt=length(n)
		N=cn[ntrt]
		ordn=order(-n, names(n), decreasing=FALSE)
		trt=ordered(trt, levels=names(n)[ordn])
		trtc=levels(trt)
		tabtrt = n[ordn]
	}
	
	tmp=mrpp.weight.trt(weight.trt, as.factor(trt))
	wtmethod=tmp$wtmethod; weight.trt=tmp$weight.trt[trtc]
}
	if(N != NROW(y)) stop('NROW(y) != length(trt)')
	R = NCOL(y)
	structure(list(distObj=distFunc(y), 
					n=tabtrt, 
					ntrt=ntrt,
					nparts=nparts(tabtrt),
					nobs = N, 
					R=R, 
					B.requested=B, 
					trt=trt,
					permutedTrt.env=permutedTrt.env, 
					idxOnly = idxOnly,
					weight.trt=structure(weight.trt, 'method'=wtmethod), 
					distFunc = distFunc
			  ),
			  class='mrpp'
	)
}

