
.mrpp.expr=quote({
	if(missing(trt)&&missing(permutedTrt))
		stop('"trt" and "permutedTrt" cannot be both missing')
    if(missing(trt)) {  ## recovering trt from the first permutation
      trt=trt.permutedTrt(permutedTrt)
    }
	if(missing(B)){
		B=if(missing(permutedTrt)) {
			as.integer(min(nparts(table(trt)), 1e4L)) 
		}else nperms.permutedTrt(permutedTrt)
	}
	permutedTrt.env=new.env(hash=FALSE, parent=constEnv)
    if(missing(permutedTrt)) {
		delayedAssign('permutedTrt', permuteTrt(trt,B), eval.env=environment(), assign.env=permutedTrt.env)
	}else{
		permutedTrt.env$permutedTrt = permutedTrt
		#idxOnly = !is.na(attr(permutedTrt, 'idx')[1L])
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

mrpp <- function(y, trt, B, permutedTrt, weight.trt='df',...) UseMethod('mrpp')

mrpp.matrix <- eval(bquote(function(y, trt, B, permutedTrt, weight.trt='df', distFunc=dist,...)
{	stopifnot(is.numeric(y))
	.(.mrpp.expr)
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
					#idxOnly = idxOnly,
					weight.trt=structure(weight.trt, 'method'=wtmethod), 
					distFunc = distFunc
			  ),
			  class='mrpp'
	)
}, # of function
) # of bquote
) # of eval 

mrpp.dist=eval(bquote(function(y, trt, B, permutedTrt, weight.trt='df', distFunc=dist, R, ...)
{
	.(.mrpp.expr)
	if(N != attr(y,'Size')) stop('attr(y,"Size") != length(trt)')
	if(missing(R)) R = attr(y,'Size')-1L
	structure(list(distObj=y, 
					n=tabtrt, 
					ntrt=ntrt,
					nparts=nparts(tabtrt),
					nobs = N, 
					R=R, 
					B.requested=B, 
					trt=trt,
					permutedTrt.env=permutedTrt.env, 
					#idxOnly = idxOnly,
					weight.trt=structure(weight.trt, 'method'=wtmethod), 
					distFunc = distFunc
			  ),
			  class='mrpp'
	)
} # of function
) # of bquote
) # of eval
