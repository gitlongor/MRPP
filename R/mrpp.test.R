mrpp.weight.trt=function(weight.trt, trt)
{
    tabtrt=table(trt)
    ntrt=length(tabtrt)
	if(is.character(weight.trt)){
		wtmethod=match.arg(weight.trt, c('df','n','equal','pairs'))
		weight.trt=switch(wtmethod, 
			'df'=(tabtrt-1)/(sum(tabtrt)-ntrt), 
			'n' =tabtrt/sum(tabtrt), 
			'equal' = 1/ntrt, 
			'pairs' = local({ num = tabtrt*(tabtrt-1); ans = num/sum(num); ans[num==0]=0; ans})
		)
	}else wtmethod='custom'
	wtmethod=factor(wtmethod, levels=c('df','n','equal','pairs','custom'))
	if(!is.double(weight.trt)) weight.trt = as.double(weight.trt)
	if(length(weight.trt)!=ntrt) weight.trt = rep(weight.trt, length=ntrt)
	names(weight.trt)=names(tabtrt)
	list(wtmethod=wtmethod, weight.trt=weight.trt)
}
.pdfmethods=c('pearson3','pearson3gca','gammagca','permutation')
.pdfnmoment=structure(c(3L, 4L, 4L, 0L), names=.pdfmethods)
mrpp.test.mrpp = function(y, method='pearson3gca', eps=1e-8, ... )
{
	if(!missing(...)).NotYetUsed(..., error=FALSE)
	dname = paste('"mrpp" object',dQuote(substr(as.character(match.call()[['y']]), 1L, 20L)))
	
	methods=unlist(strsplit(method[1L], ".", fixed=TRUE))
	if(length(methods)>2L) stop('unknown method')
	pdfmethod=.pdfmethods[pmatch(methods[1L], .pdfmethods)]
	nmoms=.pdfnmoment[pdfmethod]
	if(pdfmethod=='permutation') {
		kernel=NULL
		if(length(methods)==2L) pdfmethod='permutation.midp'
	}else{
		if(length(methods)==2L) kernel=.kernels[pmatch(methods[2L], .kernels)] else kernel=NULL
	}
	if( nmoms > 0L){# need moments
		cums=cumulant(y, order=seq_len(nmoms))
		cums[-2:-1]=cums[-2:-1]/cums[2L]^((2+seq_len(nmoms-2L))/2)
		cums[2L]=sqrt(cums[2L])
		# cums =c(mean, sd, skew, exkurt)
	}
	if(pdfmethod%in%c('pearson3','pearson3gca','gammagca') && is.null(kernel)) {# no permutation
		tmpPermutedTrt=permuteTrt(y$trt, B=1L)
		stats=.Call(mrppstats,y$distObj,tmpPermutedTrt, as.numeric(y$weight.trt), PACKAGE='MRPP')
		B=1L
	}else{
		stats=.Call(mrppstats,y$distObj,y$permutedTrt.env$permutedTrt, as.numeric(y$weight.trt), PACKAGE='MRPP')
		B=y$permutedTrt.env$B
	}

	pval=switch(pdfmethod,
		pearson3=ppearson3(stats[1L], cums[1L], cums[2L], cums[3L]), 
		pearson3gca=ppearson3gca(stats[1L], cums[1L], cums[2L], cums[3L], cums[4L]), 
		gammagca=pgammagca(stats[1L], cums[1L], cums[2L], cums[3L], cums[4L]), 
		permutation=,
		permutation.midp=mean(stats[1L]-stats>=-eps),
		NA_real_
	)
	if(pdfmethod=='permutation' || pdfmethod=='permutation.midp'){
		midp=local({
			adj=0.5/y$nparts
			if(pval > adj) as.numeric(pval - adj) else midp.empirical(stats, eps)
		})
		pval=if(pdfmethod=='permutation') structure(pval, midp=midp) else structure(midp, raw=pval)
	}
	if(!is.null(kernel)){
		bw=bw.matchpdf(stats, kernel, pdf=switch(pdfmethod,
			pearson3=dpearson3(stats, cums[1L], cums[2L], cums[3L]), 
			pearson3gca=dpearson3gca(stats, cums[1L], cums[2L], cums[3L], cums[4L]), 
			gammagca=dgammagca(stats, cums[1L], cums[2L], cums[3L], cums[4L])
			)
		)
		pval=pkde(stats, bw, kernel)(stats[1L])
	}
	pval.method.string=gsub(" +", " ", paste0(collapse=" ", c(
		sprintf("%d-sample MRPP test with", y$ntrt),
		switch(pdfmethod, 
			pearson3=, pearson3gca=, gammagca=local({c(
			'approximate p-value method', 
			dQuote(if(is.null(kernel)) pdfmethod else paste0(c(pdfmethod, kernel), collapse="."))
			)}),
			permutation = if(B == y$nparts) 'exact p-value' else 'Monte Carlo p-value', 
			permutation.midp = "mid p-value"), 
		sprintf("and weighting method %s", dQuote(attr(y$weight.trt, 'method')))
		)))
	
	
    ans=list(statistic=c("MRPP statistic"=stats[1L]), 
			 all.statistics=stats, 
             p.value=pval,
			 parameter=structure(
				c("number of permutations"=B), 
				weight.trt=y$weight.trt, 
				pdfmethod=pdfmethod, kernel=as.character(kernel), exact= B==y$nparts && pdfmethod=="permutation"
			 ),
             data.name=dname, 
             method=pval.method.string
	)
    class(ans)=c('mrpp.test','htest')
    ans	
}

.fix.mrpp.test.data.name=expression({
	this.call=match.call()
	if(!is.null(this.call$trt)) {
		dname=paste(dname, 'and treatment factor', deparse(this.call$trt))
	}else if(!is.null(this.call$permutedTrt)){
		dname=paste(dname, 'and permuted treatment', deparse(this.call$permutedTrt))
	}
})
mrpp.test.default <-
function(y, ...) {
    dname=paste("Response data", deparse(substitute(y)))
	eval(.fix.mrpp.test.data.name)
	ddd=list(...)
	ddd$y=y
	idx=names(ddd)%in%names(formals(mrpp))
	y=do.call('mrpp', ddd[idx])
	ans=do.call('mrpp.test.mrpp', c(list(y=y), ddd[!idx]))
	ans$data.name = dname
    ans
}

mrpp.test.formula <-
function(y, ...) 
{
    if ((length(y) != 3L) || (length(attr(terms(y[-2L]), "term.labels")) != 1L)) 
        stop("'formula' not yet supported")
	this.call=match.call()
	ddd=list(...)
	idx=names(ddd)%in%names(formals(model.frame.default))
	mf.list=ddd[idx]; mf.list$formula=y
    mf = do.call('model.frame.default', mf.list)
	mrpp.list=ddd[!idx]
	
    resp.col <- attr(attr(mf, 'terms'), "response")
	mrpp.list$y=mf[[ resp.col ]]
	mrpp.list$trt=mf[[ -resp.col ]]

    ans=do.call('mrpp.test', mrpp.list)
    ans$data.name=paste("'formula'", deparse(y))
    if(!is.null(this.call$data)) ans$data.name = paste(ans$data.name, "for dataset", sQuote(deparse(this.call$data)))
    ans
}

mrpp.test.dist <- function(y, ...) 
{
	dname=paste("'dist' object", deparse(substitute(y)))
	eval(.fix.mrpp.test.data.name)
	resp=suppressWarnings(cmdscale(y, attr(y, 'Size')-1L))
	ddd=list(...); ddd$y=resp
	ans=do.call('mrpp.test.default', ddd)
	ans$data.name=dname
    ans
}

mrpp.test <-
function(y,...) UseMethod("mrpp.test")

p.value.mrpp.test=function(x, type=c('raw','midp'),...)
{
	if(!missing(...)).NotYetUsed(..., error=FALSE)
	type=match.arg(type)
	attrs=attributes(x$p.value)
	midraw=if(is.null(attrs$midp)) {
				if(is.null(attrs$raw)){
						rep(x$p.value, 2L)
				}else c(x$p.value, attrs$raw) 
			}else c(attrs$midp, x$p.value)
	switch(type, 
		raw = midraw[2L],
		midp = midraw[1L]
	)
}
