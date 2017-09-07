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
.pdfmethods=c('p3tlnormmix','pearson3','pearson3gca','tlnorm','gammagca','permutation')
.pdfnmoment=structure(c(4L, 3L, 4L, 3L, 4L, 0L), names=.pdfmethods)
mrpp.test.mrpp = function(y, test.method='pearson3gca', eps=1e-8, ... )
{
	if(!missing(...)).NotYetUsed(..., error=FALSE)
	dname.env=new.env(hash=FALSE, parent=globalenv())
	dname = delayedAssign('data.name', 
		paste('"mrpp" object',dQuote({
			tmp=deparse(substitute(y),width.cutoff=20L)
			if(length(tmp)>1L) paste(tmp[1L], '...')  else tmp[1L]
		}))
		,assign.env=dname.env)
	
	methods=unlist(strsplit(test.method[1L], ".", fixed=TRUE))
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
		cums=cumulant.mrpp(y, order=seq_len(nmoms))
		cums[2L]=sqrt(cums[2L])
		cums[-2:-1]=cums[-2:-1]/cums[2L]^(2+seq_len(nmoms-2L))
		# cums =c(mean, sd, skew, exkurt)
	}
	if(pdfmethod%in%c('pearson3','pearson3gca','tlnorm','p3tlnormmix','gammagca') && is.null(kernel)) {# no permutation
		tmpPermutedTrt=permuteTrt1(y$trt) # avoids evaluation fo permutedTrt
		stats=.Call(mrppstats,y$distObj,tmpPermutedTrt, as.numeric(y$weight.trt), PACKAGE='MRPP')
		B=1L
	}else{
		stats=.Call(mrppstats,y$distObj,y$permutedTrt, as.numeric(y$weight.trt), PACKAGE='MRPP')
		B=y[['permutedTrt.env']]$B
	}

	pval=switch(pdfmethod,
		pearson3=ppearson3(stats[1L], cums[1L], cums[2L], cums[3L]), 
		pearson3gca=ppearson3gca(stats[1L], cums[1L], cums[2L], cums[3L], cums[4L]),
		tlnorm=ptlnorm(stats[1L], cums[1L], cums[2L], cums[3L]), 
		p3tlnormmix=pMixP3Tln(stats[1L], cums[1L], cums[2L], cums[3L], cums[4L],proper=FALSE), 
		gammagca=pgammagca(stats[1L], cums[1L], cums[2L], cums[3L], cums[4L]), 
		permutation=,
		permutation.midp=mean(stats[1L]-stats>=-eps),
		NA_real_
	)
	if(pdfmethod=='permutation' || pdfmethod=='permutation.midp'){
		midp=local({
			adj=0.5/y$nparts
			if(pval > adj) as.numeric(pval - adj) else p.empirical(stats, eps)
		})
		pval=if(pdfmethod=='permutation') structure(pval, midp=midp) else structure(midp, raw=pval)
	}
	if(!is.null(kernel)){
		bw=bw.matchpdf(stats, kernel, pdf=switch(pdfmethod,
			pearson3=dpearson3(stats, cums[1L], cums[2L], cums[3L]), 
			pearson3gca=dpearson3gca(stats, cums[1L], cums[2L], cums[3L], cums[4L]), 
			p3tlnormmix=dMixP3Tln(stats, cums[1L], cums[2L], cums[3L], cums[4L],proper=FALSE), 
			tlnorm=dtlnorm(stats, cums[1L], cums[2L], cums[3L], cums[4L]),
			gammagca=dgammagca(stats, cums[1L], cums[2L], cums[3L], cums[4L])
			)
		)
		pval=pkde(stats, bw, kernel)(stats[1L])
	}
	pval.method.string=gsub(" +", " ", paste0(collapse=" ", c(
		sprintf("%d-sample MRPP test with", y$ntrt),
		switch(pdfmethod, 
			pearson3=, pearson3gca=, tlnorm=, p3tlnormmix=, gammagca=local({c(
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
             data.name=dname.env, 
             method=pval.method.string
	)
    class(ans)=c('mrpp.test','htest')
    ans	
}


`$.mrpp.test`=function(x, name)
{# proper extraction of data.name component
	x=unclass(x)
	#ans=getElement(x, substitute(name))
	ans=x[[substitute(name), exact=FALSE]]
	if(is.environment(ans)) ans$data.name else ans
}

.fix.mrpp.test.data.name=quote({
	#this.call=match.call()
	if(!is.null(this.call$trt)) {
		paste(tmpdname, 'and treatment factor', deparse(this.call$trt))
	}else if(!is.null(this.call$permutedTrt)){
		paste(tmpdname, 'and permuted treatment', deparse(this.call$permutedTrt))
	}
})
mrpp.test.default <- eval(bquote(
function(y, ...) {
	this.call=match.call()
	
	ddd=list(...)
	ddd$y=y
	idx=names(ddd)%in%names(formals(mrpp.matrix))
	mrpp.obj=do.call('mrpp', ddd[idx])
	ans=do.call('mrpp.test.mrpp', c(list(y=mrpp.obj), ddd[!idx]))
	delayedAssign('data.name', {
			tmpdp=deparse(substitute(y), width.cutoff=20L)
			tmpdname=paste("Response data", if(length(tmpdp)>1L) paste(tmpdp[1L],'...') else tmpdp)
			.(.fix.mrpp.test.data.name)
		}
		,assign.env=ans[['data.name']])
    ans
}))
attr(mrpp.test.default, 'srcref')=NULL

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
	delayedAssign('data.name', {
		tmpdname=paste("'formula'", deparse(y))
		if(!is.null(this.call$data)) {
			paste(tmpdname, "for dataset", sQuote({
				tmpdp=deparse(this.call$data, width.cutoff=20L)
				if(length(tmpdp)>1L) paste(tmpdp[1L],'...') else tmpdp})
			)
		}else tmpdname
	},assign.env=ans[['data.name']])
    ans
}

mrpp.test.dist <- eval(bquote(function(y, ...) 
{
#	dname=paste("'dist' object", deparse(substitute(y)))
#	eval(.fix.mrpp.test.data.name)
#	resp=suppressWarnings(cmdscale(y, attr(y, 'Size')-1L))
#	ddd=list(...); ddd$y=resp
#	ans=do.call('mrpp.test.mrpp', ddd)
#	ans$data.name=dname
#   ans
	this.call=match.call()
	this.call[[1L]]=as.symbol('mrpp')
	idx = c(1L, which(names(this.call)[-1L]%in%names(formals(mrpp.dist)))+1L)
	mrpp.obj=eval.parent(this.call[idx])
	ddd=list(...); ddd$y=mrpp.obj
	idx=names(ddd)%in%names(formals(mrpp.test.mrpp))
	ans=do.call('mrpp.test.mrpp', ddd[idx])
	delayedAssign('data.name',{
			tmpdname = paste("'dist' object", deparse(substitute(y)))
			.(.fix.mrpp.test.data.name)
		}, assign.env=ans[['data.name']])
	ans
}))
attr(mrpp.test.dist, 'srcref')=NULL

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
