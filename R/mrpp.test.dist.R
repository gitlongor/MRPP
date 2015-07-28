mrpp.weight.trt=function(weight.trt, trt)
{
    tabtrt=table(trt)
    ntrt=length(tabtrt)
	if(is.character(weight.trt)){
		wtmethod=match.arg(weight.trt, c('df','n','equal','terms'))
		weight.trt=switch(wtmethod, 
			'df'=(tabtrt-1)/(sum(tabtrt)-ntrt), 
			'n' =tabtrt/sum(tabtrt), 
			'equal' = 1/ntrt, 
			'terms' = local({ num = tabtrt*(tabtrt-1); ans = num/sum(num); ans[num==0]=0; ans})
		)
	}else wtmethod='custom'
	wtmethod=factor(wtmethod, levels=c('df','n','equal','terms','custom'))
	if(!is.double(weight.trt)) weight.trt = as.double(weight.trt)
	if(length(weight.trt)!=ntrt) weight.trt = rep(weight.trt, length=ntrt)
	names(weight.trt)=names(tabtrt)
	list(wtmethod=wtmethod, weight.trt=weight.trt)
}
mrpp.test.dist <-
function(y, trt, B=as.integer(min(nparts(table(trt)), 1e4L)), permutedTrt, weight.trt='df', eps=1e-8, ...) ## this uses C code
## y is a dist object; weight.trt: 0=sample size-1; 1=sample size
{
    if(missing(y) || !inherits(y,'dist')) stop('dist object missing or incorrect')
    N= as.integer( ( 1 + sqrt(1.0 + 8.0 * length(y))) * .5   +.5);          ### i.e.,   attr(y,'Size'), however, attributes might be lost during subsetting. 
    if(missing(trt)) {  ## recovering trt from the first permutation
      trt=trt.permutedTrt(permutedTrt)
    }
    if(missing(permutedTrt)) {
        permutedTrt=permuteTrt(trt,B, ...)
        dname=paste('"dist" object',deparse(substitute(y)), 
                             'and treatment group', deparse(substitute(trt)))
    }else dname=paste('"dist" object',deparse(substitute(y)), 
                             'and permuted treatment', deparse(substitute(permutedTrt)))
    B=nperms.permutedTrt(permutedTrt)

    tabtrt=table(trt)[names(permutedTrt)]
    ntrt=length(tabtrt)
	tmp=mrpp.weight.trt(weight.trt, as.factor(trt))
	wtmethod=tmp$wtmethod; weight.trt=tmp$weight.trt[names(permutedTrt)]
	
    stats=.Call(mrppstats,y,permutedTrt, as.numeric(weight.trt), PACKAGE='MRPP')

    ans=list(statistic=c("MRPP statistic"=stats[1L]), 
			 all.statistics=stats, 
             p.value=structure(mean(stats[1L]-stats>=-eps), midp = midp(stats, eps)),
			 parameter=structure(c("number of permutations"=B, 'weight method'=wtmethod), 'weight methods'=levels(wtmethod), weight.trt=weight.trt),
             data.name=dname, 
			 #  .Random.seed=attr(permutedTrt,'.Random.seed'),  ## this is commented out since the random number seed in gmp:::urand.bigz will be used. 
             method=sprintf('%d-sample MRPP test',ntrt)
             )
    class(ans)=c('mrpp.test','htest')
    ans
}

p.value.mrpp.test=function(x, type=c('raw','midp'),...)
{
	type=match.arg(type)
	switch(type, 
		raw = x$p.value[[1L]], 
		midp = attr(x$p.value, 'midp')
	)
}

mrpp.test.default <-
function(y, ...) {
    repl.text=paste("Response data ", deparse(substitute(y)), sep='')
	y=dist(y)
	.Class=c(.Class, 'dist')
    ans=NextMethod(.Generic, y)
    ans$data.name=gsub("\"dist\" object dist(y)", repl.text, ans$data.name, fixed=TRUE)
    ans
}

mrpp.test.formula <-
function(y, data, ...) 
{
    if (missing(y) || (length(y) != 3L) || (length(attr(terms(y[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    if(missing(data)) mf = model.frame(y) else mf <- model.frame(y, data=data)
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    res=mf[[response]]
    ans=mrpp.test(dist(res),trt=g, ...)
    ans$data.name=paste("'formula object:'", paste(names(mf), collapse = " ~ "))
    if(!missing(data)) ans$data.name = paste(ans$data.name, "in data.frame", substitute(data))
    ans
}

mrpp.test <-
function(y,...) UseMethod("mrpp.test")

