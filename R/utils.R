safeseq =
function (from = 1L, to = 1L, by = 1L, ...) 
{
    disc = by * (from - to)
    if (disc > 0) {
        vector(class(disc), 0L)
    }
    else seq(from = from, to = to, by = by, ...)
}


midp = function(x, eps=1e-8)
{
	mean(x<x[1L]-eps) + .5* mean(x>=x[1L]-eps & x<=x[1L]+eps)
}


factorial.rising=function(start, nterms)
{
	stopifnot(nterms>0)
	if(start!=as.integer(start)).NotYetImplemented()
	ans = prod.bigz(as.bigz(seq(from=start, length=nterms)))
	tryCatch(as.integer(ans), warning=function(w)ans)
}

`[[.dist`=function(x, i,j)
{
	if(i==j) return(0)
	if(i>j){tmp=i; i=j; j=tmp}
	unclass(x)[attr(x,'Size')*(i-1L) - i*(i-1L)/2L + j-i]
}
